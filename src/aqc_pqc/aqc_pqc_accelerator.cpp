#define _USE_MATH_DEFINES
#include <cmath>

#include "aqc_pqc/AqcPqcAccelerator.h"
#include "logger.h"
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <nlopt.h>

namespace FastVQA{

const long double PI =   3.14159265358979323846264338327950288419716939937510L;
const long double PI_2 = 1.57079632679489661923132169163975144209858469968755L;

AqcPqcAccelerator::AqcPqcAccelerator(AqcPqcAcceleratorOptions options){

	if(options.accelerator_type != "quest")
		throw_runtime_error("No other accelerator than QuEST implemented");

	if(options.printGroundStateOverlap && options.initialGroundState == InitialGroundState::None)
		throw_runtime_error("printGroundStateOverlap=true but initial ground state not specified");

	logi("Using " + options.ansatz_name + " ansatz");

	this->options = options;
	this->env = createQuESTEnv();

}

AqcPqcAccelerator::~AqcPqcAccelerator(){

	if(initialized){
		this->finalize();
		initialized = false;
	}
}

void AqcPqcAccelerator::initialize(PauliHamiltonian* h0, PauliHamiltonian* h1){

	if(initialized)
		throw_runtime_error("AqcPqcAccelerator has already been initialized!");

	this->nbQubits = h0->nbQubits;
	if(this->nbQubits != h1->nbQubits)
		throw_runtime_error("Dimensions of initial and final Hamiltonian do not match");

	hamil_int.nbQubits = nbQubits;
	hamil_int.coeffs0 = h0->coeffs;
	hamil_int.pauliOpts = h0->pauliOpts;
	hamil_int.initial_type = h0->type;

	std::vector<std::string> sequences;

	for(std::vector<double>::size_type i = 0; i < hamil_int.coeffs0.size(); ++i){
		std::string s;

		for(int j = 0; j < nbQubits; ++j)
			s+=std::to_string(h0->pauliOpts[nbQubits*i+j]);

		std::vector<std::string>::iterator itr = std::find(sequences.begin(), sequences.end(), s);
		if (itr != sequences.cend())
			throw_runtime_error("Hamiltonian contains multiple equivalent terms");
		sequences.push_back(s);
		hamil_int.coeffs1.push_back(0);
	}

	for(std::vector<double>::size_type i = 0; i < h1->coeffs.size(); ++i){
		std::string s;

		for(int j = 0; j < nbQubits; ++j)
			s+=std::to_string(h1->pauliOpts[nbQubits*i+j]);

		std::vector<std::string>::iterator itr = std::find(sequences.begin(), sequences.end(), s);
		if (itr != sequences.cend()){
			int index = std::distance(sequences.begin(), itr);
			hamil_int.coeffs1[index]+=h1->coeffs[i];
		}
		else{
			sequences.push_back(s);
			for(auto &c:s)
				hamil_int.pauliOpts.push_back(c-'0');
			hamil_int.coeffs0.push_back(0);
			hamil_int.coeffs1.push_back(h1->coeffs[i]);
		}
	}
	if(hamil_int.coeffs0.size() != hamil_int.coeffs1.size())
		throw_runtime_error("Problem with mixer hamiltonian calculation");

	unsigned long int keys[1];
	keys[0] = 1997;
	seedQuEST(&env, keys, 1);
	logd("Setting seed to " + std::to_string(keys[0]), options.log_level);
	this->qureg = createQureg(this->nbQubits, env);
	workspace = createQureg(this->nbQubits, env);


	if(options.outputLogToFile){
		logFile.open(options.logFileName);
		if(!logFile.is_open())
			throw_runtime_error("Error opening log file " + options.logFileName);

		logFile << std::fixed << std::setprecision(10);

	}

	initialized = true;

}

void AqcPqcAccelerator::finalize(){
	destroyQureg(this->qureg, env);
	destroyQureg(this->workspace, env);
	logFile.close();
}

typedef struct {
	Eigen::Vector<qreal, Eigen::Dynamic> Xi;
    Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> N;
} OptData;

typedef struct {
	std::vector<std::shared_ptr<Parameter>> *parameters;
	AqcPqcAccelerator *acc;
	PauliHamiltonian *h;
	OptData *optData;
} ConstrData;

template <class MatrixT>
bool isPsd(const MatrixT& A) {
  if (!A.isApprox(A.transpose())) {
    return false;
  }
  const auto ldlt = A.template selfadjointView<Eigen::Upper>().ldlt();
  if (ldlt.info() == Eigen::NumericalIssue || !ldlt.isPositive()) {
    return false;
  }
  return true;
}

double lin_system_f(unsigned n, const double *z, double *grad, void *data){

	OptData *d = (OptData *) data;

    if (grad) {
    	//UNIMPLEMENTED
        /*grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);*/
    }

    qreal ret=0;
    Eigen::Vector<qreal, Eigen::Dynamic> z_vect(n);
    for(unsigned int i = 0; i < n; ++i)
    	z_vect(i)=z[i];
    Eigen::Vector<qreal, Eigen::Dynamic> x = d->Xi+d->N*z_vect;
    for(int i = 0; i < d->Xi.rows(); ++i)
    	ret += x[i]*x[i];
    return ret;
}

double eq_constraint(unsigned n, const double *x, double *grad, void *data){
    //my_constraint_data *d = (my_constraint_data *) data;
    //double a = d->a, b = d->b;
    if (grad) {
    	//UNIMPLEMENTED
    }

    ConstrData *d = (ConstrData *) data;
    Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> H(d->parameters->size(), d->parameters->size());

    Eigen::Vector<qreal, Eigen::Dynamic> eps_vect(n);
	for(unsigned int i = 0; i < n; ++i)
		eps_vect(i)=x[i];
	Eigen::Vector<qreal, Eigen::Dynamic> res_eps = d->optData->Xi+d->optData->N*eps_vect;

	for(unsigned int i = 0; i < d->parameters->size(); ++i){
		(*d->parameters)[i]->value += res_eps[i];
		qreal original_i = (*d->parameters)[i]->value;

		for(unsigned int j = 0; j <= i; ++j){

			qreal original_j = (*d->parameters)[j]->value;

			(*d->parameters)[i]->value += PI_2;
			(*d->parameters)[j]->value += PI_2;
			qreal a = d->acc->_calc_expectation(d->h);

			(*d->parameters)[j]->value -= PI;
			qreal b = d->acc->_calc_expectation(d->h);

			(*d->parameters)[i]->value -= PI;
			(*d->parameters)[j]->value += PI;
			qreal c = d->acc->_calc_expectation(d->h);

			(*d->parameters)[j]->value -= PI;
			qreal dd = d->acc->_calc_expectation(d->h);

			(*d->parameters)[i]->value = original_i;
			(*d->parameters)[j]->value = original_j;

			H(i,j) = 0.25 * (a-b-c+dd);

			if(i != j)
				H(j,i) = H(i,j);
		}
	}

	for(unsigned int i = 0; i < d->parameters->size(); ++i){
		(*d->parameters)[i]->value -= res_eps[i];
	}

	return isPsd(H)?0:1000;
}


void AqcPqcAccelerator::run(){

	int nbSteps = options.nbSteps;
	this->ansatz = getAnsatz(options.ansatz_name, hamil_int.nbQubits);
	initOptimalParamsForMinusSigmaXHamiltonian(&ansatz);

	std::vector<std::shared_ptr<Parameter>> parameters = ansatz.circuit.getParamsPtrs();

	logi(std::to_string(nbSteps) + " steps");

	//double *first_order_terms = (double*) malloc(parameters.size() * sizeof(double));
	//double *second_order_terms[parameters.size()];
	//for(int i = 0; i < parameters.size(); i++)
	//	second_order_terms[i] = (double*)malloc((i+1) * sizeof(double));

	Eigen::Vector<qreal, Eigen::Dynamic> minus_q(parameters.size());
	Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> A(parameters.size(), parameters.size());
	nlopt_opt opt/*, lopt*/;

	for(int k = 1; k < nbSteps+1; ++k){

		/*for(int i = 0; i < parameters.size(); ++i)
			std::cerr<<parameters[i]->value<<" ";
		std::cerr<<"\n";*/

		double theta = (double)(k)/(nbSteps);

		PauliHamiltonian h = this->_calc_intermediate_hamiltonian(theta);
		logd("theta="+std::to_string(theta), options.log_level);
		logd(h.getPauliHamiltonianString(2), options.log_level);

		for(std::vector<std::shared_ptr<FastVQA::Parameter>>::size_type i = 0; i < parameters.size(); ++i){

			qreal original_i = parameters[i]->value;

			parameters[i]->value += PI_2;
			qreal a = _calc_expectation(&h);
			parameters[i]->value -= PI;
			qreal b = _calc_expectation(&h);
			parameters[i]->value = original_i;

			minus_q(i)=-0.5*(a-b);

			//if(fabsl(minus_q(i)) > 6e-17)
			//	minus_q(i)=0;


			parameters[i]->value = original_i;

			for(unsigned int j = 0; j <= i; ++j){

				qreal original_j = parameters[j]->value;

				/*if(i == j && i == 1){
					std::cerr<<original_i<<" "<<original_j<<"\n";
				}else
					continue;*/

				parameters[i]->value += PI_2;
				parameters[j]->value += PI_2;

				//for(int p = 0; p < parameters.size(); ++p)
				//	std::cerr<<parameters[p]->value<<"\n";

				qreal a = _calc_expectation(&h);

				parameters[j]->value -= PI;
				qreal b = _calc_expectation(&h);

				parameters[i]->value -= PI;
				parameters[j]->value += PI;
				qreal c = _calc_expectation(&h);

				parameters[j]->value -= PI;
				qreal d = _calc_expectation(&h);

				parameters[i]->value = original_i;
				parameters[j]->value = original_j;

				A(i,j) = 0.25 * (a-b-c+d);
				if(i == j && i == 1){
					std::cerr<<a<<" "<<b<<" "<<c<<" "<<d<<"\n";
				}

				if(i != j)
					A(j,i) = A(i,j);
			}
		}

		/*std::cerr<<-minus_q<<std::endl;
		std::cerr<<A<<std::endl;throw;*/

		Eigen::Vector<qreal, Eigen::Dynamic> Xi(parameters.size());
		Xi = A.fullPivHouseholderQr().solve(minus_q);

		//std::cerr<<"xi"<<Xi<<std::endl<<std::endl;

		Eigen::FullPivLU<Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic>> lu(A);
		Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> A_null_space = lu.kernel();
		//for faster way see https://stackoverflow.com/questions/34662940/how-to-compute-basis-of-nullspace-with-eigen-library

		int opt_dim = A_null_space.cols();

		double *lb, *ub;
		OptData data {Xi, A_null_space};

		if(options.checkHessian){
			opt = nlopt_create(NLOPT_LN_COBYLA, opt_dim);
			//opt = nlopt_create(NLOPT_LN_AUGLAG_EQ, opt_dim);
			//lopt = nlopt_create(NLOPT_LN_COBYLA, opt_dim);
			//nlopt_set_local_optimizer(opt, lopt);
			ConstrData constr_data {&parameters, this, &h, &data};
			nlopt_add_equality_constraint(opt, eq_constraint, &constr_data, 10e-1);

			lb = (double*) malloc(opt_dim * sizeof(double));
			ub = (double*) malloc(opt_dim * sizeof(double));
			for(int i = 0; i < opt_dim; ++i){
				lb[i] = -HUGE_VAL;ub[i]=-lb[i];
			}
			nlopt_set_lower_bounds(opt, lb);
			nlopt_set_upper_bounds(opt, ub);
			nlopt_set_min_objective(opt, lin_system_f, &data);
			nlopt_set_xtol_rel(opt, 1e-4);
			nlopt_set_xtol_abs1(opt, 1e-4);

		}else{
			opt = nlopt_create(NLOPT_LN_COBYLA, opt_dim);
			lb = (double*) malloc(opt_dim * sizeof(double));
			ub = (double*) malloc(opt_dim * sizeof(double));
			for(int i = 0; i < opt_dim; ++i){
				lb[i] = -HUGE_VAL;ub[i]=-lb[i];
			}
			nlopt_set_lower_bounds(opt, lb);
			nlopt_set_upper_bounds(opt, ub);
			nlopt_set_min_objective(opt, lin_system_f, &data);
			nlopt_set_xtol_rel(opt, 1e-4);
			nlopt_set_xtol_abs1(opt, 1e-4);
		}

		double *eps = (double*) malloc(opt_dim * sizeof(double));
		for(int i = 0; i < opt_dim; ++i)
			eps[i] = 0.5;
		double minf;

		int opt_res = nlopt_optimize(opt, eps, &minf);
		if (opt_res < 0) {
		    std::cerr<<"nlopt failed! error code: "<< opt_res <<std::endl;
		}
		else {
			 /*std::cerr<<"found minimum at f("<<eps[0];
			 for(int i = 1; i < opt_dim; ++i)
				 std::cerr<<","<<eps[i];
			 std::cerr<< ")= "<<minf<<"\n";*/

			Eigen::Vector<qreal, Eigen::Dynamic> eps_vect(opt_dim);
			for(int i = 0; i < opt_dim; ++i)
				eps_vect(i)=eps[i];


			//std::cerr<<eps_vect;
			//std::cerr<<"A_null:"<<A_null_space;
			Eigen::Vector<qreal, Eigen::Dynamic> res_eps = Xi+A_null_space*eps_vect;

			//std::cerr<<res_eps;//throw;

			for(unsigned int i = 0; i < parameters.size(); ++i){
				parameters[i]->value += res_eps[i];
			}
			//std::cerr<<"x"<<A*res_eps-minus_q<<std::endl;
		}



		qreal expectation = _calc_expectation(&h);
		//logi("Energy: " + std::to_string(expectation));

		free(eps);
		free(lb);free(ub);
		nlopt_destroy(opt);//nlopt_destroy(lopt);
		for(int i = 0; i < qureg.numAmpsTotal; ++i)
			std::cerr<<qureg.stateVec.real[i]*qureg.stateVec.real[i]+qureg.stateVec.imag[i]*qureg.stateVec.imag[i]<<" ";
		std::cerr<<"\n";
		if(options.compareWithClassicalEigenSolver){

			double id_term;
			Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> m = getIntermediateMatrixRepresentation(&h, &id_term);
			Eigen::Vector<std::complex<qreal>, Eigen::Dynamic> evals_vect = m.eigenvalues();

			std::vector<double> evals(evals_vect.size());
			for(int i = 0; i < evals_vect.size(); ++i){
				std::complex c = evals_vect.coeff(i);
				if(c.imag() != 0)
					loge("Imaginary part of eigenvalue is non-zero: " + std::to_string(c.imag()));
				evals[i] = c.real();
			}

			std::sort(evals.begin(), evals.end());

			std::string str_groundStateOverlap = "";
			std::string str_finalGroundStateOverlap = "";

			double gsOverlap, fgsOverlap;

			if(options.printGroundStateOverlap){

				if(options.solution == -1)
					throw_runtime_error("printGroundStateOverlap=true but final ground state is not set.");

				switch(options.initialGroundState){
				case PlusState:{
					double overlap_val_real = 0;
					double overlap_val_imag = 0;
					double p = (1-theta)*(1/sqrt(qureg.numAmpsTotal));
					for(long long int i = 0; i < qureg.numAmpsTotal; ++i){

						//std::cout << qureg.stateVec.real[i] << "+" << qureg.stateVec.imag[i] << "i  ";

						if(i == options.solution){
							overlap_val_real += (p + theta) * qureg.stateVec.real[i];
							overlap_val_imag += (p + theta) * qureg.stateVec.imag[i];
						}else{
							overlap_val_real += p * qureg.stateVec.real[i];
							overlap_val_imag += p * qureg.stateVec.imag[i];
						}
					}
					gsOverlap = pow(overlap_val_real,2)+pow(overlap_val_imag, 2);
					fgsOverlap = pow(qureg.stateVec.real[options.solution],2)+pow(qureg.stateVec.imag[options.solution], 2);
					str_groundStateOverlap = " GS overlap= " + std::to_string(gsOverlap);
					str_finalGroundStateOverlap = " GS overlap= " + std::to_string(fgsOverlap);
					break;}
				case None:
				default:
					throw_runtime_error("Not implemented. Err code: 11");
					break;
				}
			}

			std::ostringstream oss;
			oss << std::fixed;
			oss << std::setprecision(2);
			oss << ((int)(theta*10000))/100.;

			logi(oss.str()+"% "+"Exact ground state: " + std::to_string(evals[0]) + " calculated: " + std::to_string(expectation) + "    (diff="+std::to_string(std::abs(evals[0]-expectation))+")  " + str_groundStateOverlap + " " + str_finalGroundStateOverlap);

			if(options.outputLogToFile){
				logFile << evals[0] << " " << expectation << " " << (options.printGroundStateOverlap ? std::to_string(gsOverlap) + " " : "") << (options.printGroundStateOverlap ? std::to_string(fgsOverlap) : "") << std::endl;
			}

			int i_max = 0;double maxd = 0;

			for(int i = 0; i < qureg.numAmpsTotal;++i)
				if(pow(qureg.stateVec.real[i],2)+pow(qureg.stateVec.imag[i], 2) > maxd){
					maxd=pow(qureg.stateVec.real[i],2)+pow(qureg.stateVec.imag[i], 2);
					i_max = i;
				}
			std::cerr<< /*pow(qureg.stateVec.real[i],2)+pow(qureg.stateVec.real[i], 2)*/i_max <<" ";

		}


	}
	/*free(first_order_terms);
    for (int i = 0; i < parameters.size(); i++)
    	free(second_order_terms[i]);*/

}

PauliHamiltonian AqcPqcAccelerator::_calc_intermediate_hamiltonian(double lambda){

	std::vector<qreal> coeffs;

	for(std::vector<double>::size_type i = 0; i < hamil_int.coeffs0.size(); ++i){
		coeffs.push_back((1-lambda)*hamil_int.coeffs0[i]+lambda*hamil_int.coeffs1[i]);
	}

	PauliHamiltonian res(nbQubits, coeffs, hamil_int.pauliOpts);

	return res;
}

qreal AqcPqcAccelerator::calc_intermediate_expectation(ExperimentBuffer* buffer, double lambda, bool init_zero_state){

	if(lambda < 0 || lambda > 1)
		throw_runtime_error("Invalid lambda value");

	PauliHamiltonian h = _calc_intermediate_hamiltonian(lambda);
	return _calc_expectation(&h);

}

qreal AqcPqcAccelerator::_calc_expectation(PauliHamiltonian *h){
	initZeroState(qureg);
	run_circuit(ansatz.circuit, false);
	PauliHamil hamil;
	h->toQuestPauliHamil(&hamil);
	return calcExpecPauliHamil(qureg, hamil, workspace);
}

Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> AqcPqcAccelerator::getIntermediateMatrixRepresentation(PauliHamiltonian* h, double* id_coeff){

	if(hamil_int.initial_type == PauliHamiltonianType::General){
		throw_runtime_error("getIntermediateMatrixRepresentation not implemented for general h0 type");
	}else if(hamil_int.initial_type == PauliHamiltonianType::MinusSigmaX){

	Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> m = Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic>::Zero((1LL<<h->nbQubits), (1LL<<h->nbQubits));
	DiagonalOp diagOp = createDiagonalOp(h->nbQubits, env, true);

	PauliHamil h_cpy;
	h_cpy.numQubits = h->nbQubits;
	std::vector<qreal> coeffs;
	std::vector<int> codes;

	qreal x_coeff=0;
	*id_coeff=0;

	for(std::vector<double>::size_type i = 0; i < h->coeffs.size(); ++i){

		if(h->coeffs[i] == 0)
			continue;

		bool all_x=true, all_i=true;
		bool x_coeff_seen=false;bool id_coeff_seen=false;

		for(int q = 0; q < h->nbQubits; ++q){

			int opt = h->pauliOpts[i*h->nbQubits + q];
			if(opt != 1){
				all_x = false;
			}
			if(opt != 0){
				all_i = false;
			}
		}

		if(!all_x && !all_i){
			coeffs.push_back(h->coeffs[i]);
			for(int q = 0; q < h->nbQubits; ++q)
				codes.push_back(h->pauliOpts[i*h->nbQubits + q]);
		}else if(!all_i){
			if(x_coeff_seen)
				throw_runtime_error("Two all x terms");
			x_coeff = h->coeffs[i];
			x_coeff_seen=true;
		}else{
			if(id_coeff_seen)
				throw_runtime_error("Two id coeffs terms");
			*id_coeff=h->coeffs[i];
			id_coeff_seen = true;
		}
	}

	h_cpy.numSumTerms = coeffs.size();
	h_cpy.termCoeffs = &coeffs[0]; //conversion to c array
	h_cpy.pauliCodes = (enum pauliOpType*)(&codes[0]);

	if(h_cpy.numSumTerms > 0)
		initDiagonalOpFromPauliHamil(diagOp, h_cpy);

	//logw("Not distributed", options.log_level);
	for(long long i = 0; i < m.rows(); ++i){

		if(h_cpy.numSumTerms > 0)
			m.diagonal()[i] = diagOp.real[i];
		m(i, m.rows()-1-i) += x_coeff;
	}
	destroyDiagonalOp(diagOp, env);
	return m;
	}else if(hamil_int.initial_type == PauliHamiltonianType::SumMinusSigmaX){
		Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> m = Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic>::Zero(/*(1LL<<h->nbQubits)*/4, 4/*(1LL<<h->nbQubits)*/);
		std::cerr<<m<<std::endl;
		std::function<void(int, int, int)> f;
		f= [&m, &f](int q, int a, int b)->void {
			if(q == 1){
				m(a+1, b) = -1;
				m(a, b+1) = -1;
			}
			else{
				f(q-1, a, b);
				f(q-1, a, b+(1LL<<(q-1)));
				f(q-1, a+(1LL<<(q-1)), b);
				f(q-1, a+(1LL<<(q-1)), b+(1LL<<(q-1)));
			};
		};

		f(/*h->nbQubits-1*/2, 0, 0);
		/*f(h->nbQubits-1, 0, 0);
		f(h->nbQubits-1, 0, (1LL<<(h->nbQubits-1)));
		f(h->nbQubits-1, (1LL<<(h->nbQubits-1)), 0);
		f(h->nbQubits-1, (1LL<<(h->nbQubits-1)), (1LL<<(h->nbQubits-1)));*/
		std::cerr<<m<<"\n";
		return m;
	}
	throw_runtime_error("Not implemented");
}

}
