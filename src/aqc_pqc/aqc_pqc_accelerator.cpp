#define _USE_MATH_DEFINES
#include <cmath>

#include "aqc_pqc/AqcPqcAccelerator.h"
#include "logger.h"
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <nlopt.h>

namespace FastVQA{

AqcPqcAccelerator::AqcPqcAccelerator(AqcPqcAcceleratorOptions options){

	if(options.accelerator_type != "quest"){
		loge("No other accelerator than QuEST implemented");
		throw;
	}

	logi("Using " + options.ansatz_name + " ansatz");

	this->options = options;
	this->env = createQuESTEnv();

}

void AqcPqcAccelerator::initialize(PauliHamiltonian* h0, PauliHamiltonian* h1){

	this->nbQubits = h0->nbQubits;
	if(this->nbQubits != h1->nbQubits)
		throw_runtime_error("Dimensions of initial and final Hamiltonian do not match");

	hamil_int.nbQubits = nbQubits;
	hamil_int.coeffs0 = h0->coeffs;
	hamil_int.pauliOpts = h0->pauliOpts;

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

}

typedef struct {
	Eigen::VectorXd Xi;
    Eigen::MatrixXd N;
} OptData;

double lin_system_f(unsigned n, const double *z, double *grad, void *data){

	OptData *d = (OptData *) data;

    if (grad) {
    	//UNIMPLEMENTED
        /*grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);*/
    }

    double ret=0;
    Eigen::VectorXd z_vect(n);
    for(unsigned int i = 0; i < n; ++i)
    	z_vect(i)=z[i];
    Eigen::VectorXd x = d->Xi+d->N*z_vect;
    for(int i = 0; i < d->Xi.rows(); ++i)
    	ret += x[i]*x[i];
    return ret;//sqrt(x[1]);
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

	Eigen::VectorXd minus_q(parameters.size());
	Eigen::MatrixXd A(parameters.size(), parameters.size());
	nlopt_opt opt;

	for(int k = 0; k < nbSteps; ++k){

		PauliHamiltonian h = this->_calc_intermediate_hamiltonian((double)k/(nbSteps-1));

		for(std::vector<std::shared_ptr<FastVQA::Parameter>>::size_type i = 0; i < parameters.size(); ++i){

			double original_i = parameters[i]->value;
			std::cerr<<original_i<<" ";

			parameters[i]->value += M_PI_2;
			double a = _calc_expectation(&h);
			parameters[i]->value -= M_PI;
			double b = _calc_expectation(&h);
			parameters[i]->value = original_i;

			//std::cerr<<a<<" "<<b<<"\n";
			minus_q(i)=-0.5*(a-b);

			for(unsigned int j = 0; j <= i; ++j){

				double original_j = parameters[j]->value;

				parameters[i]->value += M_PI_2;
				parameters[j]->value += M_PI_2;
				double a = _calc_expectation(&h);

				parameters[j]->value -= M_PI;
				double b = _calc_expectation(&h);

				parameters[i]->value -= M_PI;
				parameters[j]->value += M_PI;
				double c = _calc_expectation(&h);

				parameters[j]->value -= M_PI;
				double d = _calc_expectation(&h);

				parameters[i]->value = original_i;
				parameters[j]->value = original_j;

				A(i,j) = 0.25 * (a-b-c+d);

				if(i != j)
					A(j,i) = A(i,j);
			}
		}std::cerr<<"\n";

		Eigen::VectorXd Xi(parameters.size());
		Xi = A.fullPivHouseholderQr().solve(minus_q);

		Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
		Eigen::MatrixXd A_null_space = lu.kernel();
		//for faster way see https://stackoverflow.com/questions/34662940/how-to-compute-basis-of-nullspace-with-eigen-library

		int opt_dim = A_null_space.cols();

		opt = nlopt_create(NLOPT_LN_COBYLA, opt_dim);
		OptData data {Xi, A_null_space};
		double *lb = (double*) malloc(opt_dim * sizeof(double));
		double *ub = (double*) malloc(opt_dim * sizeof(double));
		for(int i = 0; i < opt_dim; ++i){
			lb[i] = -HUGE_VAL;ub[i]=-lb[i];
		}
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
		nlopt_set_min_objective(opt, lin_system_f, &data);
		nlopt_set_xtol_rel(opt, 1e-4);
		nlopt_set_xtol_abs1(opt, 1e-4);

		double *eps = (double*) malloc(opt_dim * sizeof(double));
		for(int i = 0; i < opt_dim; ++i)
			eps[i] = 0.5;
		double minf;

		int opt_res = nlopt_optimize(opt, eps, &minf);
		if (opt_res < 0) {
		    std::cerr<<"nlopt failed! error code: "<< opt_res <<std::endl;
		}
		else {
			 std::cerr<<"found minimum at f("<<eps[0];
			 for(int i = 1; i < opt_dim; ++i)
				 std::cerr<<","<<eps[i];
			 std::cerr<< ")= "<<minf<<"\n";

			 for(int i = 1; i < opt_dim; ++i){
				 parameters[i]->value += eps[i];
			 }
		}

		std::cerr<<"Energy: " << _calc_expectation(&h) <<"\n";

		free(eps);
		free(lb);free(ub);
		nlopt_destroy(opt);

		if(options.compareWithClassicalEigenSolver){

			double id_term;
			Eigen::MatrixXd m = getIntermediateMatrixRepresentation(&h, &id_term);
			Eigen::VectorXcd evals = m.eigenvalues();

			std::vector<std::complex>(evals.size);
			for(int i = 0; i < evals.size(); ++i){

			}

		}


	}
	/*free(first_order_terms);
    for (int i = 0; i < parameters.size(); i++)
    	free(second_order_terms[i]);*/

}

PauliHamiltonian AqcPqcAccelerator::_calc_intermediate_hamiltonian(double lambda){

	std::vector<double> coeffs;

	for(std::vector<double>::size_type i = 0; i < hamil_int.coeffs0.size(); ++i){
		coeffs.push_back((1-lambda)*hamil_int.coeffs0[i]+lambda*hamil_int.coeffs1[i]);
	}

	PauliHamiltonian res(nbQubits, coeffs, hamil_int.pauliOpts);

	return res;
}

double AqcPqcAccelerator::calc_intermediate_expectation(ExperimentBuffer* buffer, double lambda, bool init_zero_state){

	if(lambda < 0 || lambda > 1)
		throw_runtime_error("Invalid lambda value");

	PauliHamiltonian h = _calc_intermediate_hamiltonian(lambda);
	return _calc_expectation(&h);

}

double AqcPqcAccelerator::_calc_expectation(PauliHamiltonian *h){
	initZeroState(qureg);
	run_circuit(ansatz.circuit);
	Qureg workspace = createQureg(nbQubits, env);
	PauliHamil hamil;
	h->toQuestPauliHamil(&hamil);
	return calcExpecPauliHamil(qureg, hamil, workspace);
}

Eigen::MatrixXd AqcPqcAccelerator::getIntermediateMatrixRepresentation(PauliHamiltonian* h, double* id_coeff){

	Eigen::MatrixXd m = Eigen::MatrixXd::Zero((1LL<<h->nbQubits), (1LL<<h->nbQubits));
	DiagonalOp diagOp = createDiagonalOp(h->nbQubits, env, true);

	PauliHamil h_cpy;
	h_cpy.numQubits = h->nbQubits;
	std::vector<double> coeffs;
	std::vector<int> codes;

	double x_coeff=0;
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
}

}
