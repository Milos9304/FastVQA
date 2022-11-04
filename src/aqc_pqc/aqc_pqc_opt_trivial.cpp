#define _USE_MATH_DEFINES
#include <cmath>

#include "aqc_pqc/AqcPqcAccelerator.h"
#include <nlopt.h>
#include <iostream>

namespace FastVQA{

typedef struct {
		Eigen::Vector<qreal, Eigen::Dynamic> *Q;
	    Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic>* A;
} OptData;

typedef struct {
	std::vector<std::shared_ptr<Parameter>> *parameters;
	AqcPqcAccelerator *acc;
	PauliHamiltonian *h;
	OptData *optData;
} ConstrData;

	double ineq_constraint_trivial(unsigned n, const double *x, double *grad, void *data){

		ConstrData *d = (ConstrData *) data;
		Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> H(d->parameters->size(), d->parameters->size());

		Eigen::Vector<qreal, Eigen::Dynamic> res_eps(n);
		for(unsigned int i = 0; i < n; ++i)
			res_eps(i)=x[i];

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


		Eigen::SelfAdjointEigenSolver<Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic>> solver(H.rows());
		solver.compute(H);
		Eigen::Vector<qreal, Eigen::Dynamic> lambda = solver.eigenvalues().reverse();
		auto X = solver.eigenvalues();

		//std::cerr<<"pass: "<<X.col(0)[0]<<"\n";
		//std::cerr<<"min eval: "<<X.col(0)[0]<<"\n";


		if (grad) {
			for(unsigned int i = 0; i < d->parameters->size(); ++i){

				double e0_aij = 0;


			}
			std::cerr<<solver.eigenvectors().col(0)<<std::endl;
			throw;
			//UNIMPLEMENTED

		}

		return -X.col(0)[0];//pass?-1:1;
	}

	double lin_system_f_trivial(unsigned n, const double *z, double *grad, void *data){
			OptData *dt = (OptData *) data;

			Eigen::Vector<qreal, Eigen::Dynamic> z_vect(n);
			for(unsigned int i = 0; i < n; ++i)
				z_vect(i)=z[i];

			Eigen::Vector<qreal, Eigen::Dynamic> x = *(dt->Q)+*(dt->A)*z_vect;

			Eigen::Vector<qreal, Eigen::Dynamic> gram_m = (*(dt->A)).transpose() * (*(dt->A));

			//in trivial, A is symmetric, so nxn
			if (grad) {
				for(int d = 0; d < n; ++d){
					double g1 = 0, g2 = 0;
					for(int k = 0; k < n; ++k){
						g1+=z[k]*gram_m(d,k);
						g2+=(*(dt->A))(k,d);
					}
					grad[d]=2*(g1+(*(dt->Q))[d]*g2);
				}
			}
			return (x.transpose()*x)(0,0);
		}

	Eigen::Vector<qreal, Eigen::Dynamic> AqcPqcAccelerator::_optimize_trivially(PauliHamiltonian *h, Eigen::Vector<qreal, Eigen::Dynamic> *Q, Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> *A, std::vector<std::shared_ptr<Parameter>> *parameters){

		int opt_dim = Q->rows();
		nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, opt_dim);
		//nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, opt_dim);
		OptData data {Q, A};

		//std::cerr<<"A: " << *A << "\n" << "q: " << -(*Q)<<std::endl;throw;

	    ConstrData constr_data {parameters, this, h, &data};
		//nlopt_add_equality_constraint(opt, eq_constraint, &constr_data, 0);
		nlopt_add_inequality_constraint(opt, ineq_constraint_trivial, &constr_data, 0.0002);
		double* lb = (double*) malloc(opt_dim * sizeof(double));
		double* ub = (double*) malloc(opt_dim * sizeof(double));
		for(int i = 0; i < opt_dim; ++i){
			lb[i] = -HUGE_VAL;ub[i]=-lb[i];
		}
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
		nlopt_set_min_objective(opt, lin_system_f_trivial, &data);

		//nlopt_set_ftol_rel(opt, options.xtol);
		//nlopt_set_ftol_abs(opt, options.xtol);

		nlopt_set_xtol_rel(opt, options.xtol);
		nlopt_set_xtol_abs1(opt, options.xtol);

		nlopt_set_maxtime(opt, 90);
		//nlopt_set_maxeval(opt, 500);
		double *eps = (double*) malloc(opt_dim * sizeof(double));
		for(int i = 0; i < opt_dim; ++i)
			eps[i] = 0;
		double minf;
		int opt_res = nlopt_optimize(opt, eps, &minf);
		if (opt_res < 0) {
			std::cerr<<"nlopt failed! error code: "<< opt_res <<std::endl;
		}
		else {
			std::cerr<<"nlopt success code: "<< opt_res << std::endl;
			 //std::cerr<<"found minimum at f("<<eps[0];
			 //for(int i = 1; i < opt_dim; ++i)
			//	 std::cerr<<","<<eps[i];
			 //std::cerr<< ")= "<<minf<<"\n";

			Eigen::Vector<qreal, Eigen::Dynamic> eps_vect(opt_dim);
			for(int i = 0; i < opt_dim; ++i)
				eps_vect(i)=eps[i];

			//std::cerr<<eps_vect<<std::endl;
			//std::cerr<<"A_null:"<<A_null_space;

			free(eps);
			free(lb);free(ub);
			nlopt_destroy(opt);
			return eps_vect;

			//std::cerr<<res_eps;//throw;

			//std::cerr<<"x"<<A*res_eps-Q<<std::endl;
		}


	}

}
