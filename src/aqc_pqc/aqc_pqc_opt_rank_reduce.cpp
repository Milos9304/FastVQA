#define _USE_MATH_DEFINES
#include <cmath>

#include "aqc_pqc/AqcPqcAccelerator.h"
#include <nlopt.h>
#include <iostream>

namespace FastVQA{

double ineq_constraint_rank_reduce(unsigned n, const double *x, double *grad, void *data){
		//std::cerr<<"e";
		//my_constraint_data *d = (my_constraint_data *) data;
		//double a = d->a, b = d->b;

		ConstrData_rank_reduce *d = (ConstrData_rank_reduce *) data;
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


		Eigen::SelfAdjointEigenSolver<Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic>> solver(H.rows());
		solver.compute(H);
		Eigen::Vector<qreal, Eigen::Dynamic> lambda = solver.eigenvalues().reverse();
		auto X = solver.eigenvalues();
		auto E = solver.eigenvectors();

		//std::cerr<<"pass: "<<X.col(0)[0]<<"\n";
		//std::cerr<<"min eval: "<<X.col(0)[0]<<"\n";

		if (grad) {
			std::cerr<<"e_grad";
			std::cerr<<E.col(0)<<std::endl;throw;
			//for(unsigned int k = 0; k < n; ++k){
			//	//grad[k] =
			//}

		}

		return -X.col(0)[0];//pass?-1:1;
	}

	double lin_system_f_rank_reduce(unsigned n, const double *z, double *grad, void *data){
		OptData_rank_reduce *d = (OptData_rank_reduce *) data;
		//std::cerr<< "probably wrong derivative";
		//std::cerr<<"g";
		if (grad) {
			std::cerr<<"grad";
			for(unsigned int k = 0; k < n; ++k){
				qreal s = 0;
				for(int j = 0; j < d->Xi.rows(); ++j){
					double s1 = 0;
					for(unsigned int l = 0; l < n; ++l){
						s1 += d->N(j,l)*z[l];
					}
					s += d->Xi[j]*d->N(j,k) + 2*s1*d->N(j,k)*z[k];
				}
				grad[k] = (double)s;
			}
		}

		qreal ret=0;
		Eigen::Vector<qreal, Eigen::Dynamic> z_vect(n);
		for(unsigned int i = 0; i < n; ++i)
			z_vect(i)=z[i];
		//std::cerr<<z_vect<<" ";
		//std::cerr<<d->Xi<<" ";
		//std::cerr<<"N "<<d->N<<" koniec";
		Eigen::Vector<qreal, Eigen::Dynamic> x = d->Xi+d->N*z_vect;
		for(int i = 0; i < d->Xi.rows(); ++i)
			ret += x[i]*x[i];//std::cerr<<"-"<<ret<<"-";
		//std::cerr<<std::setprecision (50)<<" "<<ret<<"\n";
		//std::cerr<<ret<<" ";
		return ret;
	}

	/*double ineq_constraint(unsigned n, const double *x, double *grad, void *data){
		return AqcPqcAccelerator::ineq_constraint(n,x,grad,data);
	};*/

	Eigen::Vector<qreal, Eigen::Dynamic> AqcPqcAccelerator::_optimize_with_rank_reduction(PauliHamiltonian *h, Eigen::Vector<qreal, Eigen::Dynamic> *Q, Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> *A, std::vector<std::shared_ptr<Parameter>> *parameters){



		Eigen::Vector<qreal, Eigen::Dynamic> Xi(Q->rows());
		Xi = A->fullPivHouseholderQr().solve(*Q);
		//std::cerr<<"A: "<<A<< ""<<" -Q: "<<Q<<" Xi: "<<Xi<<std::endl;
		//bool solution_exists = (A*Xi).isApprox(Q, 10e-4);
		//std::cerr<<"exists:"<<solution_exists<<" ";

		//Eigen::FullPivLU<Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic>> lu(*A);


		//lu.setThreshold(this->options.roundDecimalPlaces == -1 ? 10e-6 : pow(10, -(this->options.roundDecimalPlaces+1)));
		//lu.setThreshold(threshold);
		//std::cerr<<"Rank = " << lu.rank()<<std::endl;
		//Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> A_null_space = lu.kernel();
		//for faster way see https://stackoverflow.com/questions/34662940/how-to-compute-basis-of-nullspace-with-eigen-library

		Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic>> cod;

		double rank;
		double threshold = 0;
		do{
			threshold += 0.01;
			cod.setThreshold(threshold);
			cod.compute(*A);
			rank = cod.rank();
		}while(rank == 10);
		std::cerr<<"Rank = " << rank <<std::endl;


		// Find URV^T
		Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> V = cod.matrixZ().transpose();
		//std::cerr<<"V:"<<V<<std::endl;std::cerr<<cod.rank()<<" "<<V.rows() <<  V.cols() <<V.cols() - cod.rank() <<std::endl;
		Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> A_null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
		Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> P = cod.colsPermutation();
		//std::cout << "The null space: \n" << A_null_space << "\n" ;

		A_null_space = P * A_null_space; // Unpermute the columns
		//std::cerr<<"A:"<<*A<<"\n"<<"V:"<<V<<"\n"<<"AN:"<<A_null_space<<"\nP:"<<P<<"\n";
		// The Null space:
		//std::cout << "The null space: \n" << A_null_space << "\n" ;
		// Check that it is the null-space:
		//std::cout << "mat37 * Null_space = \n" << (*A) * A_null_space  << '\n';

		int opt_dim = A_null_space.cols();
		//std::cerr<<A_null_space.rows() << " " << A_null_space.cols()<<"\n";
		//std::cerr<<A_null_space;
		std::cerr<<"d:"<<opt_dim<<std::endl;

		double *lb, *ub;
		OptData_rank_reduce data {Xi, A_null_space};

		nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, opt_dim);

		//opt = nlopt_create(NLOPT_LN_AUGLAG_EQ, opt_dim);
		//lopt = nlopt_create(NLOPT_LN_COBYLA, opt_dim);
		//nlopt_set_local_optimizer(opt, lopt);
		ConstrData_rank_reduce constr_data {parameters, this, h, &data};
		//nlopt_add_equality_constraint(opt, eq_constraint, &constr_data, 0);
		nlopt_add_inequality_constraint(opt, ineq_constraint_rank_reduce, &constr_data, 0);
		lb = (double*) malloc(opt_dim * sizeof(double));
		ub = (double*) malloc(opt_dim * sizeof(double));
		for(int i = 0; i < opt_dim; ++i){
			lb[i] = -HUGE_VAL;ub[i]=-lb[i];
		}
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
		nlopt_set_min_objective(opt, lin_system_f_rank_reduce, &data);
		nlopt_set_xtol_rel(opt, options.xtol);
		nlopt_set_xtol_abs1(opt, options.xtol);
		nlopt_set_maxtime(opt, 90);

		//		}else
		//			throw;

		double *eps = (double*) malloc(opt_dim * sizeof(double));
		for(int i = 0; i < opt_dim; ++i)
			eps[i] = 0;
		double minf;
		int opt_res = nlopt_optimize(opt, eps, &minf);
		if (opt_res < 0) {
			std::cerr<<"nlopt failed! error code: "<< opt_res <<std::endl;
			return Xi;
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

			//std::cerr<<"Solution found: " << Xi+A_null_space*eps_vect << "\n";
			//std::cerr<<"Min eval: " << ineq_constraint_rank_reduce(opt_dim, eps, NULL, &constr_data)<<"\n";
			//std::cerr<<"This should be zero: " << (*A)*(Xi+A_null_space*eps_vect)-(*Q) << "\n";
			//std::cerr<<"A_null:"<<A_null_space;

			free(eps);
			free(lb);free(ub);
			nlopt_destroy(opt);
			return Xi+A_null_space*eps_vect;

			//std::cerr<<res_eps;//throw;

			//std::cerr<<"x"<<A*res_eps-Q<<std::endl;
		}
	}

}
