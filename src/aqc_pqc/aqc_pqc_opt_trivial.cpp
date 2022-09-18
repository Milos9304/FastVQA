#define _USE_MATH_DEFINES
#include <cmath>

#include "aqc_pqc/AqcPqcAccelerator.h"
#include <nlopt.h>
#include <iostream>

namespace FastVQA{

	Eigen::Vector<qreal, Eigen::Dynamic> AqcPqcAccelerator::_optimize_trivially(PauliHamiltonian *h, Eigen::Vector<qreal, Eigen::Dynamic> *minus_q, Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> *A, std::vector<std::shared_ptr<Parameter>> *parameters){

		int opt_dim = minus_q->rows();
		nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, opt_dim);
		OptData data {Xi, A_null_space};
	    ConstrData constr_data {parameters, this, h, &data};
		//nlopt_add_equality_constraint(opt, eq_constraint, &constr_data, 0);
		nlopt_add_inequality_constraint(opt, ineq_constraint, &constr_data, 0);
		/*lb = (double*) malloc(opt_dim * sizeof(double));
		ub = (double*) malloc(opt_dim * sizeof(double));
		for(int i = 0; i < opt_dim; ++i){
			lb[i] = -HUGE_VAL;ub[i]=-lb[i];
		}
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
		nlopt_set_min_objective(opt, lin_system_f, &data);
		nlopt_set_xtol_rel(opt, 1e-1);
		nlopt_set_xtol_abs1(opt, 1e-1);*/

	}

}
