/*
 * nlopt_test.cpp
 *
 *  Created on: September 5, 2022
 *      Author: Milos Prokop
 */

#include <gtest/gtest.h>
#include <cmath>
#include "fastVQA.h"
#include <nlopt.h>


#define GTEST_COUT std::cerr << "[          ] [ INFO ]"
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
	std::cerr<<"a";

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    std::cerr<<"b";
    return sqrt(x[1]);
}
typedef struct {
    double a, b;
} my_constraint_data;
double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
 }
TEST(nlopt_test, initializeMinusSigmaXHamiltonian){

	double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LN_COBYLA, 2); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc, NULL);
	my_constraint_data data[2] = { {2,0}, {-1,1} };
	nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
	nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);
	nlopt_set_xtol_rel(opt, 1e-4);
	double x[2] = { 1.234, 5.678 };  /* `*`some` `initial` `guess`*` */
	double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */
	if (nlopt_optimize(opt, x, &minf) < 0) {
	   std::cerr<<"nlopt failed!\n";
	}
	else {
		std::cerr<<"found minimum at "<< x[0] <<","<<x[1]<<"  : "<< minf<<"\n";
	}
}

int main(int argc, char **argv) {

	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	return ret;


}
