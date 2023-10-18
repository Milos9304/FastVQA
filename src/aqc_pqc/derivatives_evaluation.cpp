#define _USE_MATH_DEFINES
#include <cmath>

#include "aqc_pqc/AqcPqcAccelerator.h"
#include <iostream>

namespace FastVQA{

	double AqcPqcAccelerator::calc_third_derivative(int i, int j, int k, void* constr_data){

		ConstrData_trivial *d = (ConstrData_trivial *) constr_data;
		qreal original_i = (*d->parameters)[i]->value;
		qreal original_j = (*d->parameters)[j]->value;
		qreal original_k = (*d->parameters)[k]->value;

		(*d->parameters)[i]->value += PI_2;
		(*d->parameters)[j]->value += PI_2;
		(*d->parameters)[k]->value += PI_2;
		qreal a = d->acc->_calc_expectation(d->h);

		(*d->parameters)[k]->value -= PI;
		qreal b = d->acc->_calc_expectation(d->h);

		(*d->parameters)[i]->value -= PI;
		(*d->parameters)[k]->value += PI;
		qreal c = d->acc->_calc_expectation(d->h);

		(*d->parameters)[k]->value -= PI;
		qreal dd = d->acc->_calc_expectation(d->h);

		(*d->parameters)[i]->value += PI;
		(*d->parameters)[j]->value -= PI;
		(*d->parameters)[k]->value += PI;
		qreal e = d->acc->_calc_expectation(d->h);

		(*d->parameters)[k]->value -= PI;
		qreal f = d->acc->_calc_expectation(d->h);

		(*d->parameters)[i]->value -= PI;
		(*d->parameters)[k]->value += PI;
		qreal g = d->acc->_calc_expectation(d->h);

		(*d->parameters)[k]->value -= PI;
		qreal h = d->acc->_calc_expectation(d->h);

		(*d->parameters)[i]->value = original_i;
		(*d->parameters)[j]->value = original_j;
		(*d->parameters)[k]->value = original_k;

		return (a-b-c+dd-e+f+g-h)/8;

	};

}
