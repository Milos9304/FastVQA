/*
 * utils.cpp
 *
 *  Created on: Apr 14, 2021
 *      Author: Milos Prokop
 */

#include "utils.hpp"
#include <cmath>

void OnlineMeanVarianceCalculator::update(double new_sample){

	mean = 1./(n+1)*(n*mean+new_sample);
	variance = double(n)/(n+1)*variance+pow((new_sample-mean),2)/double(n);

}
