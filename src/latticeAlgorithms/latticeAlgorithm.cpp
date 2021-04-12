/*
 * latticeAlgorithm.cpp
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#include "latticeAlgorithm.h"
#include "fplll.h"

void LatticeAlgorithm::performLLLonLattice(){

	//lll_reduction(lattice->get_orig_lattice(), delta, eta, method, floatType, precision, flags)
	lll_reduction(*lattice->get_orig_lattice(), /*new lattice here*/0.9, 0.9, LLLMethod::LM_FAST, FloatType::FT_DOUBLE);

}


