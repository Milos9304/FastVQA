/*
 * iterativeLatticeReduction.cpp
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#include "iterativeLatticeReduction.h"

void IterativeLatticeReduction::run(){

	performLLLonLattice();

	for(int i = 0; i < n_iters; ++i){
		run_quantum();
		performLLLonLattice();
	}

}

void IterativeLatticeReduction::run_quantum(){
	this->quantum_oracle();
}
