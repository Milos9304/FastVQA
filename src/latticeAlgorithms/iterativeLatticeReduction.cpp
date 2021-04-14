/*
 * iterativeLatticeReduction.cpp
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#include "iterativeLatticeReduction.h"

void IterativeLatticeReduction::run(){

	for(int i = 0; i < n_iters; ++i){
		performLLLonLattice();
		run_quantum();
	}

}

void IterativeLatticeReduction::run_quantum(){

	xacc::qbit* buffer;
	this->quantum_oracle(&buffer, lattice->toHamiltonianString(options), lattice->name);

	std::cout << "Min QUBO: " << (*buffer)["opt-val"].as<double>() << "\n";
	std::cerr << "Opt config: " << (*buffer)["opt-config"].as<std::string>() << "\n";
	throw;
}
