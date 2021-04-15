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
		std::pair<std::string, double> opt_config = run_quantum();




	}

}

void IterativeLatticeReduction::run_test(){

	for(int i = 0; i < n_iters; ++i){

		//performLLLonLattice();
		std::pair<std::string, double> opt_config = run_quantum();

		lattice->quboToXvector(opt_config.first);

	}

}


std::pair<std::string, double> IterativeLatticeReduction::run_quantum(){

	xacc::qbit* buffer;
	this->quantum_oracle(&buffer, lattice->toHamiltonianString(options), lattice->name);

	//std::cerr << "Min QUBO: " << (*buffer)["opt-val"].as<double>() << "\n";
	//std::cerr << "Opt config: " << (*buffer)["opt-config"].as<std::string>() << "\n";

	return std::pair<std::string, double>((*buffer)["opt-config"].as<std::string>(),
			(*buffer)["opt-val"].as<double>());

	//std::cout << "Min QUBO: " << (*buffer)["opt-val"].as<double>() << "\n";
	//std::cerr << "Opt config: " << (*buffer)["opt-config"].as<std::string>() << "\n";
}
