/*
 * iterativeLatticeReduction.cpp
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#include "iterativeLatticeReduction.h"
#include "fplll.h"

void IterativeLatticeReduction::run(){

	for(int i = 0; i < n_iters; ++i){

		performLLLonLattice(0.99, 0.51);
		std::pair<std::string, double> opt_config = run_quantum();

		double short_vector_len_sq = opt_config.second;


	}

}

//x_vect, energy
std::pair<VectorInt, double> IterativeLatticeReduction::run_test(){

	std::pair<std::string, double> opt_config = run_quantum();
	VectorInt x_vect = lattice->quboToXvector(opt_config.first);

	return std::pair<VectorInt, double>(x_vect, opt_config.second);

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
