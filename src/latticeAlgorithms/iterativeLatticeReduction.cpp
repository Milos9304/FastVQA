/*
 * iterativeLatticeReduction.cpp
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#include "iterativeLatticeReduction.h"
#include "fplll.h"
#include <cmath>

void IterativeLatticeReduction::run(){

	for(int i = 0; i < n_iters; ++i){

		std::pair<std::string, double> opt_config = run_quantum();

		VectorInt x_vect = lattice->quboToXvector(opt_config.first);
		loge("NOW FOLLOWS LLL X VECTOR");
		for(auto &x:x_vect)
			std::cerr<<x<<" ";
		std::cerr<<"\n";

		double short_vector_len_sq = opt_config.second;
		logw("short vector len sq: " + std::to_string(short_vector_len_sq));

		std::cerr<< "gh^2 = " << lattice->get_orig_gh() << "\n";

		int sum = 0;
		for(auto &x: lattice->get_current_lattice()->matrix[0])
				sum += x.get_si()*x.get_si();
		std::cerr<< "LLL |b1|: norm/gh = " << sqrt(sum) / sqrt(lattice->get_orig_gh().get_d()) << "\n";
		std::cerr<< "final |b1|: norm/gh = " << sqrt(short_vector_len_sq) / sqrt(   lattice->get_orig_gh().get_d()  ) << "\n";

		loge("NOW FOLLOWS FINAL X VECTOR");

		for(int c = 0; c < lattice->lll_transformation->c; ++c){
			int res = 0;
			for(unsigned int r = 0; r < x_vect.size(); ++r){
				int a = lattice->lll_transformation->matrix[r][c].get_si();
				int b = x_vect[r].get_si();
				res += b * a;
			}
			std::cerr<<res<<" ";
		}std::cerr<<"\n";


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
