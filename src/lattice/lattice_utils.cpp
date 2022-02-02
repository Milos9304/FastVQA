/*
 * quboToXvector.cpp
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#include "lattice.h"
#include "fplll/pruner/pruner.h"
#include <cmath>

#define log_pi 1.1447298858494

Hamiltonian Lattice::getHamiltonian(MapOptions* options){

	Hamiltonian result;
	this->calcHamiltonian(options, options->verbose);

	if(!qubo_generated){
		loge("Hamiltonian referenced but not yet generated!");
		throw;
	}

	result.nbQubits = expression_qubo->getIdMapSize()-1; //-1 bc of identity

	logw("Beware of overflows!! Create gmp qiskit or normalize to lower vals");
	logw("Normalization is probably best solution");

	for(auto &term : expression_qubo->polynomial){

		if(term.second == 0)
			continue;

		int id1 = term.first.first;
		int id2 = term.first.second;

		int pos1, pos2;

		result.coeffs.push_back(term.second.get_d());

		if(id1 == -1 && id2 == -1){ //id
			pos1=pos2=-1; //never matches
		}
		else if(id1 == -1){
			pos1=expression_qubo->getQubit(id2);
			pos2=pos1;
		}
		else if(id2 == -1){
			pos1=expression_qubo->getQubit(id1);
			pos2=pos1;
		}else{
			pos1=expression_qubo->getQubit(id1);
			pos2=expression_qubo->getQubit(id2);
		}

		for(int i = 0; i < result.nbQubits; ++i){
			if(i == result.nbQubits-1-pos1 || i == result.nbQubits-1-pos2)
				result.pauliOpts.push_back(3); //3 is pauli op code for z in qiskit
			else
				result.pauliOpts.push_back(0);
		}
	}
	return result;
}

mpq_class Lattice::calculate_gh_squared(MatrixInt* lattice){

	if(n_rows != n_cols){
		//loge("Gaussian heuristics not implemented for low rank matrices. Returning 1.");
		return 1;
	}

	int n = n_rows;


	if(!gso_orig_initialized){
		ZZ_mat<mpz_t> blank;

		gso_orig = new MatGSO<Z_NR<mpz_t>, FP_NR<double>>(orig_lattice, blank, blank, GSO_INT_GRAM);
		gso_orig->update_gso();

		gso_orig_initialized = true;
	}

	double ball_log_vol = (n/2.) * log_pi - lgamma(n/2. + 1);
	double log_vol = gso_orig->get_log_det(0, n).get_d();
	double log_gh =  1./n * (log_vol - 2 * ball_log_vol);

	return exp(log_gh);

}

/*
 *  @dst: destination vector
 *  @measurement: optimal config measurement x_n-1 x_n-2 ... x_0
 */
VectorInt Lattice::quboToXvector(std::string measurement){

	int n = measurement.size();

	std::map<int, int> penalized_varId_map;

	for (int i = 0; i < n; ++i) {
		int varId = qbit_to_varId_map[i];

		if(expression_penalized->getName(varId)[0] != 'z'){
			penalized_varId_map[varId] = measurement[n-1-i]-'0';
		}

	}

	penalized_varId_map[-1]=1; //to account for identity
	VectorInt res;
	for(auto &x: x_ids){
		mpq_class val = 0;

		for(auto &id_val:int_to_bin_map[x]){
			val += id_val.second * penalized_varId_map[id_val.first]; //binary var times its coeff
		}

		if(val.get_den() != 1){
			stringstream ss;
			ss << "Lattice::quboToXvector: decimal " << val << " to int conversion\n";
			logw(ss.str());
		}

		res.push_back(mpz_class(val));


	}
	return res;
}

VectorInt Lattice::quboToXvector(bool* measurement, int n){

	std::map<int, int> penalized_varId_map;

	for (int i = 0; i < n; ++i) {
		int varId = qbit_to_varId_map[i];

		if(expression_penalized->getName(varId)[0] != 'z'){
			penalized_varId_map[varId] = measurement[n-1-i]/*-'0'*/;
		}
	}

	VectorInt res;
	for(auto &x: x_ids){

		mpq_class val = 0;
		for(auto &id_val:int_to_bin_map[x]){
			val += id_val.second * penalized_varId_map[id_val.first]; //binary var times its coeff
		}

		if(val.get_den() != 1){
			stringstream ss;
			ss << "Lattice::quboToXvector: decimal " << val << " to int conversion\n";
			logw(ss.str());
		}

		res.push_back(mpz_class(val));


	}
	return res;
}

