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

mpq_class Lattice::calculate_gh_squared(MatrixInt* lattice){

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

	int x1_value;

	for (int i = 0; i < n; ++i) {
		int varId = qbit_to_varId_map[i];

		if(expression_penalized->getName(varId)[0] != 'z'){
			penalized_varId_map[varId] = measurement[n-1-i]-'0';
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

