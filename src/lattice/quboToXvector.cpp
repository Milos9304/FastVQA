/*
 * quboToXvector.cpp
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#include "lattice.h"

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

		//std::cout << val << " ";

	}//std::cout<<"\n";

	//redo z0=1, z1=x1 - NOT NEEDED AS WE CARE ONLY ABOUT X
	//penalized_varId_map[z0_id] = 1;
	//penalized_varId_map[z1_id] = x1_value;

	return res;

}

//substitute back z0=1, z1=x1
/*void Lattice::quboToPenalized(int* dst){



}*/
