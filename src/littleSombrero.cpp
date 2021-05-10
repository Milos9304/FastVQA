/*
 * littleSombrero.cpp
 *
 *  Created on: May 10, 2021
 *      Author: Milos Prokop
 */

#include "littleSombrero.h"


std::vector<std::pair<MatrixInt, std::string>> LittleSombrero::loadLs(){

	std::vector<std::pair<MatrixInt, std::string>> res;

	bool success=true;

	logd("a");

	for(int i = 0; i < ls_num_instances; ++i){

		for(int r = rank_low; r <= rank_high; ++r){

			std::string file = "q_d_180_q_65537_r_"+std::to_string(r)+"_i_"+std::to_string(i);
			std::string filename = "../generated_lattices/littleSombrero/rank_"+std::to_string(r)+"/"+file;

			MatrixInt m = VqaConfig::loadLatticeFromFile(filename, &success);
			res.push_back(std::pair<MatrixInt, std::string>(m, file));

			if(!success){
				loge("Error loading " + filename);
				throw;
			}

		}

	}

	logd("All loaded successfilly");
	return res;

}



