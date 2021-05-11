/*
 * littleSombrero.cpp
 *
 *  Created on: May 10, 2021
 *      Author: Milos Prokop
 */

#include "littleSombrero.h"

std::vector<std::pair<double, double>> LittleSombrero::loadInfo(){

	std::vector<std::pair<double, double>> res;
	std::string delimiter = ", ";
	size_t pos = 0;
	std::string token;

	double lll, svp;

	for(int i = 0; i < ls_num_instances; ++i){

		std::string file = std::to_string(i) + "_info.txt";
		std::string filename = "../generated_lattices/littleSombrero/"+file;

		std::ifstream ifs(filename);

		std::string line;

		for(int r = rank_low; r <= rank_high; ++r){

			std::getline(ifs, line);
			line = line.substr(1, line.size()-2);

			pos = line.find(delimiter);
			token = line.substr(0, pos);
			lll = stod(token);
			line.erase(0, pos + delimiter.length());
			svp = stod(line);

			res.push_back(std::pair<double, double>(lll, svp));

		}
		ifs.close();

	}

	return res;

}

std::vector<std::pair<MatrixInt, std::string>> LittleSombrero::loadLs(){

	std::vector<std::pair<MatrixInt, std::string>> res;

	bool success=true;

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



