/*
 * saveProgress.cpp
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#include "../logger.h"
#include "saveProgress.hpp"

save_instance::save_instance() {} // @suppress("Class members should be properly initialized")

void saveProgress(std::string filename, std::vector<double> coefficients, double expected_energy, std::string sv, double sv_energy, double hit_rate){

	std::ofstream ofs(filename);
	boost::archive::text_oarchive ar(ofs);
	save_instance instance(coefficients, expected_energy, sv, sv_energy, hit_rate);
	ar & instance;
	ofs.close();

	std::ofstream pyofs(filename+"_py");
	pyofs << expected_energy << "\n";
	pyofs << sv << "\n";
	pyofs << sv_energy << "\n";
	pyofs << hit_rate;
	pyofs.close();

}


void saveProgress(std::string filename, std::vector<double> coefficients, double expected_energy, double sv_energy, double hit_rate){

	std::ofstream ofs(filename);
	boost::archive::text_oarchive ar(ofs);
	save_instance instance(coefficients, expected_energy, sv_energy, hit_rate);
	ar & instance;
	ofs.close();

	std::ofstream pyofs(filename+"_py");
	pyofs << expected_energy << "\n";
	pyofs << sv_energy << "\n";
	pyofs << hit_rate;
	pyofs.close();

}

bool loadProgress(std::string filename, std::vector<double>* coefficients, double* expected_energy, double* sv_energy, double* hit_rate){

	save_instance instance;
	std::ifstream ifs(filename);
	if(!ifs.is_open()){
		loge(filename + " not found");
		return false;
	}

	boost::archive::text_iarchive ar(ifs);
	ar & instance;
	ifs.close();

	*coefficients = instance.coefficients;
	*expected_energy = instance.expected_energy;
	*sv_energy = instance.sv_energy;
	*hit_rate = instance.hit_rate;

	return true;

}


