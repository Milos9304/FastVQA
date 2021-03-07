/*
 * config_io.cpp
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#include "config_io.h"
#include "logger.h"

VqaConfig::VqaConfig(std::string pathname){

	std::string hml_file;

	std::ifstream ifs(pathname);
	std::string line;
	std::istringstream line_stream;

	if(!ifs.is_open()){
		throw std::runtime_error("Unable to open config file " + pathname);
	}

	bool lookForType = true;

	/*
	 * LOAD CONFIGURATION FILE
	 * */

	while (std::getline(ifs, line)) {

		if(line[0] == '#')
				continue;

		line_stream.str(line.substr(line.find("=")+1));

		if(lookForType and line.find("type") == std::string::npos)
			throw std::runtime_error("Invalid config file type.");

		lookForType = false;

		if (line.find("type") != std::string::npos) {
			if(line_stream.str() != "qaoa")
				throw std::runtime_error("Invalid config file type.");
		}

		else if (line.find("hamiltonian_file") != std::string::npos) {
			line_stream >> hml_file;
		}

		else if (line.find("verbose") != std::string::npos) {
			std::string temp;
			line_stream >> temp;
			if(temp == "true")
				verbose = true;
		}

		line_stream.clear();
	}

	ifs.close();

	/*
	 * LOAD CONTENTS OF HAMILTONIAN FILE
	 */

	if(hml_file.empty())
		throw std::runtime_error("'hamiltonian_file' not specified in the config file");

	std::ifstream hml_ifs(hml_file);
	if(!hml_ifs.is_open()){
			throw std::runtime_error("Unable to open hamiltonian file " + hml_file);
	}

	while (std::getline(hml_ifs, line))
		hamiltonians.push_back(line);

	hml_ifs.close();

}
