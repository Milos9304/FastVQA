/*
 * config_io.cpp
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#include "config_io.h"
#include "logger.h"

extern VqeConfig loadConfigFile(std::string pathname){

	std::string ipName;
	int nR, nC;

	std::ifstream fin(pathname);
	std::string line;
	std::istringstream sin;

	if(!fin.is_open()){
		throw std::runtime_error("Unable to open file " + pathname);
	}

	while (std::getline(fin, line)) {

	 sin.str(line.substr(line.find("=")+1));

	 if (line.find("Input name") != std::string::npos) {
	  logd("Input name " + sin.str());
	  //sin >> ipName;
	 }
	 /*else if (line.find("Num. of rows") != std::string::npos) {
	  sin >> nR;
	 }
	 else if (line.find("Num. of cols") != std::string::npos) {
	  sin >> nC;
	 }*/
	 sin.clear();

	}

	return VqeConfig();

}
