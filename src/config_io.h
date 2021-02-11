/*
 * config_io.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_CONFIG_IO_H_
#define SRC_CONFIG_IO_H_

#include <vector>
#include <fstream>

struct VqeConfig{

	std::vector<std::string> hamiltonians;

};

extern VqeConfig loadConfigFile(std::string pathname);

#endif /* SRC_CONFIG_IO_H_ */
