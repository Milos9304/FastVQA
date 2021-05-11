/*
 * littleSombrero.h
 *
 *  Created on: May 10, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LITTLESOMBRERO_H_
#define SRC_LITTLESOMBRERO_H_

#include "vqaConfig.h"

class LittleSombrero{

public:

	const static int ls_num_instances = 50;
	const static int rank_low = 10;
	const static int rank_high = 45;

	std::vector<std::pair<MatrixInt, std::string>> loadLs();
	std::vector<std::pair<double, double>> loadInfo();
};

#endif /* SRC_LITTLESOMBRERO_H_ */
