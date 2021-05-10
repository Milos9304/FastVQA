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

	int ls_num_instances = 50;
	int rank_low = 10;
	int rank_high = 45;

	std::vector<std::pair<MatrixInt, std::string>> loadLs();
};

#endif /* SRC_LITTLESOMBRERO_H_ */
