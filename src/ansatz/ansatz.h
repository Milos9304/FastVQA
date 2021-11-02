/*
 * ansatz.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_ANSATZ_H_
#define SRC_ANSATZ_H_

#include <random>
#include "../circuit/Circuit.h"

struct Ansatz{
	Circuit circuit;
	int num_params;
};

Ansatz getAnsatz(std::string ansatz_type, int num_qubits, int seed);
#endif /* SRC_ANSATZ_H_ */
