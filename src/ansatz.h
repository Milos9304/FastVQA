/*
 * ansatz.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_ANSATZ_H_
#define SRC_ANSATZ_H_

#include "Circuit.h"

struct Ansatz{
	Circuit circuit;
	Parameters parameters;
};

Ansatz getAnsatz(std::string ansatz_type, int num_qubits);
#endif /* SRC_ANSATZ_H_ */
