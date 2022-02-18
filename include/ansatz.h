/*
 * ansatz.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_ANSATZ_H_
#define FASTVQA_ANSATZ_H_

#include <random>
#include "circuit.h"

namespace fastVQA{
struct Ansatz{
	Circuit circuit;
	int num_params;
};
Ansatz getAnsatz(std::string ansatz_type, int num_qubits, int seed);
}

#endif /* FASTVQA_ANSATZ_H_ */
