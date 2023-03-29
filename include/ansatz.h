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

namespace FastVQA{

class Ansatz{
public:
	std::string name;
	Circuit circuit;
	int num_params;
	int depth;

	Ansatz(){}
	Ansatz(std::string name, int depth){
		this->name = name;
		this->depth = depth;
	}

};
Ansatz getAnsatz(std::string ansatz_type, int num_qubits, int depth=1, int seed=0);
void initOptimalParamsForMinusSigmaXHamiltonian(Ansatz *ansatz);

}

#endif /* FASTVQA_ANSATZ_H_ */
