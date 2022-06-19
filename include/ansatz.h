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

	Ansatz(std::string name){
		this->name = name;
	}

};
Ansatz getAnsatz(std::string ansatz_type, int num_qubits, int seed);
void initOptimalParamsForMinusSigmaXHamiltonian(Ansatz *ansatz);

}

#endif /* FASTVQA_ANSATZ_H_ */
