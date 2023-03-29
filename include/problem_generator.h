/*
 * problem_generator.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_PROBLEM_GENERATOR_H_
#define FASTVQA_PROBLEM_GENERATOR_H_

#include "hamiltonian.h"

namespace FastVQA{

	class ProblemGenerator{

		Hamiltonian genRandomMaxCutHamiltonian(int nbQubits);

	};

}
#endif /* SRC_OPTIMIZER_H_ */
