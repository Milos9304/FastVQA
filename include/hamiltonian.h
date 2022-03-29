/*
 * hamiltonian.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_HAMILTONIAN_H_
#define FASTVQA_HAMILTONIAN_H_

#include <string>
#include <vector>

namespace fastVQA{

struct Hamiltonian{
	int nbQubits;

	std::string getHamiltonianString(){return "";}

	//quest formulation
	std::vector<double> coeffs;
	std::vector<int> pauliOpts;
};

}

#endif /* FASTVQA_HAMILTONIAN_H_ */
