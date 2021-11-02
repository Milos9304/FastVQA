/*
 * hamiltonian.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_HAMILTONIAN_H_
#define SRC_HAMILTONIAN_H_

#include <string>
#include <vector>

struct Hamiltonian{
	int nbQubits;

	std::string getHamiltonianString(){return "";}

	//quest formulation
	std::vector<double> coeffs;
	std::vector<int> pauliOpts;
};

//Hamiltonian getHamiltonian(std::string Hamiltonian_type, int num_qubits);
#endif /* SRC_HAMILTONIAN_H_ */
