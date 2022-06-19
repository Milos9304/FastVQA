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
#include "QuEST.h"

namespace FastVQA{

class Hamiltonian{

	public:

		int nbQubits;

		std::string getHamiltonianString(){return "";}

		//quest formulation
		std::vector<double> coeffs;
		std::vector<int> pauliOpts;

		void initializeMinusSigmaXHamiltonian();

		Hamiltonian(){}
		Hamiltonian(int nbQubits){
			this->nbQubits = nbQubits;
		}
		Hamiltonian(int nbQubits, std::vector<double> coeffs, std::vector<int> pauliOpts){
			this->nbQubits = nbQubits;
			this->coeffs = coeffs;
			this->pauliOpts = pauliOpts;
		}

		void toPauliHamil(PauliHamil* hamil);
};

}

#endif /* FASTVQA_HAMILTONIAN_H_ */
