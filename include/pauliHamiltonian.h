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
#include <Eigen/Dense>

namespace FastVQA{

class PauliHamiltonian{

	public:

		int nbQubits;

		std::string getPauliHamiltonianString(int double_precision=2);
		Eigen::MatrixXd getMatrixRepresentation(bool diagonalOp=false);

		//quest formulation
		std::vector<double> coeffs;
		std::vector<int> pauliOpts;

		void initializeMinusSigmaXHamiltonian();

		PauliHamiltonian(){}
		PauliHamiltonian(int nbQubits){
			this->nbQubits = nbQubits;
		}
		PauliHamiltonian(int nbQubits, std::vector<double> coeffs, std::vector<int> pauliOpts){
			this->nbQubits = nbQubits;
			this->coeffs = coeffs;
			this->pauliOpts = pauliOpts;
		}

		void toQuestPauliHamil(PauliHamil* hamil);
};

}

#endif /* FASTVQA_HAMILTONIAN_H_ */
