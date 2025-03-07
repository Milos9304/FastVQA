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
#include <Eigen/Dense>
#include "QuEST.h"

namespace FastVQA{

enum PauliHamiltonianType {General, MinusSigmaX, SumMinusSigmaX};

class PauliHamiltonian{

	public:

		int nbQubits;
		PauliHamiltonianType type = PauliHamiltonianType::General;

		std::string getPauliHamiltonianString(int double_precision=2);
		Eigen::MatrixXd getMatrixRepresentation(bool diagonalOp);
		Eigen::Matrix<qreal, Eigen::Dynamic, Eigen::Dynamic> getMatrixRepresentation2(bool diagonalOp);

		//quest formulation
		std::vector<qreal> coeffs;
		std::vector<int> pauliOpts;

		void initializeMinusSigmaXHamiltonian();
		void initializeSumMinusSigmaXHamiltonian();

		PauliHamiltonian(){}
		PauliHamiltonian(int nbQubits){
			this->nbQubits = nbQubits;
		}
		PauliHamiltonian(int nbQubits, std::vector<qreal> coeffs, std::vector<int> pauliOpts){
			this->nbQubits = nbQubits;
			this->coeffs = coeffs;
			this->pauliOpts = pauliOpts;
		}

		void toQuestPauliHamil(PauliHamil* hamil);

		void to_ising_file(std::string filename);

	    //if empty, ONLY ground states of custom_solutions are considered as valid solution
		//otherwise, the following vector contains string of the form 01.100.1. where . can be 0/1 and the string marks a valid solution
		std::vector<std::string> custom_solutions;

};

}

#endif /* FASTVQA_HAMILTONIAN_H_ */
