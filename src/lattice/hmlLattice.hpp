/*
 * hmlLattice.hpp
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICE_HMLLATTICE_HPP_
#define SRC_LATTICE_HMLLATTICE_HPP_

#include "abstractLatticeInput.hpp"

class HmlLattice : public AbstractLatticeInput{

	public:

		HmlLattice(int n, std::string hamiltonian){
			this->n = n;
			this->hamiltonian = hamiltonian;
		}

		VectorInt quboToXvector(std::string measurement){
			VectorInt r;
			return r;
		}

		std::string toHamiltonianString(){return hamiltonian;}

		std::pair<std::vector<double>, std::vector<int>> getHmlInQuestFormulation(){}



	private:
		std::string hamiltonian;

};



#endif /* SRC_LATTICE_HMLLATTICE_HPP_ */
