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

		std::string toHamiltonianString(){return hamiltonian;}

	private:
		std::string hamiltonian;

};



#endif /* SRC_LATTICE_HMLLATTICE_HPP_ */
