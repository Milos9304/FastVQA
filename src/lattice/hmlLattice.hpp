/*
 * hmlLattice.hpp
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

/*#ifndef SRC_LATTICE_HMLLATTICE_HPP_
#define SRC_LATTICE_HMLLATTICE_HPP_

#include "abstractLatticeInput.hpp"

class HmlLattice{

	public:

		HmlLattice(int n_rows, int n_cols, std::string hamiltonian){
			this->n_rows = n_rows;
			this->n_cols = n_cols;
			this->hamiltonian = hamiltonian;
		}

		VectorInt quboToXvector(std::string measurement){
			VectorInt r;
			return r;
		}

		std::string toHamiltonianString(){return hamiltonian;}

	private:
		std::string hamiltonian;

};



#endif *//* SRC_LATTICE_HMLLATTICE_HPP_ */
