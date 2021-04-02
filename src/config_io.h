/*
 * config_io.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_CONFIG_IO_H_
#define SRC_CONFIG_IO_H_

#include "lattice.h"
#include <vector>
#include <fstream>

class VqaConfig{

	public:

		bool verbose = true;

		VqaConfig(std::string pathname);

		std::vector<Lattice> getLattices(){
			return lattices;
		}

		void setCurrentHamiltonian(std::string hamiltonian){
			current_hamiltonian = hamiltonian;
		}

		std::string getCurrentHamiltonian(){
			return current_hamiltonian;
		}

	private:

		std::vector<Lattice> lattices;
		std::string current_hamiltonian = "";

		MatrixInt loadLatticeFromFile(std::string filename, bool *success);



};



#endif /* SRC_CONFIG_IO_H_ */
