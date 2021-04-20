/*
 * config_io.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_VQACONFIG_H_
#define SRC_VQACONFIG_H_

#include "lattice/lattice.h"
#include <vector>
#include <fstream>
#include <filesystem>

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

		static MatrixInt loadLatticeFromFile(std::string filename, bool *success);

	private:

		std::vector<Lattice> lattices;
		std::string current_hamiltonian = "";

};

#endif /* SRC_VQACONFIG_H_ */
