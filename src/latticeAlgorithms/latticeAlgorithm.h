/*
 * latticeAlgorithm.h
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_
#define SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_

#include "../lattice.h"
#include <functional>
#include "PauliOperator.hpp"

typedef std::function<void(std::string hamiltonian/*xacc::quantum::PauliOperator*/, std::string name)> QOracle;

class LatticeAlgorithm{

	protected:
		Lattice* lattice;
		MapOptions* options;
		xacc::quantum::PauliOperator hamiltonian;
		QOracle quantum_oracle;

	public:

		LatticeAlgorithm(Lattice* lattice, MapOptions* options, QOracle quantum_oracle){
			this->lattice = lattice;
			this->options = options;
			this->quantum_oracle = quantum_oracle;
			this->hamiltonian = lattice->getHamiltonian(options);
		}

		void setLattice(Lattice* lattice){
			this->lattice = lattice;
		}

		void performLLLonLattice();

		virtual void run() = 0;

		virtual ~LatticeAlgorithm() {}

	private:
		virtual void run_quantum() = 0;

};



#endif /* SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_ */
