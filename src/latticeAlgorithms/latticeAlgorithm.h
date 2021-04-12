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

typedef std::function<void()> QOracle;

class LatticeAlgorithm{

	Lattice* lattice;

	public:

		LatticeAlgorithm(Lattice* lattice, QOracle quantum_oracle){
			this->lattice = lattice;
			this->quantum_oracle = quantum_oracle;
		}

		void setLattice(Lattice* lattice){
			this->lattice = lattice;
		}

		void performLLLonLattice();

		virtual void run() = 0;

		virtual ~LatticeAlgorithm() {}

	protected:
		QOracle quantum_oracle;

	private:
		virtual void run_quantum() = 0;

};



#endif /* SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_ */
