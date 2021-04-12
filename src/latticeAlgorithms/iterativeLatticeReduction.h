/*
 * iterativeLatticeReduction.h
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICEALGORITHMS_ITERATIVELATTICEREDUCTION_H_
#define SRC_LATTICEALGORITHMS_ITERATIVELATTICEREDUCTION_H_

#include "latticeAlgorithm.h"

class IterativeLatticeReduction : public LatticeAlgorithm {

	int n_iters;

	public:

		IterativeLatticeReduction(Lattice* lattice, QOracle quantum_oracle, int n_iters) :
			LatticeAlgorithm(lattice,  quantum_oracle) {

				this->n_iters = n_iters;

		}
		void run();

	private:

		void run_quantum();

};



#endif /* SRC_LATTICEALGORITHMS_ITERATIVELATTICEREDUCTION_H_ */
