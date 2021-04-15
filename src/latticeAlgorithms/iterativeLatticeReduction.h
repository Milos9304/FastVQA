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
	bool testing;

	public:

		IterativeLatticeReduction(Lattice* lattice, MapOptions* options, QOracle quantum_oracle, int n_iters, bool testing=false) : // @suppress("Class members should be properly initialized")
			LatticeAlgorithm(lattice, options, quantum_oracle) {

				this->n_iters = n_iters;
				this->testing = testing;

		}
		void run();
		void run_test();

	private:

		std::pair<std::string, double> run_quantum(); //returns optimal configuration and its energy

};



#endif /* SRC_LATTICEALGORITHMS_ITERATIVELATTICEREDUCTION_H_ */
