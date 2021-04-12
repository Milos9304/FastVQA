/*
 * latticeAlgorithm.h
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_
#define SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_

#include "../lattice.h"

class LatticeAlgorithm{

	Lattice* lattice;

	public:

		void setLattice(Lattice* lattice){
					this->lattice = lattice;
		}

		void performLLLonLattice();

		//virtual void run() = 0;

		virtual ~LatticeAlgorithm() {}

};



#endif /* SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_ */
