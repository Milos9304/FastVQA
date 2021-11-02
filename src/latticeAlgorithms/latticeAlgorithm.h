/*
 * latticeAlgorithm.h
 *
 *  Created on: Apr 12, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_
#define SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_

#include "../fastVQA.h"
#include "../io/logger.h"
#include "../lattice/lattice.h"
#include <functional>
#include "fplll.h"

//#include "PauliOperator.hpp"
//#include "xacc.hpp"

typedef std::function<void(ExperimentBuffer* buffer, std::string name)> QOracle;

class LatticeAlgorithm{

	protected:
		Lattice* lattice;
		MapOptions* options;
		//xacc::quantum::PauliOperator hamiltonian;
		QOracle quantum_oracle;

	public:

		LatticeAlgorithm(Lattice* lattice, MapOptions* options, QOracle quantum_oracle){
			this->lattice = lattice;
			this->options = options;
			this->quantum_oracle = quantum_oracle;
			//this->hamiltonian = lattice->getHamiltonian(options);
		}

		void setLattice(Lattice* lattice){
			this->lattice = lattice;
		}

		void performLLLonLattice(double delta=0.99, double eta=0.51, LLLMethod method=LLLMethod::LM_PROVED);

		VectorInt xVectToShortVect(VectorInt* x_vect);

		virtual void run() = 0;
		virtual std::pair<VectorInt, double> run_test() = 0;

		virtual ~LatticeAlgorithm() {}

	private:
		virtual std::pair<std::string, double> run_quantum() = 0; //opt config and its energy

};



#endif /* SRC_LATTICEALGORITHMS_LATTICEALGORITHM_H_ */
