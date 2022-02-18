/*
 * vqe.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef VQAS_VQE_H_
#define VQAS_VQE_H_

/*#include "../lattice/lattice.h"
#include "../io/logger.h"
#include "../io/saveProgress.hpp"
#include "../io/indicators/progress_bar.hpp"*/
#include "ansatz.h"
/*#include "../optimizer/optimizer.h"*/
#include "vqeOptions.h"

namespace fastVQA{
class Vqe{

	public:

		void run_vqe(ExperimentBuffer* buffer, Hamiltonian* hamiltonian, std::string name, VQEOptions* options);

	private:

		Ansatz ansatz;

		int num_qubits;
		int num_params;
		int max_iters;
		int nbSamples_calcVarAssignment;

		void execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt, std::vector<long long unsigned int> zero_reference_states, Hamiltonian* hamiltonian);
};
}

#endif /* VQAS_VQE_H_ */
