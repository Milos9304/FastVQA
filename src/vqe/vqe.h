/*
 * vqe.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_VQE_VQE_H_
#define SRC_VQE_VQE_H_

//#include "PauliOperator.hpp"
#include "../lattice/lattice.h"
#include "../io/logger.h"
#include "../io/saveProgress.hpp"
#include "../io/indicators/progress_bar.hpp"
#include "../executionStatistics.h"
#include "../fastVQA.h"
#include "../ansatz/ansatz.h"
#include "../optimizer/optimizer.h"
#include "vqeOptions.h"

//void run_qaoa(xacc::quantum::PauliOperator, bool verbose);

class Vqe{

	public:

		void run_vqe(ExperimentBuffer* buffer, Hamiltonian* hamiltonian, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, VQEOptions* options, int seed=1997);

	private:

		Ansatz ansatz;

		int num_qubits;
		int num_params;
		int max_iters;
		/*bool initialize(std::string name, Ansatz* ansatz,
			Hamiltonian* observable,
			StatsFunction stats_function,
			Optimizer* optimizer,
			bool overlap_trick,
			int zero_reference_state);
*/
		void execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt, int zero_reference_state, Hamiltonian* hamiltonian);
};

#endif /* SRC_QAOA_QAOA_H_ */
