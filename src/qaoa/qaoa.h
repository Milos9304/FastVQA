/*
 * qaoa.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_QAOA_QAOA_H_
#define SRC_QAOA_QAOA_H_

//#include "PauliOperator.hpp"

#include "qaoaOptions.h"
#include "../fastVQA.h"

#include "../io/indicators/progress_bar.hpp"
#include "../executionStatistics.h"
#include "../io/logger.h"
#include "../io/saveProgress.hpp"


//void run_qaoa(xacc::quantum::PauliOperator, bool verbose);

class Qaoa{

	public:

		void run_qaoa(ExperimentBuffer* buffer, Hamiltonian* hamiltonian, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, QAOAOptions* options);

		//static void run_qaoa_slave_process();

	private:

		Ansatz ansatz;

		int num_qubits;
		int num_params;
		int max_iters;

		void execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt, int zero_reference_state, Hamiltonian* hamiltonian);
};



#endif /* SRC_QAOA_QAOA_H_ */
