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
/*#include "../fastVQA.h"

#include "../io/indicators/progress_bar.hpp"
#include "../executionStatistics.h"
#include "../io/logger.h"
#include "../io/saveProgress.hpp"
*/

//void run_qaoa(xacc::quantum::PauliOperator, bool verbose);
namespace FastVQA{


class Qaoa{

	public:

		/*static void run_qaoa(ExperimentBuffer* buffer, std::pair<std::vector<double>, std::vector<int>>, std::string name, ExecutionStatistics* execStats, QAOAOptions* options);
		static void run_qaoa(ExperimentBuffer* buffer, std::string hamiltonian, std::string name, ExecutionStatistics* execStats, QAOAOptions* options);
		static void run_qaoa(ExperimentBuffer* buffer, std::string hamiltonian, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, QAOAOptions* options);
		static void run_qaoa(ExperimentBuffer* buffer, std::string hamiltonian, std::pair<std::vector<double>, std::vector<int>> hamiltonian2, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, QAOAOptions* options);
		*/
		//void run_qaoa(ExperimentBuffer* buffer, CostFunction cost_f, int num_qubits, VQEOptions* options);
		void run_qaoa(ExperimentBuffer* buffer, PauliHamiltonian* hamiltonian, QAOAOptions* options);

		//static void run_qaoa_slave_process();

	private:

		int num_qubits;
		std::string instance_name;
		long long int max_iters;

		int log_level;

		void __initialize(ExperimentBuffer* buffer, QAOAOptions* options);
		void __execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* opt);

};

}

#endif /* SRC_QAOA_QAOA_H_ */
