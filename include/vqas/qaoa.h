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

		static void run_qaoa(ExperimentBuffer* buffer, std::pair<std::vector<double>, std::vector<int>>, std::string name, ExecutionStatistics* execStats, QAOAOptions* options);
		static void run_qaoa(ExperimentBuffer* buffer, std::string hamiltonian, std::string name, ExecutionStatistics* execStats, QAOAOptions* options);
		static void run_qaoa(ExperimentBuffer* buffer, std::string hamiltonian, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, QAOAOptions* options);
		static void run_qaoa(ExperimentBuffer* buffer, std::string hamiltonian, std::pair<std::vector<double>, std::vector<int>> hamiltonian2, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, QAOAOptions* options);

		static void run_qaoa_slave_process();

	private:

		static void _run_qaoa(ExperimentBuffer* buffer, std::string hamiltonian, std::pair<std::vector<double>, std::vector<int>> hamiltonian2, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, QAOAOptions* options);

};



#endif /* SRC_QAOA_QAOA_H_ */
