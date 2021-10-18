/*
 * vqe.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_VQE_VQE_H_
#define SRC_VQE_VQE_H_

//#include "PauliOperator.hpp"
#include "../indicators/progress_bar.hpp"
#include "../executionStatistics.h"
#include "xacc_observable.hpp"
#include "xacc_service.hpp"
#include "xacc.hpp"
#include "vqeOptions.h"

//void run_qaoa(xacc::quantum::PauliOperator, bool verbose);

class Vqe{

	public:

		static void run_vqe(xacc::qbit** buffer, std::pair<std::vector<double>, std::vector<int>>, std::string name, ExecutionStatistics* execStats, VQEOptions* options);
		static void run_vqe(xacc::qbit** buffer, std::string hamiltonian, std::string name, ExecutionStatistics* execStats, VQEOptions* options);
		static void run_vqe(xacc::qbit** buffer, std::string hamiltonian, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, VQEOptions* options);
		static void run_vqe(xacc::qbit** buffer, std::string hamiltonian, std::pair<std::vector<double>, std::vector<int>> hamiltonian2, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, VQEOptions* options);

		static void run_vqe_slave_process();

	private:

		static void _run_vqe(xacc::qbit** buffer, std::string hamiltonian, std::pair<std::vector<double>, std::vector<int>> hamiltonian2, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, VQEOptions* options);

};



#endif /* SRC_QAOA_QAOA_H_ */
