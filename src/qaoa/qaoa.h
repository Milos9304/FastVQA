/*
 * qaoa.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_QAOA_QAOA_H_
#define SRC_QAOA_QAOA_H_

//#include "PauliOperator.hpp"
#include "../indicators/progress_bar.hpp"
#include "../executionStatistics.h"
#include "xacc_observable.hpp"
#include "xacc_service.hpp"
#include "xacc.hpp"
#include "qaoaOptions.h"


//void run_qaoa(xacc::quantum::PauliOperator, bool verbose);



void run_qaoa(xacc::qbit** buffer, std::string hamiltonian, std::string name, ExecutionStatistics* execStats, QAOAOptions* options);
void run_qaoa(xacc::qbit** buffer, std::string hamiltonian, std::string name, indicators::ProgressBar* bar, ExecutionStatistics* execStats, QAOAOptions* options);

#endif /* SRC_QAOA_QAOA_H_ */
