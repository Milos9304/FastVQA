/*
 * stats.h
 *
 *  Created on: Apr 14, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_EXECUTIONSTATISTICS_H_
#define SRC_EXECUTIONSTATISTICS_H_

#include "utils.hpp"
#include "gmpxx.h"
#include <chrono>

class ExecutionStatistics{

	int n_quantum_oracle_calls;

	std::chrono::steady_clock::time_point quantum_time_start, optimizer_time_start;

public:

	OnlineMeanVarianceCalculator quantum_iteration_time;
	OnlineMeanVarianceCalculator optimizer_iteration_time;
	//OnlineMeanVarianceCalculator expected_energy_stats;

	void startQuantumIterLog();
	void finishQuantumIterLog();

	void startOptimizerIterLog();
	void finishOptimizerIterLog();

	void reset(){
		quantum_iteration_time.reset();
		optimizer_iteration_time.reset();
	}

};



#endif /* SRC_EXECUTIONSTATISTICS_H_ */
