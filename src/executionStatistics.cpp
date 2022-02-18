/*
 * executionStatistics.cpp
 *
 *  Created on: Apr 14, 2021
 *      Author: Milos Prokop
 */

#include "executionStatistics.h"

void ExecutionStatistics::startQuantumIterLog(){
	quantum_time_start = std::chrono::steady_clock::now();
}

void ExecutionStatistics::finishQuantumIterLog(){
	double count = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-quantum_time_start).count();
	quantum_iteration_time.update(count);
}

void ExecutionStatistics::startOptimizerIterLog(){
	optimizer_time_start = std::chrono::steady_clock::now();
}

void ExecutionStatistics::finishOptimizerIterLog(){
	double count = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-optimizer_time_start).count();
	optimizer_iteration_time.update(count);
}


