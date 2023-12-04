/*
 * experimentBuffer.h
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_EXPERIMENT_BUFFER_H_
#define FASTVQA_EXPERIMENT_BUFFER_H_

#include "QuEST.h"
#include <memory>

namespace FastVQA{

class ExperimentBuffer{
public:
	class ExperimentBufferSolution{
	public:
		ExperimentBufferSolution(std::string opt_config, double expected_energy, double hit_rate){
			this->opt_config=opt_config;
			this->expected_energy=expected_energy;
			this->hit_rate=hit_rate;
		}
		std::vector<double> opt_params;
		std::string opt_config;
		double expected_energy;
		double hit_rate;
	};
	double opt_val;
	std::vector<ExperimentBufferSolution> final_solutions;
	double getTotalHitRate(){
		double hit_rate=0;
		for(auto &fs : this->final_solutions){
			hit_rate+=fs.hit_rate;
		}
		return hit_rate;
	}

	std::vector<std::pair<std::string, double>> initial_params; //for debug purposes

	/*
	* All the return structures below are optional,
	* i.e. they will be written only if enabled by a flag
	*/

	// Reference to final stateVector
	bool storeQuregPtr = false;
	std::shared_ptr<Qureg> stateVector;

	std::vector<double> intermediateEnergies;
	std::vector<double> intermediateGroundStateOverlaps;
};
}
#endif /* FASTVQA_EXPERIMENT_BUFFER_H_ */
