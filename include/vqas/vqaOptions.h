/*
 * vqaOptions.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_VQAOPTIONS_H_
#define FASTVQA_VQAOPTIONS_H_

//#include "logger.h"
//#include "indicators/progress_bar.hpp"
/*#include "executionStatistics.h"
//#include "lattice/hmlLattice.hpp"
#include "lattice/lattice.h"
#include "accelerator/accelerator.h"
#include "optimizer/optimizer.h"
#include <functional>
#include <fstream>*/

//typedef std::function<void(int, double, double, double, std::string)> StatsFunction;
//typedef std::function<std::shared_ptr<Optimizer>(std::vector<double>, int)> OptimizerPartial;

#include "accelerator.h"
#include "optimizer.h"

namespace FastVQA{

class VQAOptions{

	public:

		//0 - debug, 1 - info, 2 - warning, 3 - error
		int log_level = 1;

		//indicators::ProgressBar* progress_bar;
		//ExecutionStatistics* execStats;

		// Maximum number of iterations. 0 for unlimited
		long long int max_iters = 0;

		// Maximum number of iterations. 0 for unlimited
		double ftol = 10e-6;

		// If a non-zero integer x is set, calculate the final assignment every x iterations.
		int detailed_log_freq = 0;

		// Log the intermediate energies to file
		bool logEnergies = false;

		// Number of measurement samples of final ansatz
		int nbSamples_calcVarAssignment = 1024;

		// Print calculated expectation after each cycle to the standard output
		bool expectationToStandardOutput = false;

		// Instance name
		std::string instance_name;

		//bool calcVarAssignment = false;
		//bool extendedParametrizedMode = false;
		//bool simplifiedSimulation = true;

		// This does not free qureg from memory after run of VQA.
		// Can be used for post run analysis. Need to destroy the qureg manually
		// before a new experiment is about to be run.
		bool keepReferenceToQureg = false;

		//bool saveIntermediate = false;
		//std::string s_intermediateName="";
		//bool loadIntermediate = false;
		//std::string l_intermediateName="";


		//bool save_ansatz = false;
		//bool load_ansatz = false;


		//bool provideHamiltonian = false;

		//penalization strategy
		//TODO: more info
	    bool overlap_trick = false;

	    //Integer representations of variable assignment corresponding to zero state
		std::vector<long long unsigned int> zero_reference_states;

		Accelerator* accelerator;
		Optimizer* optimizer;

		//bool isSetLogStats(){
		//	return logStats;
		//}


		//std::fstream outfile;

		//StatsFunction get_stats_function(){
		//	return stats_function;
		//}

		//void set_default_stats_function(ExecutionStatistics* execStats, indicators::ProgressBar* bar);
		//void set_default_stats_function(ExecutionStatistics* execStats, indicators::ProgressBar* bar, Lattice* lattice);

	private:
		std::function<void(int, double, double, double, std::string)> stats_function;
		bool logStats = false;
		bool detailedLoggingDisabled = true;

};
}
#endif /* FASTVQA_VQAOPTIONS_H_ */
