/*
 * vqaOptions.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_VQAOPTIONS_H_
#define SRC_VQAOPTIONS_H_

#include "io/logger.h"
#include "io/indicators/progress_bar.hpp"
#include "executionStatistics.h"
//#include "lattice/hmlLattice.hpp"
#include "lattice/lattice.h"
#include "accelerator/accelerator.h"
#include "optimizer/optimizer.h"
#include <functional>
#include <fstream>

typedef std::function<void(int, double, double, double, std::string)> StatsFunction;

typedef std::function<std::shared_ptr<Optimizer>(std::vector<double>, int)> OptimizerPartial;

class VQAOptions{

	public:

		indicators::ProgressBar* progress_bar;
		ExecutionStatistics* execStats;

		int max_iters = 0;
		//int p=1;
		int nbSamples_calcVarAssignment = 1024;
		int detailed_log_freq = 0;

		bool calcVarAssignment = false;
		//bool extendedParametrizedMode = false;
		bool simplifiedSimulation = true;

		std::string ansatz_name = "Ry_CNOT_all2all_Ry";

		bool saveIntermediate = false;
		std::string s_intermediateName="";
		bool loadIntermediate = false;
		std::string l_intermediateName="";
		bool logEnergies = false;

		bool save_ansatz = false;
		bool load_ansatz = false;

		bool debug = false; //print debug msgs
		bool provideHamiltonian = false;

		bool overlap_trick = false; //penalization strategy
		int zero_reference_state=0;

		Accelerator* accelerator;
		Optimizer* optimizer;

		bool isSetLogStats(){
			return logStats;
		}

		bool verbose = true;
		std::fstream outfile;

		StatsFunction get_stats_function(){
			return stats_function;
		}

		//void set_default_stats_function(ExecutionStatistics* execStats, indicators::ProgressBar* bar);
		void set_default_stats_function(ExecutionStatistics* execStats, indicators::ProgressBar* bar, Lattice* lattice);

	private:
		std::function<void(int, double, double, double, std::string)> stats_function;
		bool logStats = false;
		bool detailedLoggingDisabled = true;

};

#endif /* SRC_VQAOPTIONS_H_ */
