/*
 * qaoaOptions.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_QAOA_QAOAOPTIONS_H_
#define SRC_QAOA_QAOAOPTIONS_H_

#include "../indicators/progress_bar.hpp"
#include "../executionStatistics.h"
#include <functional>
#include <fstream>
#include <xacc.hpp>

typedef std::function<std::shared_ptr<xacc::Accelerator>(std::shared_ptr<xacc::Observable>)> AcceleratorPartial;
typedef std::function<std::shared_ptr<xacc::Optimizer>(std::vector<double>, int)> OptimizerPartial;

class QAOAOptions{

	public:

		int max_iters = 0;
		int p=1;
		int nbSamples_calcVarAssignment = 1024;

		bool calcVarAssignment = false;
		bool extendedParametrizedMode = false;
		bool simplifiedSimulation = true;

		bool saveIntermediate = false;
		std::string s_intermediateName="";
		bool loadIntermediate = false;
		std::string l_intermediateName="";
		bool logEnergies = false;

		AcceleratorPartial accelerator;
		OptimizerPartial optimizer;

		bool isSetLogStats(){
			return logStats;
		}

		std::string getParameterScheme(){
			return extendedParametrizedMode ? "Extended" : "Standard";
		}

		bool verbose = true;
		std::fstream outfile;

		std::function<void(int, double)> get_stats_function(){
			return stats_function;
		}


		void set_default_stats_function(ExecutionStatistics* execStats, indicators::ProgressBar* bar);

	private:
		std::function<void(int, double)> stats_function;
		bool logStats = false;

};

#endif /* SRC_QAOA_QAOAOPTIONS_H_ */
