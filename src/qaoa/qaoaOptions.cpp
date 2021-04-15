/*
 * qaoaOptions.cpp
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#include "qaoaOptions.h"

void QAOAOptions::set_default_stats_function(ExecutionStatistics* executionStatistics, indicators::ProgressBar* bar){

	stats_function = [this, executionStatistics, bar](int mode, double energy) {

			switch(mode){
			  case 0:
				  executionStatistics->startQuantumIterLog(); //quantum_init
				  break;
			  case 1:
				  executionStatistics->finishQuantumIterLog();
				  break;
			  case 2:
				  executionStatistics->startOptimizerIterLog();
				  break;
			  case 3:
				  executionStatistics->finishOptimizerIterLog();
				  break;
			  case 4:
			  default:
				 outfile << energy <<"\n";
				 bar->tick();
				 break;
		  }

	};

	this->logStats = true;
}


