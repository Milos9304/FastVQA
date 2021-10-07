/*
 * qaoaOptions.cpp
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#include "qaoaOptions.h"
#include "../logger.h"

void QAOAOptions::set_default_stats_function(ExecutionStatistics* executionStatistics, indicators::ProgressBar* bar){

	this->set_default_stats_function(executionStatistics, bar, new HmlLattice(0, 0, ""));
	this->detailedLoggingDisabled = true;
}

void QAOAOptions::set_default_stats_function(ExecutionStatistics* executionStatistics, indicators::ProgressBar* bar, AbstractLatticeInput* lattice){

	this->detailedLoggingDisabled = false;

	stats_function = [this, executionStatistics, bar, lattice](int mode, double energy, double opt_energy, double hit_rate, std::string opt_config) {

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
				  logw("OPT ENERGY " + std::to_string(opt_energy));
				  outfile << energy << " " << opt_energy << " " << hit_rate;

				  if(!detailedLoggingDisabled){
					  for(auto &i : lattice->quboToXvector(opt_config))
						  outfile << " " << i;
				  }
				  outfile << "\n";
				  outfile.flush();

				  break;
			  case 5:
			  default:
				 outfile << energy <<"\n";
				 bar->tick();
				 break;
		  }

	};

	this->logStats = true;
}


