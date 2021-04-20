#include "qaoa.h"

#include <iostream>
//#include <xacc.hpp>

//#include "QuEST.h"
#include <random>

#include <fstream>
#include "../logger.h"
#include "../io/saveProgress.hpp"

indicators::ProgressBar* progress_bar;
ExecutionStatistics* execStats;

void run_qaoa(xacc::qbit** buffer,
		std::string hamiltonian,
		std::string name,
		ExecutionStatistics* executionStats,
		QAOAOptions* qaoaOptions){

	run_qaoa(buffer, hamiltonian, name, nullptr, executionStats, qaoaOptions);

}


void run_qaoa(xacc::qbit** buffer,
		std::string hamiltonian,
		std::string name,
		indicators::ProgressBar* bar,
		ExecutionStatistics* executionStats,
		QAOAOptions* qaoaOptions){

   int max_iters = qaoaOptions->max_iters;
   bool verbose = qaoaOptions->verbose;

   if(qaoaOptions->logEnergies)
	   qaoaOptions->outfile.open("statsfile_"+name+".txt", std::fstream::out | std::ios_base::trunc);//std::ios_base::app

   if(bar)
	   progress_bar = bar;

   execStats = executionStats;

   //xacc::setOption("quest-verbose", "true");
   //xacc::setOption("quest-debug", "true");

   //std::cout << hamiltonian <<"\n";throw;

   // The corresponding QUBO Hamiltonian is:

   auto observable = xacc::quantum::getObservable("pauli", hamiltonian);
   logd("Hamiltonian loaded into observable");

   //std::cout << "obs1 " << observable->toString() << "\n";

   int p = qaoaOptions->p;

   int num_qubits = observable->nBits();
   int num_params_per_p = qaoaOptions->extendedParametrizedMode ? observable->getNonIdentitySubTerms().size() + observable->nBits() : 2;
   int num_params_total = p * num_params_per_p;

   //auto buffer = xacc::qalloc(num_qubits);
   *buffer = new xacc::qbit(xacc::qalloc(num_qubits));
   logd("Memory for " + std::to_string(num_qubits) + "qubits allocated");

   if(verbose){
	   logi(qaoaOptions->extendedParametrizedMode ? "parametrized mode" : "normal mode");
	   logi(   std::to_string(observable->nBits()) + " qubits, p="
			   + std::to_string(p) + ", "
			   + std::to_string(num_params_per_p) + " params per p, "
			   + std::to_string(num_params_total) + " params total"
	   );
   }

   std::vector<double> initialParams;
   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<> dis(-2.0, 2.0);


   if(qaoaOptions->loadIntermediate){
	   double expected_energy, sv_energy, hit_rate;
	   bool success = loadProgress(qaoaOptions->l_intermediateName, &initialParams, &expected_energy, &sv_energy, &hit_rate);
	   if(success){
		   if(verbose)
			   logi("Init params loaded from " + qaoaOptions->l_intermediateName);
	   }else{
		   loge("Problem loading params from " + qaoaOptions->l_intermediateName);
		   return;
	   }
   }else{
	   // Init random parameters
	   for(int i = 0; i < num_params_total; ++i){
		 initialParams.emplace_back(dis(gen));
	   }
	   if(verbose)
		   logi("Random params generated");
   }

   auto qaoa = xacc::getService<xacc::Algorithm>("QAOA");

   //std::function<void(int, double)> stats_function = stats_func;

   bool initOk;

   // Doesn't require to prepare the same circuit over and over again, but needs to clone statevect.
   auto acc = xacc::getAccelerator("quest", {{"nbQbits", observable->nBits()}, {"repeated_measurement_strategy", true}});

   auto optimizer = xacc::getOptimizer("nlopt",
		   {{"initial-parameters", initialParams}, {"nlopt-maxeval", max_iters}});



   if(qaoaOptions->isSetLogStats()){
	   initOk = qaoa->initialize({
			   std::make_pair("accelerator", qaoaOptions->accelerator(observable)),
			   std::make_pair("optimizer", qaoaOptions->optimizer(initialParams, max_iters)),
			   std::make_pair("observable", observable),
			   std::make_pair("detailed_log_freq", qaoaOptions->detailed_log_freq),
			   // number of time steps (p) param
			   std::make_pair("steps", qaoaOptions->p),
			   std::make_pair("calc-var-assignment", qaoaOptions->calcVarAssignment),
			   std::make_pair("simplified-simulation", qaoaOptions->simplifiedSimulation),
			   std::make_pair("stats_func", qaoaOptions->get_stats_function()),
			   //Number of samples to estimate optimal variable assignment
			   std::make_pair("parameter-scheme", qaoaOptions->getParameterScheme()),
			   std::make_pair("nbSamples", qaoaOptions->nbSamples_calcVarAssignment),
			   std::make_pair("debugMsgs", qaoaOptions->debug)
		});
   }
   else{
	   initOk = qaoa->initialize({
		   	   std::make_pair("accelerator", qaoaOptions->accelerator(observable)),
		   	   std::make_pair("optimizer", qaoaOptions->optimizer(initialParams, max_iters)),
	   		   std::make_pair("observable", /*static_cast<xacc::Observable*>(&*/observable/*)*/),
			   std::make_pair("detailed_log_freq", qaoaOptions->detailed_log_freq),
	   		   // number of time steps (p) param
	   		   std::make_pair("steps", qaoaOptions->p),
	   		   std::make_pair("calc-var-assignment", qaoaOptions->calcVarAssignment),
	   		   std::make_pair("simplified-simulation", qaoaOptions->simplifiedSimulation),
	   		   //Number of samples to estimate optimal variable assignment
			   std::make_pair("parameter-scheme", qaoaOptions->getParameterScheme()),
	   		   std::make_pair("nbSamples", qaoaOptions->nbSamples_calcVarAssignment),
			   std::make_pair("debugMsgs", qaoaOptions->debug)
	   		});
   }
   if(initOk)
	   logi("QAOA init successful.");
   else{
   	   loge("QAOA Init error!");
   	   return;
   }

   logd("Executing qaoa");
   qaoa->execute(**buffer);
   logd("Qaoa execution done");

   if(qaoaOptions->saveIntermediate){
	   std::vector<double> params = (**buffer)["opt-params"].as<std::vector<double>>();
	   double expected_energy = (**buffer)["expected-val"].as<double>();
	   double sv_energy = (**buffer)["opt-val"].as<double>();
	   double hit_rate = (**buffer)["hit_rate"].as<double>();

	   saveProgress(qaoaOptions->s_intermediateName, params, expected_energy, sv_energy, hit_rate);
   }


   //std::cout << "Min QUBO: " << (**buffer)["opt-val"].as<double>() << "\n";
   //std::vector<double> params = (*buffer)["opt-params"].as<std::vector<double>>();



   qaoaOptions->outfile.close();

}
