#include "vqe.h"

#include <iostream>
//#include <xacc.hpp>

//#include "QuEST.h"
#include <random>

#include <fstream>
#include "../logger.h"
#include "../io/saveProgress.hpp"
#include "mpi.h"

void Vqe::run_vqe(xacc::qbit** buffer,
		std::string hamiltonian,
		std::string name,
		ExecutionStatistics* executionStats,
		VQEOptions* vqeOptions){

	Vqe::_run_vqe(buffer, hamiltonian, std::pair<std::vector<double>, std::vector<int>>(), name, nullptr, executionStats, vqeOptions);

}

void Vqe::run_vqe(xacc::qbit** buffer,
		std::pair<std::vector<double>, std::vector<int>> hamiltonian,
		std::string name,
		ExecutionStatistics* executionStats,
		VQEOptions* vqeOptions){

	Vqe::_run_vqe(buffer, "", hamiltonian, name, nullptr, executionStats, vqeOptions);

}


void Vqe::run_vqe(xacc::qbit** buffer,
		std::string hamiltonian,
		std::string name,
		indicators::ProgressBar* bar,
		ExecutionStatistics* executionStats,
		VQEOptions* vqeOptions){

	Vqe::_run_vqe(buffer, hamiltonian, std::pair<std::vector<double>, std::vector<int>>(), name, bar, executionStats, vqeOptions);

}

void Vqe::run_vqe(xacc::qbit** buffer,
		std::string hamiltonian,
		std::pair<std::vector<double>, std::vector<int>> hamiltonian2,
		std::string name,
		indicators::ProgressBar* bar,
		ExecutionStatistics* executionStats,
		VQEOptions* vqeOptions){

	Vqe::_run_vqe(buffer, hamiltonian, hamiltonian2, name, bar, executionStats, vqeOptions);

}

void Vqe::_run_vqe(xacc::qbit** buffer,
		std::string hamiltonian,
		std::pair<std::vector<double>, std::vector<int>> hamiltonian2,
		std::string name,
		indicators::ProgressBar* bar,
		ExecutionStatistics* executionStats,
		VQEOptions* vqeOptions){

   int max_iters = vqeOptions->max_iters;
   bool verbose = vqeOptions->verbose;

   if(vqeOptions->logEnergies){
	   //vqeOptions->outfile.rdbuf()->pubsetbuf(0, 0); //disable buffer
	   vqeOptions->outfile.open("../experiment_files/statsfile_"+name+".txt", std::fstream::out | std::ios_base::trunc); //| std::ios_base::trunc);//std::ios_base::app
   }

   if(bar)
	   vqeOptions->progress_bar = bar;

   vqeOptions->execStats = executionStats;

   //xacc::setOption("quest-verbose", "true");
   //xacc::setOption("quest-debug", "true");

   //std::cout << hamiltonian <<"\n";throw;

   // The corresponding QUBO Hamiltonian is:

   auto observable = xacc::quantum::getObservable("pauli", hamiltonian);
   logd("Hamiltonian loaded into observable");

   int num_qubits = observable->nBits();
   //int num_params_per_p = vqeOptions->extendedParametrizedMode ? observable->getNonIdentitySubTerms().size() + observable->nBits() : 2;
   int num_params_total = 2;//p * num_params_per_p;

   //auto buffer = xacc::qalloc(num_qubits);
   *buffer = new xacc::qbit(xacc::qalloc(num_qubits));
   logd("Memory for " + std::to_string(num_qubits) + "qubits allocated");

   if(verbose){
	   //logi(vqeOptions->extendedParametrizedMode ? "parametrized mode" : "normal mode");
	   logi(std::to_string(observable->nBits()) + " qubits");
   }

   std::vector<double> initialParams;
   //std::random_device rd;
   std::mt19937 gen(1997); //rd() instead of 1997
   std::uniform_real_distribution<> dis(-2.0, 2.0);

   if(vqeOptions->loadIntermediate){
	   double expected_energy, sv_energy, hit_rate;
	   bool success = loadProgress(vqeOptions->l_intermediateName, &initialParams, &expected_energy, &sv_energy, &hit_rate);
	   if(success){
		   if(verbose)
			   logi("Init params loaded from " + vqeOptions->l_intermediateName);
	   }else{
		   loge("Problem loading params from " + vqeOptions->l_intermediateName);
		   return;
	   }
   }else{
	   // Init random parameters
	   for(int i = 0; i < num_params_total; ++i){
		 //initialParams.emplace_back(dis(gen));
		   double p=dis(gen);logd(std::to_string(p));
		   initialParams.emplace_back(p);
	   }
	   if(verbose)
		   logi("Random params generated");
   }

   auto vqe = xacc::getService<xacc::Algorithm>("vqe");

   //std::function<void(int, double)> stats_function = stats_func;

   bool initOk;
   std::vector<double> hamCoeffs;
   std::vector<int> hamPauliCodes;

   if(vqeOptions->provideHamiltonian){
	   hamCoeffs = hamiltonian2.first;
	   hamPauliCodes = hamiltonian2.second;
   }

   // Doesn't require to prepare the same circuit over and over again, but needs to clone statevect.
   //auto acc = xacc::getAccelerator("quest",
   //	   {{"nbQbits", observable->nBits()},
   //		{"startWithPlusState", true},
   //		{"repeated_measurement_strategy", true},
   //		{"hamiltonianProvided", vqeOptions->provideHamiltonian},
   //	    {"hamiltonianCoeffs", hamCoeffs},
   //		{"pauliCodes", hamPauliCodes}
   //	   });

   auto optimizer = xacc::getOptimizer("nlopt",
		   {{"initial-parameters", initialParams}, {"nlopt-maxeval", max_iters}});

   //if(vqeOptions->isSetLogStats()){
	   /*initOk = qaoa->initialize({
		   	   std::make_pair("name", name),
			   std::make_pair("accelerator", vqeOptions->accelerator(observable,
					   vqeOptions->provideHamiltonian,
					   hamCoeffs,
					   hamPauliCodes,
					   name)),
			   std::make_pair("optimizer", vqeOptions->optimizer(initialParams, max_iters)),
			   std::make_pair("observable", observable),
			   std::make_pair("detailed_log_freq", vqeOptions->detailed_log_freq),
			   // number of time steps (p) param
			   std::make_pair("steps", vqeOptions->p),
			   std::make_pair("load_ansatz",  vqeOptions->load_ansatz),
			   std::make_pair("calc-var-assignment", vqeOptions->calcVarAssignment),
			   std::make_pair("simplified-simulation", vqeOptions->simplifiedSimulation),
			   std::make_pair("stats_func", vqeOptions->get_stats_function()),
			   //Number of samples to estimate optimal variable assignment
			   std::make_pair("parameter-scheme", vqeOptions->getParameterScheme()),
			   std::make_pair("nbSamples", vqeOptions->nbSamples_calcVarAssignment),
			   std::make_pair("questHamExpectation", vqeOptions->provideHamiltonian),
			   std::make_pair("debugMsgs", vqeOptions->debug),
			   std::make_pair("overlapTrick", vqeOptions->overlap_trick),
			   std::make_pair("zeroRefState", vqeOptions->zero_reference_state),
		});*/

   	   //Ansatz construction
   	   	 xacc::HeterogeneousMap m;
   	   	 std::shared_ptr<xacc::CompositeInstruction> ansatz = std::dynamic_pointer_cast<xacc::CompositeInstruction>(
             xacc::getService<xacc::Instruction>("qaoa"));
         m.insert("nbQubits", num_qubits);
         m.insert("nbSteps", /*m_nbSteps*/1);
         m.insert("ref-ham", nullptr);
         m.insert("cost-ham", observable);
         m.insert("parameter-scheme", "Standard");
         /*if (m_initial_state){
             m.insert("initial-state", m_initial_state);
         }*/
         m.insert("shuffle-terms", /*m_shuffleTerms*/false);
         ansatz->expand(m);

         if(vqeOptions->logEnergies)
        	 initOk=vqe->initialize({std::make_pair("ansatz", ansatz),
	                    std::make_pair("observable", observable),
						std::make_pair("stats_func", vqeOptions->get_stats_function()),
	                    std::make_pair("accelerator", vqeOptions->accelerator(observable,
	     					   vqeOptions->provideHamiltonian,
	     					   hamCoeffs,
	     					   hamPauliCodes,
	     					   name)),
	                    std::make_pair("optimizer", optimizer),
						std::make_pair("overlapTrick", vqeOptions->overlap_trick),
						std::make_pair("zeroRefState", vqeOptions->zero_reference_state)});
         else
        	 initOk=vqe->initialize({std::make_pair("ansatz", ansatz),
        	 	                    std::make_pair("observable", observable),
        	 	                    std::make_pair("accelerator", vqeOptions->accelerator(observable,
        	 	     					   vqeOptions->provideHamiltonian,
        	 	     					   hamCoeffs,
        	 	     					   hamPauliCodes,
        	 	     					   name)),
        	 	                    std::make_pair("optimizer", optimizer),
									std::make_pair("overlapTrick", vqeOptions->overlap_trick),
									std::make_pair("zeroRefState", vqeOptions->zero_reference_state)});
   //}
	   /*
   else{
	   initOk = qaoa->initialize({
		   	   std::make_pair("name", name),
		       std::make_pair("accelerator", vqeOptions->accelerator(observable,
		    		   vqeOptions->provideHamiltonian,
					   hamCoeffs,
					   hamPauliCodes,
					   name)),
		   	   std::make_pair("optimizer", vqeOptions->optimizer(initialParams, max_iters)),


	   		   //std::make_pair("observable", static_cast<xacc::Observable*>(&observable)),
			   std::make_pair("observable", observable),

			   std::make_pair("detailed_log_freq", vqeOptions->detailed_log_freq),
	   		   // number of time steps (p) param
	   		   std::make_pair("steps", vqeOptions->p),
			   std::make_pair("load_ansatz",  vqeOptions->load_ansatz),
	   		   std::make_pair("calc-var-assignment", vqeOptions->calcVarAssignment || vqeOptions->load_ansatz),
	   		   std::make_pair("simplified-simulation", vqeOptions->simplifiedSimulation),
	   		   //Number of samples to estimate optimal variable assignment
			   std::make_pair("parameter-scheme", vqeOptions->getParameterScheme()),
	   		   std::make_pair("nbSamples", vqeOptions->nbSamples_calcVarAssignment),
			   std::make_pair("questHamExpectation", vqeOptions->provideHamiltonian),
			   std::make_pair("debugMsgs", vqeOptions->debug),
			   std::make_pair("overlapTrick", vqeOptions->overlap_trick),
			   std::make_pair("zeroRefState", vqeOptions->zero_reference_state)
	   		});
   }*/

   if(initOk)
	   logi("VQE init successful.");
   else{
   	   loge("VQE Init error!");
   	   return;
   }

   logd("Executing vqe");
   vqe->execute(**buffer);
   logd("Vqe execution done");

   if(vqeOptions->saveIntermediate){
	   std::vector<double> params = (**buffer)["opt-params"].as<std::vector<double>>();
	   double expected_energy = (**buffer)["expected-val"].as<double>();
	   double sv_energy = (**buffer)["opt-val"].as<double>();
	   double hit_rate = (**buffer)["hit_rate"].as<double>();

	   saveProgress(vqeOptions->s_intermediateName, params, expected_energy, sv_energy, hit_rate);
   }


   std::cout << "Min QUBO: " << (**buffer)["opt-val"].as<double>() << "\n";
   //std::vector<double> params = (*buffer)["opt-params"].as<std::vector<double>>();

   vqeOptions->outfile.close();

}

void Vqe::run_vqe_slave_process(){
	auto acc = xacc::getAccelerator("quest");
}
