#include "xacc.hpp"

#include "qaoa.h"

#include <iostream>
#include <xacc.hpp>
#include "xacc_observable.hpp"
#include "xacc_service.hpp"
//#include "QuEST.h"
#include <random>

#include <fstream>
#include "../logger.h"

const std::string accelerator_name = "quest";
std::ofstream outfile;

void stats_func(double energy){

  outfile << energy<<"\n";

}

void run_qaoa(VqaConfig* vqaConfig){

   outfile.open("statsfile.txt", std::fstream::out | std::ios_base::trunc);//std::ios_base::app
   xacc::Initialize();

   //xacc::setOption("quest-verbose", "true");
   //xacc::setOption("quest-debug", "true");

   // The corresponding QUBO Hamiltonian is:
   auto observable = xacc::quantum::getObservable("pauli", vqaConfig->getCurrentHamiltonian());

   int p = 1;

   int num_qubits = observable->nBits();
   int num_params_per_p = observable -> getNonIdentitySubTerms().size() + observable->nBits();
   int num_params_total = p * num_params_per_p;

   auto buffer = xacc::qalloc(num_qubits);

   if(vqaConfig->verbose){
	   logi(std::to_string(observable->nBits()) + " qubits, p="
			   + std::to_string(p) + ", "
			   + std::to_string(num_params_per_p) + " params per p"
			   + std::to_string(num_params_total) + " params total"
	   );
   }

   auto acc = xacc::getAccelerator(accelerator_name, {std::make_pair("nbQbits", observable->nBits()),
	         // Doesn't require to prepare the same circuit over and over again, but needs to clone statevect.
		   	   	   	   	   	   	   	   	   	   	   	  std::make_pair("repeated_measurement_strategy", true)});

   std::vector<double> initialParams;
   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<> dis(-2.0, 2.0);

   // Init random parameters
   for(int i = 0; i < num_params_total; ++i){
      initialParams.emplace_back(dis(gen));
   }

   auto optimizer = xacc::getOptimizer("nlopt",
      xacc::HeterogeneousMap {
         std::make_pair("initial-parameters", initialParams),
         std::make_pair("nlopt-maxeval", 0)});//num_params_total*1/*100*/) });

   auto qaoa = xacc::getService<xacc::Algorithm>("QAOA");

   std::function<void(double)> stats_function = stats_func;

   const bool initOk = qaoa->initialize({
	   std::make_pair("accelerator", acc),
	   std::make_pair("optimizer", optimizer),
	   std::make_pair("observable", observable),
	   // number of time steps (p) param
	   std::make_pair("steps", p),
	   std::make_pair("calc-var-assignment", true),
	   std::make_pair("simplified-simulation", true),
	   std::make_pair("stats_func", stats_function),
	   //Number of samples to estimate optimal variable assignment
	   std::make_pair("nbSamples", 1024/*5*/)
    });

   if(initOk)
	   logi("QAOA init successful.");
   else{
   	   loge("QAOA Init error!");
   	   return;
   }

   qaoa->execute(buffer);

   std::cout << "Min QUBO: " << (*buffer)["opt-val"].as<double>() << "\n";
   std::vector<double> params = (*buffer)["opt-params"].as<std::vector<double>>();


   xacc::Finalize();

   outfile.close();
}
