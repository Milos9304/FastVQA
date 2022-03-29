#include "vqas/vqe.h"
#include "logger.h"

#include <iostream>
//#include <xacc.hpp>

//#include "QuEST.h"
#include <random>

#include <fstream>


//#include "mpi.h"
namespace fastVQA{
void Vqe::run_vqe(ExperimentBuffer* buffer,
		Hamiltonian* hamiltonian,
		VQEOptions* options){

   max_iters = options->max_iters;
   instance_name = options->instance_name;
   nbSamples_calcVarAssignment = options->nbSamples_calcVarAssignment;
   int log_level = options->log_level;

   //if(vqeOptions->logEnergies){
	   //vqeOptions->outfile.rdbuf()->pubsetbuf(0, 0); //disable buffer
	//   vqeOptions->outfile.open("../experiment_files/statsfile_"+name+".txt", std::fstream::out | std::ios_base::trunc); //| std::ios_base::trunc);//std::ios_base::app
   //}

   //vqeOptions->execStats = executionStats;

   num_qubits = hamiltonian->nbQubits;

   logd("Running VQE with " + std::to_string(num_qubits) + " qubits", log_level);

   /*if(vqeOptions->loadIntermediate){
	   //double expected_energy, sv_energy, hit_rate;
	   bool success = loadProgress(vqeOptions->l_intermediateName, &initialParams, &expected_energy, &sv_energy, &hit_rate);
	   if(success){
		   if(verbose)
			   logi("Init params loaded from " + vqeOptions->l_intermediateName);
	   }else{
		   loge("Problem loading params from " + vqeOptions->l_intermediateName);
		   return;
	   }//
   }else{
	   // Init random parameters
	   //for(int i = 0; i < num_params_total; ++i){
		 //initialParams.emplace_back(dis(gen));
	//	   double p=dis(gen);logd(std::to_string(p));
		//   initialParams.emplace_back(p);
	   //}
	   //if(verbose)
	//	   logi("Random params generated");
   }
    */
   //std::function<void(int, double)> stats_function = stats_func;

   // Doesn't require to prepare the same circuit over and over again, but needs to clone statevect.
   //auto acc = xacc::getAccelerator("quest",
   //	   {{"nbQbits", observable->nBits()},
   //		{"startWithPlusState", true},
   //		{"repeated_measurement_strategy", true},
   //		{"hamiltonianProvided", vqeOptions->provideHamiltonian},
   //	    {"hamiltonianCoeffs", hamiltonian->coeffs},
   //		{"pauliCodes", hamiltonian->pauliOpts}
   //	   });

  // auto optimizer = xacc::getOptimizer("nlopt",
//		   {{"initial-parameters", initialParams}, {"nlopt-maxeval", max_iters}});

	ansatz = getAnsatz(options->ansatz_name, num_qubits, 1997);
	num_params = ansatz.num_params;

	logd("Ansatz generated. Executing vqe..");
	execute(buffer, options->accelerator, options->optimizer, options->zero_reference_states, hamiltonian, options->expectationToStandardOutput);
	logd("Vqe execution done");

	/*if(vqeOptions->saveIntermediate){
	   std::vector<double> params = buffer->opt_params;
	   double expected_energy = buffer->expected_energy;
	   double sv_energy = buffer->opt_val;
	   double hit_rate = buffer->hit_rate;

	   saveProgress(vqeOptions->s_intermediateName, params, expected_energy, sv_energy, hit_rate);
	}*/

	logi("Min QUBO found: " + std::to_string(buffer->opt_val), log_level);
	//std::vector<double> params = (*buffer)["opt-params"].as<std::vector<double>>();

	//vqeOptions->outfile.close();

}

void Vqe::execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* optimizer, std::vector<long long unsigned int> zero_reference_states, Hamiltonian* hamiltonian, bool logExpecStd){

	acc->options.zero_reference_states = zero_reference_states;
	acc->initialize(hamiltonian);

	std::vector<double> intermediateEnergies;
	acc->set_ansatz(&ansatz);

	if(acc->alpha_f(0,1,1,1)==0)
		logw("f: Constant Alpha: " + std::to_string(acc->options.samples_cut_ratio));
	else
		logw("f: Linear Init alpha: " + std::to_string(acc->options.samples_cut_ratio) + " Final alpha: " +	std::to_string(acc->options.final_alpha) + " Max iters: " + std::to_string(acc->options.max_alpha_iters));

	int iteration_i = 0;

	OptFunction f([&, this](const std::vector<double> &x, std::vector<double> &dx) {

		double ground_state_overlap;
		double expectation = acc->calc_expectation(buffer, x, iteration_i++, &ground_state_overlap);
		buffer->intermediateEnergies.push_back(expectation);
		buffer->intermediateGroundStateOverlaps.push_back(ground_state_overlap);
		if(logExpecStd)
			std::cout << instance_name << " "<< std::to_string(expectation) << std::endl;
		return (expectation);

	}, num_params);

	std::vector<double> initial_params;

	OptResult result;

	if(ansatz.circuit.qaoa_ansatz){

		std::vector<double> initialParams;
		//std::random_device rd;
		std::mt19937 gen(1997); //rd() instead of 1997
		std::uniform_real_distribution<> dis(-3.141592654, 3.141592654);

		for(int i = 0; i < ansatz.num_params/2; ++i){
			initial_params.push_back(dis(gen));
			initial_params.push_back(dis(gen));
		}

		std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
		std::vector<double> upperBounds(initial_params.size(), 3.141592654);

		result = optimizer->optimize(f, initial_params, 10e-6, max_iters, lowerBounds, upperBounds);

	}else{

		for(auto &gate : ansatz.circuit.gates)
			if(gate.param->name != "")
				initial_params.push_back(gate.param->value);

		std::vector<double> lowerBounds(initial_params.size(), -3.141592654);
		std::vector<double> upperBounds(initial_params.size(),  3.141592654);

		logw("are these really good ones?");

		result = optimizer->optimize(f, initial_params, 10e-6, max_iters, lowerBounds, upperBounds);
	}

	//double finalCost = result.first;

	std::string opt_config;

	acc->finalConfigEvaluator(buffer, result.second, nbSamples_calcVarAssignment);

	std::cout << "Final opt-val: " << buffer->opt_val << "\n";
	std::cout << "Final opt-config: " << buffer->opt_config << "\n";
	std::cout << "Final hit-rate: " << buffer->hit_rate << "\n";

	acc->finalize();

}
}
