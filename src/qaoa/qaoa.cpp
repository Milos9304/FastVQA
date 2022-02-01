#include "qaoa.h"

#include <iostream>
//#include <xacc.hpp>

//#include "QuEST.h"
#include <random>

#include <fstream>
#include "mpi.h"

void Qaoa::run_qaoa(ExperimentBuffer* buffer,
		Hamiltonian* hamiltonian,
		std::string name,
		indicators::ProgressBar* bar,
		ExecutionStatistics* executionStats,
		QAOAOptions* qaoaOptions){

   int max_iters = qaoaOptions->max_iters;
   bool verbose = qaoaOptions->verbose;

   if(qaoaOptions->logEnergies){
	   //qaoaOptions->outfile.rdbuf()->pubsetbuf(0, 0); //disable buffer
	   qaoaOptions->outfile.open("../experiment_files/statsfile_"+name+".txt", std::fstream::out | std::ios_base::trunc); //| std::ios_base::trunc);//std::ios_base::app
   }

   if(bar)
	   qaoaOptions->progress_bar = bar;

   qaoaOptions->execStats = executionStats;

   loge("kokot");
   num_qubits = hamiltonian->nbQubits;
   ansatz = getAnsatz("qaoa", num_qubits, 1997);
   num_params = ansatz.num_params;


   logd("Executing qaoa");
   execute(buffer, qaoaOptions->accelerator, qaoaOptions->optimizer, qaoaOptions->zero_reference_state, hamiltonian);
   logd("Vqe execution done");

   if(qaoaOptions->saveIntermediate){
   std::vector<double> params = buffer->opt_params;
   double expected_energy = buffer->expected_energy;
   double sv_energy = buffer->opt_val;
   double hit_rate = buffer->hit_rate;

   saveProgress(qaoaOptions->s_intermediateName, params, expected_energy, sv_energy, hit_rate);
   }

	std::cout << "Min QUBO: " << buffer->opt_val << "\n";
	//std::vector<double> params = (*buffer)["opt-params"].as<std::vector<double>>();

	qaoaOptions->outfile.close();


   //xacc::setOption("quest-verbose", "true");
   //xacc::setOption("quest-debug", "true");

   //std::cout << hamiltonian <<"\n";throw;

   // The corresponding QUBO Hamiltonian is:



   /*


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
   //std::random_device rd;
   std::mt19937 gen(1997);
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
		 //initialParams.emplace_back(dis(gen));
		   double p=dis(gen);logd(std::to_string(p));
		   initialParams.emplace_back(p);
	   }
	   if(verbose)
		   logi("Random params generated");
   }

   auto qaoa = xacc::getService<xacc::Algorithm>("QAOA");

   //std::function<void(int, double)> stats_function = stats_func;

   bool initOk;
   std::vector<double> hamCoeffs;
   std::vector<int> hamPauliCodes;

   if(qaoaOptions->provideHamiltonian){
	   hamCoeffs = hamiltonian2.first;
	   hamPauliCodes = hamiltonian2.second;
   }



   auto optimizer = xacc::getOptimizer("nlopt",
		   {{"initial-parameters", initialParams}, {"nlopt-maxeval", max_iters}});

   if(qaoaOptions->isSetLogStats()){
	   initOk = qaoa->initialize({
		   	   std::make_pair("name", name),
			   std::make_pair("accelerator", qaoaOptions->accelerator(observable,
					   qaoaOptions->provideHamiltonian,
					   hamCoeffs,
					   hamPauliCodes,
					   name)),
			   std::make_pair("optimizer", qaoaOptions->optimizer(initialParams, max_iters)),
			   std::make_pair("observable", observable),
			   std::make_pair("detailed_log_freq", qaoaOptions->detailed_log_freq),
			   // number of time steps (p) param
			   std::make_pair("steps", qaoaOptions->p),
			   std::make_pair("load_ansatz",  qaoaOptions->load_ansatz),
			   std::make_pair("calc-var-assignment", qaoaOptions->calcVarAssignment),
			   std::make_pair("simplified-simulation", qaoaOptions->simplifiedSimulation),
			   std::make_pair("stats_func", qaoaOptions->get_stats_function()),
			   //Number of samples to estimate optimal variable assignment
			   std::make_pair("parameter-scheme", qaoaOptions->getParameterScheme()),
			   std::make_pair("nbSamples", qaoaOptions->nbSamples_calcVarAssignment),
			   std::make_pair("questHamExpectation", qaoaOptions->provideHamiltonian),
			   std::make_pair("debugMsgs", qaoaOptions->debug),
			   std::make_pair("overlapTrick", qaoaOptions->overlap_trick),
			   std::make_pair("zeroRefState", qaoaOptions->zero_reference_state),
		});
   }
   else{
	   initOk = qaoa->initialize({
		   	   std::make_pair("name", name),
		       std::make_pair("accelerator", qaoaOptions->accelerator(observable,
		    		   qaoaOptions->provideHamiltonian,
					   hamCoeffs,
					   hamPauliCodes,
					   name)),
		   	   std::make_pair("optimizer", qaoaOptions->optimizer(initialParams, max_iters)),
			   std::make_pair("detailed_log_freq", qaoaOptions->detailed_log_freq),
	   		   // number of time steps (p) param
	   		   std::make_pair("steps", qaoaOptions->p),
			   std::make_pair("load_ansatz",  qaoaOptions->load_ansatz),
	   		   std::make_pair("calc-var-assignment", qaoaOptions->calcVarAssignment || qaoaOptions->load_ansatz),
	   		   std::make_pair("simplified-simulation", qaoaOptions->simplifiedSimulation),
	   		   //Number of samples to estimate optimal variable assignment
			   std::make_pair("parameter-scheme", qaoaOptions->getParameterScheme()),
	   		   std::make_pair("nbSamples", qaoaOptions->nbSamples_calcVarAssignment),
			   std::make_pair("questHamExpectation", qaoaOptions->provideHamiltonian),
			   std::make_pair("debugMsgs", qaoaOptions->debug),
			   std::make_pair("overlapTrick", qaoaOptions->overlap_trick),
			   std::make_pair("zeroRefState", qaoaOptions->zero_reference_state),
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
*/
}

void Qaoa::execute(ExperimentBuffer* buffer, Accelerator* acc, Optimizer* optimizer, int zero_reference_state, Hamiltonian* hamiltonian){

	acc->options.zero_reference_state = zero_reference_state;
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
		//logw(std::to_string(expectation));
		return (expectation);

	}, num_params);

	std::vector<double> initial_params;

	OptResult result;

		std::vector<double> initialParams;
		//std::random_device rd;
		std::mt19937 gen(1997); //rd() instead of 1997
		std::uniform_real_distribution<> dis(0, 2*3.141592654);

		for(int i = 0; i < ansatz.num_params/2; ++i){
			initial_params.push_back(dis(gen));
			initial_params.push_back(dis(gen));
		}

		std::vector<double> lowerBounds(initial_params.size(), 0);
		std::vector<double> upperBounds(initial_params.size(), 2*3.141592654);

		result = optimizer->optimize(f, initial_params, 10e-6, max_iters, lowerBounds, upperBounds);

	double finalCost = result.first;

	std::string opt_config;

	/*VQEOptimalConfigEvaluator::evaluate(
	  			  observable,
	  			  kernel->operator()(result.second),
	  			  accelerator,
	  			  nbSamples,
	  			  nbQubits,
	  			  m_maximize,
	  			  buffer,
	  			  &optimal_energy,
	  			  &opt_config,
	  			  &hit_rate);
	 */
	acc->finalConfigEvaluator(buffer, result.second, /*1024*/5000);

	std::cout << "Final opt-val: " << buffer->opt_val << "\n";
	std::cout << "Final opt-config: " << buffer->opt_config << "\n";
	std::cout << "Final hit-rate: " << buffer->hit_rate << "\n";

	acc->finalize();

}


/*void Qaoa::run_qaoa_slave_process(){
	//auto acc = xacc::getAccelerator("quest");
	loge("NOT IMPLEMENTED");
	throw;
}*/
