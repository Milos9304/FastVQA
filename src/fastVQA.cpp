#include "fastVQA.h"

#include <iostream>
#include <iterator>
#include <exception>
#include <chrono>

#include <fstream>
#include <iomanip>

#include "popl.hpp"
#include "executionStatistics.h"
#include "qaoa/qaoa.h"
#include "vqe/vqe.h"
#include "vqaConfig.h"
#include "latticeAlgorithms/iterativeLatticeReduction.h"
#include "lattice/hmlLattice.hpp"

#include "run_paper_experiment.h"

#include <mpi.h>

using namespace popl;

std::vector<std::pair<double, double>> ls_info;

std::string process(std::string const& s){
    std::string::size_type pos = s.find('_');
    if (pos != std::string::npos)
    {
        return s.substr(0, pos);
    }
    else
    {
        return s;
    }
}

int main(int ac, char** av){

	std::cout << fixed;
	std::cerr << fixed;

	int rank, numRanks;

	MPI_Init(&ac,&av);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int seed = 1997;

	//--------------------------------RANK ZERO CODE------------------------------------------
	if(rank == 0 /*|| strcmp(av[1],"qaoa")*/){

		OptionParser op("Allowed options");
		auto help_option     = op.add<Switch>("h", "help", "produce help message");
		auto qaoa 		     = op.add<Switch>("", "qaoa", "run qaoa algorithm");
		auto vqe 		     = op.add<Switch>("", "vqe", "run vqe algorithm");
		auto seed_option 	 = op.add<Value<int>>("", "seed", "seed for the experiments", seed);
		//auto enumeration     = op.add<Switch>("", "enum", "enumerate all qubo configurations");
		auto config 	     = op.add<Value<std::string>>("", "config", "config file location", "");
		auto lattice_file    = op.add<Value<std::string>>("", "lattice", "lattice file location", "");
		auto niters          = op.add<Value<int>>("i", "iters", "max num of iterations", 0);
		auto nbSamples 		 = op.add<Value<int>>("n", "nbSamples", "number of samples in var assigmnent", 1024);
		auto save_hml        = op.add<Value<std::string>>("", "savehml", "save hamiltonian to file", "");
		auto load_hml        = op.add<Value<std::string>>("", "loadhml", "save hamiltonian to file", "");
		auto debug           = op.add<Switch>("d", "debug", "print debug messages");
		auto qubits_per_x    = op.add<Value<int>>("q", "", "qubits per x", 1);
		auto overlap_trick   = op.add<Switch>("o", "", "perform overlap trick");
		auto overlap_penalty = op.add<Value<int>>("p", "", "overlap penalty", 0);
		auto lll_preprocess  = op.add<Switch>("", "lll", "perform LLL preprocessing on the lattice");
		auto second_eigval   = op.add<Switch>("", "second", "pick second lowest energy");

		auto paper_exp		 = op.add<Switch>("e", "paperexp", "perform experiment as in the paper");
		auto rank_reduce 	 = op.add<Value<int>>("r", "", "rank truncation for paperexp", 0);
		auto circ_dir_prefix = op.add<Value<std::string>>("c", "circ-dir-prefix", "", "../experiment_files");

		auto save_ansatz	 = op.add<Switch>("s", "saveAnsatz", "save ansatz files");
		auto load_ansatz	 = op.add<Switch>("l", "loadAnsatz", "load ansatz files");

		auto save_interm  = op.add<Value<std::string>>("", "si", "save intermediate results (for specific experiments only)", "");
		auto load_interm  = op.add<Value<std::string>>("", "li", "load intermediate results (for specific experiments only)", "");

		op.parse(ac, av);

		if (help_option->is_set()){
			std::cout << op << "\n";
		    MPI_Finalize();
			return 0;
		}

		seed = seed_option->value();
		logd("Using seed " + std::to_string(seed));

		int num_lattices;
		bool hml_lattice_mode=false;
		//HmlLattice* hmlLat;
		Lattice* loadL;
		std::vector<Lattice> lattices_in;
		std::vector</*AbstractLatticeInput**/Lattice*> lattices;

		if(paper_exp->is_set()){

			loge("Gaussian heuristics not implemented for low rank matrices. Returning 1.");

			SolutionDataset solutionDataset = run_paper_exp(11, 10, 180);
			std::pair<std::vector<MatrixInt>, std::vector<Solution>> dataset = solutionDataset.getMatricexAndDataset();
			std::vector<MatrixInt> matrices = std::get<0>(dataset);
			std::vector<Solution> solutions = std::get<1>(dataset);

			int i = rank_reduce->value() - solutionDataset.rank_min;
			if(i < 0 || i >= solutionDataset.num_ranks){
				logw("Invalid rank_reduce value. Will fetch the first matrix");
				i = 0;
			}

			for(auto &m: matrices){
				Lattice* new_lattice = new Lattice(m, std::to_string(solutions[i].lattice_id)+"_"+std::to_string(solutions[i].rank));
				new_lattice->reduce_rank(solutions[i].rank);
				lattices.push_back(new_lattice);
				std::cout << (std::to_string(solutions[i].lattice_id)+"_"+std::to_string(solutions[i].rank)) << " ";

				if(rank_reduce->value() == 0)
					i++;
				else
					i+=solutionDataset.num_ranks-1;

				//logw("Loading only one lattice");
				//break;

			}std::cout << "\n";
			num_lattices = lattices.size();
		}

		logd("Accelerator partial not yet returned");
		Accelerator accelerator("quest");
		accelerator.env=createQuESTEnv();
		logd("Accelerator partial returned");

		NLOptimizer optimizer;

		QAOAOptions *qaoaOptions;
		VQEOptions *vqeOptions;

		VQAOptions *vqaOptions = qaoa->is_set() ? new QAOAOptions() : (true? new VQEOptions() : new VQAOptions());
		vqaOptions->max_iters = niters->is_set() ? (niters->value() == 0 ? 5000 : niters->value()): 5000;
		vqaOptions->detailed_log_freq = 0;//50;
		vqaOptions->verbose = debug->is_set();
		vqaOptions->debug = debug->is_set();
		vqaOptions->verbose |= true;
		vqaOptions->optimizer = &optimizer;
		vqaOptions->accelerator = &accelerator;
		vqaOptions->logEnergies = true;
		vqaOptions->calcVarAssignment = true;
		vqaOptions->provideHamiltonian = true;
		vqaOptions->saveIntermediate = save_interm->is_set() ? (save_interm->value() == "" ? false : true) : false;
		vqaOptions->s_intermediateName = vqaOptions->saveIntermediate ? save_interm->value() : "";
		vqaOptions->loadIntermediate = load_interm->is_set() ? (load_interm->value() == "" ? false : true) : false;
		vqaOptions->l_intermediateName = vqaOptions->loadIntermediate ? load_interm->value() : "";
		vqaOptions->overlap_trick = overlap_trick->is_set();
		vqaOptions->nbSamples_calcVarAssignment = nbSamples->value();
		vqaOptions->save_ansatz = save_ansatz->is_set();
		vqaOptions->load_ansatz = load_ansatz->is_set();

		logd("VQAOptions set");

		if(qaoa->is_set()){
			qaoaOptions = static_cast<QAOAOptions*>(vqaOptions);
			qaoaOptions->simplifiedSimulation = true;
			qaoaOptions->extendedParametrizedMode = false;//true;
		}else if(vqe->is_set()){
			vqeOptions = static_cast<VQEOptions*>(vqaOptions);
		}

		MapOptions* mapOptions = new MapOptions();
		mapOptions->verbose = debug->is_set();
		mapOptions->num_qbits_per_x=qubits_per_x->value();
		mapOptions->penalty=overlap_penalty->value();
		if(overlap_trick->is_set() || overlap_penalty->value() == 0)
			mapOptions->pen_mode = MapOptions::no_hml_penalization;
		else
			mapOptions->pen_mode = MapOptions::penalty_all;

		ExecutionStatistics* execStats = new ExecutionStatistics();

		int counter = 0;
		int prev_lattice_id=-1;
		for(auto &lattice_abs : lattices){

			int new_id = std::stoi(process(lattice_abs->name));
			if(prev_lattice_id == new_id)
				continue;

			prev_lattice_id = new_id;

			logi("Running " + lattice_abs->name);

			Lattice *lattice = static_cast<Lattice*>(lattice_abs);

			if(qaoa->is_set()){

			}else if(vqe->is_set()){

				ProgressBar bar{bar_opts(counter, num_lattices, lattice->name, vqeOptions)};

				vqeOptions->set_default_stats_function(execStats, &bar, lattice);
				if(vqeOptions->overlap_trick)
					vqeOptions->zero_reference_state = lattice->getZeroReferenceState();

				Vqe vqe_instance;
				ExperimentBuffer buffer;

				Hamiltonian hamiltonian = lattice->getHamiltonian(mapOptions);

				logw("Before VQE run");
				vqe_instance.run_vqe(&buffer, &hamiltonian, lattice->name, &bar, execStats, vqeOptions);
				logw("VQE run");

				MatrixInt* current_lattice;
				int cols;

				if(lattice->lll_preprocessed){

				}else{
					current_lattice = lattice->get_orig_lattice();
					cols = current_lattice->c;
				}

				VectorInt x_vect = lattice->quboToXvector(buffer.opt_config);

				logd("NOW FOLLOWS X VECTOR");
					for(auto &x:x_vect)
						std::cerr<<x<<" ";
					std::cerr<<"\n";
				//loge("NOW FOLLOWS SHORTEST FOUND VECTOR");

				int double_check=0;
				for(int c = 0; c < cols; ++c){
					int res = 0;
					for(unsigned int r = 0; r < x_vect.size(); ++r){
						int a = current_lattice->matrix[r][c].get_si();
						int b = x_vect[r].get_si();
						res += b * a;
					}
					double_check += res*res;
					//std::cerr<<res<<" ";
				}//std::cerr<<"\n";

				double passed=true;
				passed = passed && (double_check == buffer.opt_val);

				logw("Double check: " + std::to_string(double_check) + " passed: " + (passed ? "TRUE" : "FALSE"));
				//if(!passed)
				//	throw;

				std::ofstream output_file("../experiment_files/statsfile_"+lattice->name+".txt");
				output_file << fixed << showpoint;
				output_file << setprecision(10);

				std::ostream_iterator<double> output_iterator(output_file, "\n");
				std::copy(buffer.intermediateEnergies.begin(), buffer.intermediateEnergies.end(), output_iterator);

				output_file << buffer.opt_val << " " << buffer.hit_rate << " ";
				for(auto &x:x_vect)
					output_file << x << " ";


			}else{
				throw;
			}

			counter++;
		}				//for(int i = 1; i < numRanks; ++i){
				//	MPI_Send(&control_val_exit, 1, MPI_INT, i, control_tag , MPI_COMM_WORLD);
				//}

				//xacc::Finalize();
				destroyQuESTEnv(accelerator.env);
				int flag;
				MPI_Finalized(&flag);
				if(!flag)
					MPI_Finalize();
				return 0;
			} else{
				//--------------------------------RANK ZERO CODE FINISH------------------------
				//----------------------------------MSG RECEIVERS------------------------------
				logd("Hello I am rank " + std::to_string(rank) +"out of " + std::to_string(numRanks) + " and I am waiting for messages");
				//xacc::Initialize();
				Vqe vqe_instance;

				Accelerator accelerator("quest");
				accelerator.env=createQuESTEnv();
				accelerator.run_vqe_slave_process();
				destroyQuESTEnv(accelerator.env);

				//xacc::Finalize();
			}
		MPI_Finalize();
		return 0;

}


