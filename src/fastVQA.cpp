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

std::string to_string_with_precision(double a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

int main(int ac, char** av){

	std::cout << fixed;
	std::cerr << fixed;

	int rank, numRanks;

	MPI_Init(&ac,&av);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int seed = 1997;

	AcceleratorOptions acceleratorOptions;
	acceleratorOptions.accelerator_type = "quest";


	//--------------------------------RANK ZERO CODE------------------------------------------
	if(rank == 0 /*|| strcmp(av[1],"qaoa")*/){

		OptionParser op("Allowed options");
		auto help_option     = op.add<Switch>("h", "help", "produce help message");
		auto qaoa 		     = op.add<Switch>("", "qaoa", "run qaoa algorithm");
		auto vqe 		     = op.add<Switch>("", "vqe", "run vqe algorithm");
		auto super_instance  = op.add<Switch>("", "superInstance", "");
		auto ansatz_name     = op.add<Value<std::string>>("a", "ansatz","Ry_CNOT_all2all_Ry/qaoa/EfficientSU2", "Ry_CNOT_all2all_Ry");
		auto seed_option 	 = op.add<Value<int>>("", "seed", "seed for the experiments", seed);
		//auto enumeration     = op.add<Switch>("", "enum", "enumerate all qubo configurations");
		auto config 	     = op.add<Value<std::string>>("", "config", "config file location", "");
		auto lattice_file    = op.add<Value<std::string>>("", "lattice", "lattice file location", "");
		auto niters          = op.add<Value<int>>("i", "iters", "max num of iterations", 1000);
		auto nbSamples 		 = op.add<Value<int>>("n", "nbSamples", "number of samples in var assigmnent", 1024);
		auto save_hml        = op.add<Value<std::string>>("", "savehml", "save hamiltonian to file", "");
		auto load_hml        = op.add<Value<std::string>>("", "loadhml", "save hamiltonian to file", "");
		auto debug           = op.add<Switch>("d", "debug", "print debug messages");
		auto qubits_per_x    = op.add<Value<int>>("q", "", "qubits per x", 1);
		auto overlap_trick   = op.add<Switch>("o", "", "perform overlap trick");
		auto overlap_penalty = op.add<Value<int>>("p", "", "overlap penalty", 0);
		auto lll_preprocess  = op.add<Switch>("", "lll", "perform LLL preprocessing on the lattice");
		auto second_eigval   = op.add<Switch>("", "second", "pick second lowest energy");

		auto initial_alpha   = op.add<Value<double>>("x", "alpha", "initial alpha value", 1);
		auto linear_alpha    = op.add<Switch>("l", "linear", "linear alpha");
		auto final_alpha     = op.add<Value<double>>("f", "final_alpha", "final alpha value", 0.5);
		auto max_alpha_iters = op.add<Value<int>>("m", "max_alpha_iters", "max alpha iters", 1000);

		auto paper_exp		 = op.add<Switch>("e", "paperexp", "perform experiment as in the paper");
		auto rank_reduce 	 = op.add<Value<int>>("r", "", "rank truncation for paperexp", 0);
		auto circ_dir_prefix = op.add<Value<std::string>>("c", "circ-dir-prefix", "", "../experiment_files");

		auto save_ansatz	 = op.add<Switch>("", "saveAnsatz", "save ansatz files");
		auto load_ansatz	 = op.add<Switch>("", "loadAnsatz", "load ansatz files");

		auto save_interm  = op.add<Value<std::string>>("", "si", "save intermediate results (for specific experiments only)", "");
		auto load_interm  = op.add<Value<std::string>>("", "li", "load intermediate results (for specific experiments only)", "");

        auto statsfile_prefix = op.add<Value<std::string>>("", "prefix", "statsfile prefix", "rank");

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

		std::string qubit_assignment_str;
		std::vector<int> qubit_assignment;

		if(paper_exp->is_set() || super_instance->is_set()){

			acceleratorOptions.samples_cut_ratio=initial_alpha->value();
			acceleratorOptions.final_alpha=final_alpha->value();
			acceleratorOptions.max_alpha_iters=max_alpha_iters->value();
			acceleratorOptions.alpha_f = linear_alpha->is_set() ? "linear" : "constant";

			logw("Setting accelerator ratio to " + std::to_string(acceleratorOptions.samples_cut_ratio));
			loge("This won't work for distributed. Change code!");

			loge("Gaussian heuristics not implemented for low rank matrices. Returning 1.");

			std::shared_ptr<SolutionDataset> solutionDataset;
			if(paper_exp-> is_set())
				solutionDataset = std::make_shared<SolutionDataset>(run_paper_exp(16, 10, 180, "out_higher_dims_small_dims"));
			else{
				//SolutionDataset dataset_temp =
				qubit_assignment_str = "     0    2     2    0    0    2    0    2    0    2    2    0     2    2    0    0    2    0    0    2    2    2    0    2    0    0    2    0    0    2";
				//                        (0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

				solutionDataset = std::make_shared<SolutionDataset>(run_paper_exp(30, 1, 180, "super_instance"));
				int num_qubits_assigned = 0;

				for(int i = 0; i < qubit_assignment_str.size(); ++i){
					switch(qubit_assignment_str[i]){
						case '0':
							qubit_assignment.push_back(0);
							break;
						case '1':
							num_qubits_assigned += 1;
							qubit_assignment.push_back(1);
							break;
						case '2':
							num_qubits_assigned += 2;
							qubit_assignment.push_back(2);
							break;
						case ' ':
							break;
						default:
							loge("Invalid qubit assignment string");
							throw;
					}
				}

				assert(num_qubits_assigned == 30);

			}

			std::pair<std::vector<MatrixInt>, std::vector<Solution>> dataset = solutionDataset->getMatricexAndDataset();
			std::vector<MatrixInt> matrices = std::get<0>(dataset);
			std::vector<Solution> solutions = std::get<1>(dataset);

			int i = rank_reduce->value() - solutionDataset->rank_min;
			if(i < 0 || i >= solutionDataset->num_ranks){
				logw("Invalid rank_reduce value. Will fetch the first matrix");
				i = 0;
			}

			int counter = 0;
			int max_instances = 128;

			for(auto &m: matrices){

				if(counter++>=max_instances){
					logw("Limiting to lower number of instances!");
					break;
				}
				Lattice* new_lattice = new Lattice(m, std::to_string(solutions[i].lattice_id)+"_"+std::to_string(solutions[i].rank));
				new_lattice->reduce_rank(solutions[i].rank);
				lattices.push_back(new_lattice);
				std::cout << (std::to_string(solutions[i].lattice_id)+"_"+std::to_string(solutions[i].rank)) << " ";

				if(rank_reduce->value() == 0)
					i++;
				else
					i+=solutionDataset->num_ranks;//-1;

				//logw("Loading only one lattice");
				//break;

			}std::cout << "\n";
			num_lattices = lattices.size();
		}

		logd("Accelerator partial not yet returned");

		Accelerator accelerator(acceleratorOptions);
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
		vqaOptions->ansatz_name = ansatz_name->value();

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
		if(qubit_assignment.size() > 0)
			mapOptions->qubit_assignment_optional = qubit_assignment;

		if(!overlap_trick->is_set()){
			mapOptions->pen_mode = MapOptions::no_hml_penalization;
			logw("no hml penalization");
		}
		else{
			mapOptions->pen_mode = MapOptions::penalty_all;
			logw("penalty all");
		}
		ExecutionStatistics* execStats = new ExecutionStatistics();

		int counter = 0;
		int prev_lattice_id=-1;

		//int skip[25] = {6, 17, 27, 32, 35, 40, 43, 44, 51, 59, 61, 65, 66, 70, 77, 82, 85, 87, 100, 112, 117, 122, 124, 125, 127};

		for(auto &lattice_abs : lattices){

			int new_id = std::stoi(process(lattice_abs->name));
			if(prev_lattice_id == new_id)
				continue;

			prev_lattice_id = new_id;

			logi("Running " + lattice_abs->name);

			/*if(counter < 5){
				counter++;
				continue;
			}else if(counter > 5)
				break;*/

			bool excluded = false;
			/*for(int i = 0; i < 25; ++i){
				if(counter == skip[i])
					excluded = true;
			}*/

			if(excluded){
				logw("Excluded");
				counter++;
				continue;
			}

			Lattice *lattice = static_cast<Lattice*>(lattice_abs);

			std::string statsfile="../experiment_files/rank_"+std::to_string(rank_reduce->value())+"/"+statsfile_prefix->value()+"_statsfile_"+lattice->name+".txt";
			loge(statsfile);
			ifstream ifile;
			ifile.open(statsfile);
			if(ifile){
				logw(statsfile + " exists");
				counter++;
				continue;
			}
			ifile.close();

			if(qaoa->is_set()){
				logi("QAOAs");
				Qaoa qaoa_instance;
				ExperimentBuffer buffer;
				ProgressBar bar{bar_opts(counter, num_lattices, lattice->name, qaoaOptions)};

				int first_vect_sq_len = 0;
				for(int c = 0; c < lattice->get_orig_lattice()->c; ++c){
					int a = lattice->get_orig_lattice()->matrix[0][c].get_si();
					first_vect_sq_len += a*a;
				}

				mapOptions->penalty = 2*first_vect_sq_len;
				logw("Penalty set to " + std::to_string(mapOptions->penalty));

				Hamiltonian hamiltonian = lattice->getHamiltonian(mapOptions);

				qaoa_instance.run_qaoa(&buffer, &hamiltonian, lattice->name, &bar, execStats, qaoaOptions);


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

								std::ofstream output_file("../experiment_files/rank_"+std::to_string(rank_reduce->value())+"/"+statsfile_prefix->value()+"_statsfile_"+lattice->name+".txt");
								output_file << fixed << showpoint;
								output_file << setprecision(10);

								//std::ostream_iterator<double> output_iterator(output_file, "\n");
								//std::copy(buffer.intermediateEnergies.begin(), buffer.intermediateEnergies.end(), output_iterator);
								for(int i = 0; i < buffer.intermediateEnergies.size(); ++i){
									output_file << buffer.intermediateEnergies[i] << " " << buffer.intermediateGroundStateOverlaps[i] << std::endl;
								}

								output_file << buffer.opt_val << " " << buffer.hit_rate << " ";
								for(auto &x:x_vect)
									output_file << x << " ";

								output_file.close();




			}else if(vqe->is_set()){
				logi("VQE");
				ProgressBar bar{bar_opts(counter, num_lattices, lattice->name, vqeOptions)};

				vqeOptions->set_default_stats_function(execStats, &bar, lattice);
				//if(vqeOptions->overlap_trick)

				Vqe vqe_instance;
				ExperimentBuffer buffer;

				Hamiltonian hamiltonian = lattice->getHamiltonian(mapOptions);
				vqeOptions->zero_reference_state = lattice->getZeroReferenceState();

				logw("Before VQE run");
				vqe_instance.run_vqe(&buffer, &hamiltonian, lattice->name, &bar, execStats, vqeOptions, seed);
				logw("VQE run");

				MatrixInt* current_lattice;
				int cols;

				if(lattice->lll_preprocessed){

				}else{
					current_lattice = lattice->get_orig_lattice();
					cols = current_lattice->c;
				}
				//qubit_assignment_str = "     0    2     2    0    0    2    0    2    0    2    2    0     2    2    0    0    2    0    0    2    2    2    0    2    0    0    2    0    0    2";
				//buffer.opt_config = "000000000000000000000000000000";


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

				std::ofstream output_file("../experiment_files/rank_"+std::to_string(rank_reduce->value())+"/"+statsfile_prefix->value()+"_statsfile_"+lattice->name+".txt");
				output_file << fixed << showpoint;
				output_file << setprecision(10);

				//std::ostream_iterator<double> output_iterator(output_file, "\n");
				//std::copy(buffer.intermediateEnergies.begin(), buffer.intermediateEnergies.end(), output_iterator);
				for(int i = 0; i < buffer.intermediateEnergies.size(); ++i){
					output_file << buffer.intermediateEnergies[i] << " " << buffer.intermediateGroundStateOverlaps[i] << std::endl;
				}
				std::cerr<<setprecision(15)<<"Final: " << buffer.opt_val << " " << buffer.hit_rate << " ";
				output_file << buffer.opt_val << " " << buffer.hit_rate << " ";
				for(auto &x:x_vect)
					output_file << x << " ";

				output_file.close();

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

				Accelerator accelerator(acceleratorOptions);

				accelerator.env=createQuESTEnv();
				accelerator.run_vqe_slave_process();
				destroyQuESTEnv(accelerator.env);

				//xacc::Finalize();
			}
		MPI_Finalize();
		return 0;

}


