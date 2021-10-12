#include <iostream>
#include <iterator>
#include <exception>
#include <chrono>

#include "popl.hpp"
#include "logger.h"
#include "executionStatistics.h"
#include "qaoa/qaoa.h"
#include "vqaConfig.h"
#include "indicators/progress_bar.hpp"
#include "latticeAlgorithms/iterativeLatticeReduction.h"
#include "lattice/hmlLattice.hpp"
#include <xacc.hpp>

#include "fastVQA.h"

#include "xacc_observable.hpp" //del
#include "PauliOperator.hpp" //del

#include "enumeration.h"
#include "run_paper_experiment.h"
#include "littleSombrero.h"

#include <mpi.h>

using namespace popl;

std::vector<std::pair<double, double>> ls_info;

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

		auto paper_exp		 = op.add<Switch>("e", "paperexp", "perform experiment as in the paper");
		auto rank_reduce 	 = op.add<Value<int>>("r", "", "rank truncation for paperexp", 0);
		auto circ_dir_prefix = op.add<Value<std::string>>("c", "circ-dir-prefix", "", "../experiment_files");

		auto save_ansatz	 = op.add<Switch>("s", "saveAnsatz", "save ansatz files");
		auto load_ansatz	 = op.add<Switch>("l", "loadAnsatz", "load ansatz files");

		auto littleSombrero = op.add<Switch>("", "littleSombrero", "perform little sombrero experiment");

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
		HmlLattice* hmlLat;
		Lattice* loadL;
		std::vector<Lattice> lattices_in;
		std::vector<AbstractLatticeInput*> lattices;

		if(qaoa -> is_set() /*|| enumeration->is_set()*/){

			if(paper_exp->is_set()){

				loge("Gaussian heuristics not implemented for low rank matrices. Returning 1.");

				SolutionDataset solutionDataset = run_paper_exp(25, 25, 180);
				std::pair<std::vector<MatrixInt>, std::vector<Solution>> dataset = solutionDataset.getMatricexAndDataset();
				std::vector<MatrixInt> matrices = std::get<0>(dataset);
				std::vector<Solution> solutions = std::get<1>(dataset);

				int i = rank_reduce->value() - solutionDataset.rank_min;
				if(i < 0 || i >= solutionDataset.num_ranks){
					logw("Invalid rank_reduce value. Will fetch the first matrix");
					i = 0;
				}

				for(auto &m: matrices){
					lattices.push_back(new Lattice(m, std::to_string(solutions[i].lattice_id)+"_"+std::to_string(solutions[i].rank)));
					std::cout << (std::to_string(solutions[i].lattice_id)+"_"+std::to_string(solutions[i].rank)) << " ";

					if(rank_reduce->value() == 0)
						i++;
					else
						i+=solutionDataset.num_ranks+1;
				}std::cout << "\n";
				num_lattices = lattices.size();
			}
			else if(!load_hml->is_set()){

				VqaConfig* vqaConfig;

				if(!config->is_set() && !littleSombrero->is_set()){

					if(!lattice_file->is_set()){
						loge("Neither config nor lattice file specified");
						return 1;
					}

					bool success;
					MatrixInt m = VqaConfig::loadLatticeFromFile(lattice_file->value(), &success);
					loadL = new Lattice(m, lattice_file->value());
					lattices.push_back(loadL);
					num_lattices = 1;

					if(success)
						logi(lattice_file->value() + " loaded");
					else
						loge("Problem loading " + lattice_file->value());

				}
				else if(littleSombrero->is_set()){

					LittleSombrero ls;
					for(std::pair<MatrixInt, std::string> &m: ls.loadLs())
						lattices.push_back(new Lattice(m.first, m.second));

					ls_info = ls.loadInfo();

				}else{

					//config file load
					try{
						vqaConfig = new VqaConfig(config->value());
					}
					catch(std::exception &e){
						loge(e.what());
						MPI_Finalize();
						return 1;
					}

					logi("Running QAOA. Configuration loaded from " + config->value());
					lattices_in = vqaConfig->getLattices();
					for(auto &l : lattices_in)
						lattices.push_back(&l);

					num_lattices = lattices.size();

					logi(std::to_string(num_lattices) + " lattice(s) in the config file");

				}

			}else{ //load from hml file

				loge("Not implemented error.");
				/*hml_lattice_mode = true;

				std::ifstream ifs(load_hml->value(), std::ios::binary | std::ios::in);
				if(!ifs.is_open()){
					loge("Cannot open " + load_hml->value());
				    MPI_Finalize();
					return 1;
				}

				int n;
				std::string n_str, hamiltonian;
				std::getline(ifs, n_str);
				std::getline(ifs, hamiltonian);
				std::istringstream iss(n_str);
				iss >> n;
				ifs.close();

				num_lattices = 1;
				hmlLat = new HmlLattice(n, hamiltonian);
				hmlLat->name = load_hml->value();*/
			}

			if(qaoa->is_set()){

				AcceleratorPartial accelerator = [overlap_penalty, overlap_trick, nbSamples, save_ansatz, load_ansatz, seed, circ_dir_prefix](std::shared_ptr<xacc::Observable> observable,
						bool hamiltonianExpectation,
						std::vector<double> hamCoeffs,
						std::vector<int>hamPauliCodes,
						std::string name){
					return xacc::getAccelerator("quest", {
							 std::make_pair("seed", seed),
							 std::make_pair("nbQbits", observable->nBits()),
							 // Doesn't require to prepare the same circuit over and over again, but needs to clone statevect.
							 std::make_pair("repeated_measurement_strategy", true),
							 std::make_pair("startWithPlusState", true),
							 std::make_pair("repeated_measurement_strategy", true),
							 std::make_pair("hamiltonianProvided", hamiltonianExpectation),
							 std::make_pair("hamiltonianCoeffs", hamCoeffs),
							 std::make_pair("pauliCodes", hamPauliCodes),
							 std::make_pair("overlapPenalty", overlap_penalty->value()),
							 std::make_pair("nbSamples", nbSamples->value()),
							 std::make_pair("name", name),
							 std::make_pair("overlapTrick", overlap_trick->is_set()),
   						     std::make_pair("saveAnsatz", save_ansatz->is_set()),
							 std::make_pair("loadAnsatz", load_ansatz->is_set()),
							 std::make_pair("circ_dir_prefix", circ_dir_prefix->value())
							 //std::make_pair("zero_config_statevect_index", overlap_trick->is_set() ? 1 : 0)
					});
				};

				OptimizerPartial optimizer = [](std::vector<double> initialParams, int max_iters) {
					return xacc::getOptimizer("nlopt", xacc::HeterogeneousMap {std::make_pair("initial-parameters", initialParams),
																			   std::make_pair("nlopt-maxeval", max_iters),
																			   std::make_pair("nlopt-ftol", 10e-9)});
				};

				QAOAOptions qaoaOptions;
				qaoaOptions.max_iters = niters->is_set() ? (niters->value() == 0 ? 5000 : niters->value()): 5000;
				qaoaOptions.detailed_log_freq = 0;//50;
				qaoaOptions.verbose = debug->is_set();
				qaoaOptions.debug = debug->is_set();
				qaoaOptions.verbose |= true;
				qaoaOptions.optimizer = optimizer;
				qaoaOptions.accelerator = accelerator;
				qaoaOptions.simplifiedSimulation = true;
				qaoaOptions.logEnergies = true;
				qaoaOptions.extendedParametrizedMode = false;//true;
				qaoaOptions.calcVarAssignment = true;
				qaoaOptions.provideHamiltonian = true;
				qaoaOptions.saveIntermediate = save_interm->is_set() ? (save_interm->value() == "" ? false : true) : false;
				qaoaOptions.s_intermediateName = qaoaOptions.saveIntermediate ? save_interm->value() : "";
				qaoaOptions.loadIntermediate = load_interm->is_set() ? (load_interm->value() == "" ? false : true) : false;
				qaoaOptions.l_intermediateName = qaoaOptions.loadIntermediate ? load_interm->value() : "";
				qaoaOptions.overlap_trick = overlap_trick->is_set();
				qaoaOptions.nbSamples_calcVarAssignment = nbSamples->value();
				qaoaOptions.save_ansatz = save_ansatz->is_set();
				qaoaOptions.load_ansatz = load_ansatz->is_set();

				MapOptions* mapOptions = new MapOptions();
				mapOptions->verbose = debug->is_set();
				mapOptions->num_qbits_per_x=qubits_per_x->value();
				if(overlap_trick->is_set() || overlap_penalty->value() == 0)
					mapOptions->pen_mode = MapOptions::no_hml_penalization;
				else
					mapOptions->pen_mode = MapOptions::penalty_all;

				ExecutionStatistics* execStats = new ExecutionStatistics();
				xacc::Initialize();
				xacc::setOption("quest-verbose", "true");
				xacc::setOption("quest-debug", "true");

				if(hml_lattice_mode){
					xacc::qbit* buffer;
					ProgressBar bar{bar_opts(0, 1, hmlLat->name)};
					qaoaOptions.set_default_stats_function(execStats, &bar, hmlLat);
					Qaoa::run_qaoa(&buffer, hmlLat->toHamiltonianString(), load_hml->value(), &bar, execStats, &qaoaOptions);
					logd("Hml mode");
				}else{
					int counter = 0;
					int prev_lattice_id=-1;
					for(auto &lattice_abs : lattices){

						int new_id = std::stoi(process(lattice_abs->name));
						if(prev_lattice_id == new_id)
							continue;

						prev_lattice_id = new_id;

						logi("Running " + lattice_abs->name);

						Lattice *lattice = static_cast<Lattice*>(lattice_abs);

						//THIS IS DONE A LITTLE BIT ABOVE
						//if(paper_exp->is_set() && rank_reduce->value())
						//	lattice->reduce_rank(rank_reduce->value());

						if(lll_preprocess->is_set()){
							lattice->lll_transformation = new MatrixInt(lattice->n_rows, lattice->n_cols);
							lattice->lll_preprocessed=true;
							lll_reduction(*(lattice->get_current_lattice()), *(lattice->lll_transformation), 0.99, 0.51, LLLMethod::LM_PROVED, FloatType::FT_DOUBLE);
						}

						lattice->toHamiltonianString(mapOptions); //remake, keep it here

						logw(lattice->toHamiltonianString(mapOptions));

						std::pair<std::vector<double>, std::vector<int>> hamiltonian2 = lattice->getHmlInQuestFormulation();

						ProgressBar bar{bar_opts(counter, num_lattices, lattice->name)};

						qaoaOptions.set_default_stats_function(execStats, &bar, lattice);
						if(qaoaOptions.overlap_trick)
							qaoaOptions.zero_reference_state = lattice->getZeroReferenceState();

						QOracle quantum_oracle = [&bar, execStats, &qaoaOptions, hamiltonian2]
												  (xacc::qbit** buffer, std::string hamiltonian, std::string name) {
							Qaoa::run_qaoa(buffer, hamiltonian, hamiltonian2, name, &bar, execStats, &qaoaOptions);
						};

						IterativeLatticeReduction ilr(lattice, mapOptions, quantum_oracle, 1);

						if(save_hml->is_set() && save_hml->value() != ""){
							std::ofstream ofs(save_hml->value(), std::ios::binary | std::ios::out);
							ofs << lattice->n_rows << " " << lattice->n_cols << "\n";
							ofs << lattice->toHamiltonianString(mapOptions);
							ofs.close();
						}

						ilr.run();

						counter++;
					}
				}

				//for(int i = 1; i < numRanks; ++i){
				//	MPI_Send(&control_val_exit, 1, MPI_INT, i, control_tag , MPI_COMM_WORLD);
				//}

				xacc::Finalize();
			} else{ //enum

				MapOptions* mapOptions = new MapOptions();
				mapOptions->num_qbits_per_x=qubits_per_x->value();
				mapOptions->verbose = false;
				mapOptions->penalty = 0;

				int qubits_fixed = 0;
				int nranks = numRanks;
				while (nranks >>= 1) ++qubits_fixed;

				if((1ULL<<qubits_fixed) != numRanks){
					loge("Num ranks must be power of 2. Exiting.");
					MPI_Finalize();
					return 1;
				}


				std::ofstream ofs;
				std::ofstream ofs_littleSombrero;

				if(rank == 0){
					ofs.open("enum_file.txt", std::ofstream::out | std::ofstream::trunc);
					ofs_littleSombrero.open("little_sombrero_out.txt", std::ofstream::out | std::ofstream::app);
				}

				int counter = 0;

				int ls_num_ranks = LittleSombrero::rank_high-LittleSombrero::rank_low+1;

				for(auto &lattice_abs : lattices){

					if(rank == 0)
						logi("Running " + lattice_abs->name);

					Lattice *lattice = static_cast<Lattice*>(lattice_abs);

					if(littleSombrero->is_set()){

						int ls_instance = counter / ls_num_ranks;
						int ls_rank = counter % ls_num_ranks + LittleSombrero::rank_low;

						logd(std::to_string(ls_instance) + " : " + std::to_string(ls_rank));

						if(counter % ls_num_ranks == 0)
							mapOptions->num_qbits_per_x = 1; //start with 1 qubit
						else if(ls_info[counter-1].second>=ls_info[counter].second-0.001){
							counter++;
							logd("Skip");
							ofs_littleSombrero << counter << " " << ls_instance << " " << ls_rank << " " << mapOptions->num_qbits_per_x << "\n";
							ofs_littleSombrero.flush();
							continue;
						}

						if(ls_rank * mapOptions->num_qbits_per_x > 30){

							ofs_littleSombrero << counter << " " << ls_instance << " " << ls_rank << " " << "skip" << "\n";
							ofs_littleSombrero.flush();
							counter++;
							continue;
						}

						while(1){

							lattice->toHamiltonianString(mapOptions); //remake, keep it here
							std::pair<std::vector<double>, std::vector<int>> hml = lattice->getHmlInQuestFormulation();

							int n = lattice->getNumQubits();
							int num_qubits = n;

							for(int j = n-1; j >= 0; --j){

								bool qbit_found = false;
								for(int i = 0; i < hml.first.size(); ++i){
									if(hml.second[i*n+j]==3){
										qbit_found = true;
										break;
									}
								}
								if(!qbit_found)
									num_qubits = j;
								else
									break;

							}

							bool* arr = (bool*)malloc((num_qubits-qubits_fixed)*sizeof(bool));

							if(rank == 0)
								logi("Enumeration over " + std::to_string(num_qubits) + " qubits");

							bool* fixed = (bool*)malloc(qubits_fixed*sizeof(bool));

							for(int i = 0; i < qubits_fixed; ++i){
								fixed[i] = rank & (1<<i);
							}

							EnumerationC e;
							e.n = n;
							e.rank = rank;
							e.fixed = fixed;
							e.n_qubits = num_qubits;
							e.qubits_fixed=qubits_fixed;
							e.lattice = lattice;
							e.coeffs = hml.first;
							e.pauliOpts = hml.second;
							e.ignoreZero = (mapOptions->penalty == 0);

							if(rank==0){
								ProgressBar bar{
								option::BarWidth{50},\
								option::Start{"["},\
								option::Fill{"="},\
								option::Lead{">"},\
								option::Remainder{" "},\
								option::End{"]"},\
								option::PrefixText{std::to_string(counter+1) + "/" + std::to_string(num_lattices)},\
								option::ForegroundColor{colors[counter % 7]},\
								option::ShowElapsedTime{true},\
								option::ShowRemainingTime{true},\
								option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},\
								option::MaxProgress{(1ULL<<(num_qubits-qubits_fixed))}};

								e.bar = &bar;
								//e.bar_tick = (1ULL<<(num_qubits-qubits_fixed)) / 100;

							}

							e.generateAllBinaryStrings(arr, 0);

							double shortestVect;

							MPI_Reduce(&e.shortestVectLen, &shortestVect, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

							double minFound = sqrt(shortestVect);

							if(rank == 0){

								//std::cerr<<"\nMinimum is: "  << shortestVect << "\n";
								std::cerr<<"\nMinimum_sqrt: "  << minFound << "\n";

								ofs_littleSombrero << counter << " " << ls_instance << " " << ls_rank << " " << mapOptions->num_qbits_per_x << "\n";
								ofs_littleSombrero.flush();

								//std::cerr<< "gh^2 = " << lattice->get_orig_gh().get_d() << "\n";
								int sum = 0;
								for(auto &x: lattice->get_current_lattice()->matrix[0])
										sum += x.get_si()*x.get_si();

								std::cerr<< "final |b1|: norm/gh = " << sqrt(shortestVect) / sqrt(   lattice->get_orig_gh().get_d()  ) << "\n";

								ofs << lattice->name << " " << shortestVect << " " << sqrt(sum) / sqrt(lattice->get_orig_gh().get_d()) << " " << sqrt(shortestVect) / sqrt( lattice->get_orig_gh().get_d()  ) << "\n";

							}

							logd("Comparing " + std::to_string(minFound) + " to " + std::to_string(ls_info[counter].second));

							if(minFound <= ls_info[counter].second + 0.001){
								logd("Success. Progressing to higher rank");
								break;
							}
throw;
							mapOptions->num_qbits_per_x++;
							logd("Increasing number of qubits to " + std::to_string(mapOptions->num_qbits_per_x));

						}

						counter++;

					}else{

						lattice->lll_transformation = new MatrixInt(lattice->n_rows, lattice->n_cols);
						lll_reduction(*(lattice->get_current_lattice()), *(lattice->lll_transformation), 0.99, 0.51, LLLMethod::LM_PROVED, FloatType::FT_DOUBLE);

						lattice->toHamiltonianString(mapOptions); //remake, keep it here
						//logw(lattice->toHamiltonianString(mapOptions));

						std::pair<std::vector<double>, std::vector<int>> hml = lattice->getHmlInQuestFormulation();

						int n = lattice->getNumQubits();
						int num_qubits = n;

						for(int j = n-1; j >= 0; --j){

							bool qbit_found = false;
							for(int i = 0; i < hml.first.size(); ++i){
								if(hml.second[i*n+j]==3){
									qbit_found = true;
									break;
								}
							}
							if(!qbit_found)
								num_qubits = j;
							else
								break;

						}

						bool* arr = (bool*)malloc((num_qubits-qubits_fixed)*sizeof(bool));

						if(rank == 0)
							logi("Enumeration over " + std::to_string(num_qubits) + " qubits");

						bool* fixed = (bool*)malloc(qubits_fixed*sizeof(bool));

						for(int i = 0; i < qubits_fixed; ++i){
							fixed[i] = rank & (1<<i);
						}

						EnumerationC e;
						e.n = n;
						e.rank = rank;
						e.fixed = fixed;
						e.n_qubits = num_qubits;
						e.qubits_fixed=qubits_fixed;
						e.lattice = lattice;
						e.coeffs = hml.first;
						e.pauliOpts = hml.second;
						e.ignoreZero = (mapOptions->penalty == 0);

						if(rank==0){
							ProgressBar bar{
							option::BarWidth{50},\
							option::Start{"["},\
							option::Fill{"="},\
							option::Lead{">"},\
							option::Remainder{" "},\
							option::End{"]"},\
							option::PrefixText{std::to_string(counter+1) + "/" + std::to_string(num_lattices)},\
							option::ForegroundColor{colors[counter % 7]},\
							option::ShowElapsedTime{true},\
							option::ShowRemainingTime{true},\
							option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},\
							option::MaxProgress{(1ULL<<(num_qubits-qubits_fixed))}};

							e.bar = &bar;
							//e.bar_tick = (1ULL<<(num_qubits-qubits_fixed)) / 100;

						}
						counter++;

						e.generateAllBinaryStrings(arr, 0);

						double shortestVect;

						MPI_Reduce(&e.shortestVectLen, &shortestVect, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

						if(rank == 0){

							std::cerr<<"\nMinimum is: "  << shortestVect << "\n";
							std::cerr<<"\nMinimum_sqrt: "  << sqrt(shortestVect) << "\n";

							std::cerr<< "gh^2 = " << lattice->get_orig_gh().get_d() << "\n";
							int sum = 0;
							for(auto &x: lattice->get_current_lattice()->matrix[0])
									sum += x.get_si()*x.get_si();
							if(!littleSombrero->is_set())
								std::cerr<< "LLL |b1|: norm/gh = " << sqrt(sum) / sqrt(lattice->get_orig_gh().get_d()) << "\n";
							std::cerr<< "final |b1|: norm/gh = " << sqrt(shortestVect) / sqrt(   lattice->get_orig_gh().get_d()  ) << "\n";

							ofs << lattice->name << " " << shortestVect << " " << sqrt(sum) / sqrt(lattice->get_orig_gh().get_d()) << " " << sqrt(shortestVect) / sqrt( lattice->get_orig_gh().get_d()  ) << "\n";

						}

					}

				}

				if(rank == 0){
					ofs.close();
					ofs_littleSombrero.close();
				}

			}
			int flag;
			MPI_Finalized(&flag);
		    if(!flag)
		    	MPI_Finalize();
			return 0;
		}

		else{
			std::cout << "Invalid argument. Use ./fastVQA -h to see help." << "\n";
		    MPI_Finalize();
			return 1;
		}
	}
	//--------------------------------RANK ZERO CODE FINISH------------------------
	//----------------------------------MSG RECEIVERS------------------------------
	else{
		logd("Hello I am rank " + std::to_string(rank) +"out of " + std::to_string(numRanks) + " and I am waiting for messages");
		xacc::Initialize();
		Qaoa::run_qaoa_slave_process();
		xacc::Finalize();
	}
    MPI_Finalize();
    return 0;
}
