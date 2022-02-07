#include "accelerator.h"

double Accelerator::evaluate_assignment(PauliHamil isingHam, std::string measurement){

	int numQubits = isingHam.numQubits;
	int z_index[2];

	double result = 0;//isingHam.termCoeffs[0];
	for(int term_i = 0; term_i < isingHam.numSumTerms; ++term_i){//auto &term: hamiltonian->getNonIdentitySubTerms()){

		//get the real part, imag expected to be 0
		double coeff = isingHam.termCoeffs[term_i];
		int num_zs = 0;

		for(int qbit_i = 0; qbit_i < numQubits; ++qbit_i){

			if(isingHam.pauliCodes[term_i*numQubits + qbit_i] == 3)
				z_index[num_zs++] = qbit_i;
			else if(isingHam.pauliCodes[term_i*numQubits + qbit_i] != 0){
				printf("Not a valid Ising hamiltonian");
				exit(1);
			}
		}

		if(num_zs == 0){
			result += coeff;
			continue;
		}
		else if(num_zs > 2){
			printf("Not a valid Ising hamiltonian");
			exit(1);
		}

		int multiplier = 1;

		int num_qubits=measurement.size();
		 //CRITICAL PART; STRING QUBIT ORDERING IS REVERSED!!
		multiplier *= (measurement[/*num_qubits-1-*/z_index[0]] == '0') ? 1 : -1;
		if(num_zs == 2)
			multiplier *= (measurement[/*num_qubits-1-*/z_index[1]] == '0') ? 1 : -1;

		result += multiplier * coeff;
	}

	return result;
}

void Accelerator::finalConfigEvaluator(ExperimentBuffer* buffer, std::vector<double> final_params, int nbSamples){

	bool second_eigenergy = true;

	/*typedef std::pair<int, double> meas_freq_eval;
	typedef std::pair<std::string, meas_freq_eval> measurement;

	std::vector<measurement> measurements;

	Qureg qureg_cache = createQureg(qureg.numQubitsInStateVec, env);

	//std::string classicalRefState_str(qureg.numQubitsRepresented, '0');
	//classicalRefState_str[qureg.numQubitsRepresented-1]='1';

	logd("Performing " + std::to_string(nbSamples) + " samples");

	for(int i = 0; i < nbSamples; ++i){

		cloneQureg(qureg_cache, qureg);
		std::string measurementStr = "";
		for(int measureQubit = qureg.numQubitsRepresented-1;
				measureQubit >= 0 ; --measureQubit){
			measurementStr += std::to_string(measure(qureg_cache, measureQubit));
		}

		bool found = false;
		int index = 0;
		for(auto &instance : measurements){
			if(instance.first == measurementStr){
				found = true;
				break;
			}
			index++;
		}

		//loge(classicalRefState_str);throw;

		//if(!overlapPenalization || measurementStr != classicalRefState_str){
		if(true)
			if(!found){

			  measurements.push_back(measurement(measurementStr, // bit string
									  meas_freq_eval(1, //frequency
											  evaluate_assignment(hamiltonian, measurementStr))));
			}else{
				measurements[index].second.first++;
			}
		}
	}
	if(false)
		  std::sort(measurements.begin(), measurements.end(),
			  [](const std::pair<std::string, std::pair<int, double>>& a, const std::pair<std::string, std::pair<int,int>>& b) {
				  //sort by global value
				  return a.second.second > b.second.second;
		  });

	else
		  std::sort(measurements.begin(), measurements.end(),
			  [](const std::pair<std::string, std::pair<int, double>>& a, const std::pair<std::string, std::pair<int,int>>& b) {
				  //sort by global value
				  return a.second.second < b.second.second;
		});

	for(auto &m : measurements){
		logw(m.first + " " + std::to_string(m.second.second));
	}

	measurement optimalMeasurement;
	if(second_eigenergy && measurements.size() > 1 && measurements[0].second.second == 0.){
		optimalMeasurement = measurements[1];
		logd("Returning second highest energy level");

		if(optimalMeasurement.second.second == 0){
			loge("Second highest energy is still 0!");
		}

	}else
		optimalMeasurement = measurements[0];

	int hits = optimalMeasurement.second.first;

	loge("Measurements size = " + std::to_string(measurements.size()));

	buffer->opt_config=optimalMeasurement.first;
	loge(optimalMeasurement.first);
	buffer->opt_val=optimalMeasurement.second.second;
	logw("Hits " + std::to_string(hits) + " / " + std::to_string(nbSamples));
	buffer->hit_rate= hits / double(nbSamples);

*/

	long long int i = 0;
	while(ref_hamil_energies[i++].first == 0);

	RefEnergy ground_state = ref_hamil_energies[i-1];
	i = ground_state.second;

	const int max_qubits=40;

	if(qureg.numQubitsInStateVec > max_qubits){
		loge("Not enough binary places to hold final opt_config. You are simulating more than" + std::to_string(max_qubits) + " qubits. Increase the number in the code.");
	}

	std::string opt_config = std::bitset<max_qubits>(i).to_string();
	opt_config=opt_config.substr(max_qubits-qureg.numQubitsInStateVec,qureg.numQubitsInStateVec);
	std::reverse(opt_config.begin(), opt_config.end());
	buffer->opt_config=opt_config;
	buffer->opt_val=ground_state.first;
	//logw("Hits " + std::to_string(hits) + " / " + std::to_string(nbSamples));
	buffer->hit_rate= qureg.stateVec.real[i]*qureg.stateVec.real[i]+qureg.stateVec.imag[i]*qureg.stateVec.imag[i];

	//std::sort(amps.begin(), amps.end();




}


void Accelerator::run(Circuit circuit, const std::vector<double> &x){

	int i = 0;
	double param_val=0;

	for(auto &gate : circuit.gates){

		if(gate.param->name != ""){
			param_val = x[i++];
		}

		apply_gate(gate, param_val);
	}

}

void Accelerator::set_ansatz(Ansatz* ansatz){

	this -> ansatz = *ansatz;

	int circuit_size = this->ansatz.circuit.gates.size();

	for(int i = 1; i < env.numRanks; ++i){
		MPI_Send(&circuit_size, 1, MPI_INT, i, new_circuit_tag , MPI_COMM_WORLD);
	}

	std::vector<double> gateCodes;
	for(auto &gate : this->ansatz.circuit.gates){
		gateCodes.push_back(gate.code);
		gateCodes.push_back(gate.qubit1);
		gateCodes.push_back(gate.qubit2);
	}

	MPI_Bcast(&gateCodes[0], circuit_size*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

double Accelerator::calc_expectation(ExperimentBuffer* buffer, const std::vector<double> &x, int iteration_i, double* ground_state_overlap_out){
	int x_size = x.size();
	if(env.numRanks > 1){
		loge("UNIMPLEMENTED");
		throw;
	}

	for(int i = 1; i < env.numRanks; ++i){
		MPI_Send(&x_size, 1, MPI_INT, i, calc_exp_val_tag , MPI_COMM_WORLD);
	}

	std::vector<double> x_copy(x);

	MPI_Bcast(&x_copy[0], x_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(ansatz.circuit.qaoa_ansatz){

		if(x.size() != ansatz.num_params){
			loge("Wrong number of parameters");
		}

		int p = ansatz.num_params/2;

		initPlusState(qureg);

		for(int i = 0; i < p; ++i){
			applyTrotterCircuit(qureg, hamiltonian,	x[2*i], 1, 1);
			multiRotatePauli(qureg, qubits_list, all_x_list, qureg.numQubitsInStateVec, x[2*i+1]);
		}

	}else{
		initZeroState(qureg);
		run(ansatz.circuit, x);
	}

	long long int ground_index = ref_hamil_energies[0].second;
	int ground_energy = ref_hamil_energies[0].first;
	if(ground_energy == 0)
		loge("Zero not excluded properly!");

	//code below calculate overlap with the shortest vector as a sum of overlaps over all feasible solutions
	int j = 0;
	*ground_state_overlap_out = 0;
	while(ref_hamil_energies[j].first <= ground_energy){
		long long int index = ref_hamil_energies[j].second;
		*ground_state_overlap_out += qureg.stateVec.real[index]*qureg.stateVec.real[index]+qureg.stateVec.imag[index]*qureg.stateVec.imag[index];
		j++;
	}
	std::cerr<<"This time we have " << j << " shortest vectors in the search space\n";

	/*qureg.stateVec.real[1<<zero_reference_state] = 0;
	qureg.stateVec.imag[1<<zero_reference_state] = 0;

	for(unsigned long i = 0; i < qureg.numAmpsPerChunk; ++i){
		norm+=qureg.stateVec.real[i]*qureg.stateVec.real[i]+qureg.stateVec.imag[0]*qureg.stateVec.imag[0];
	}*/

	/*for(auto &refEnergy : reference_energies_indexes){

		Complex amp = getAmp(qureg, refEnergy.second);
		energy += refEnergy.first * (amp.real*amp.real + amp.imag*amp.imag);

	}*/

	double energy=0;
	double alpha_sum=0;
	double amp;

	long long int i = 0;

	long long zero_state_index = options.zero_reference_state;

	double zero_state_amp = qureg.stateVec.real[zero_state_index]*qureg.stateVec.real[zero_state_index]+qureg.stateVec.imag[zero_state_index]*qureg.stateVec.imag[zero_state_index];
	double zero_state_amp_per_elem = zero_state_amp / (qureg.numAmpsPerChunk-1);

	double alpha = alpha_f(options.samples_cut_ratio, options.final_alpha, iteration_i, options.max_alpha_iters);
	//loge(std::to_string(alpha));

	while(alpha_sum < alpha){

		 if(i >= ref_hamil_energies.size()){
			 //logw("Probably something ain't alright");
			 break;
		 }

		long long int q_index = ref_hamil_energies[i].second;
		amp = zero_state_amp_per_elem + qureg.stateVec.real[q_index]*qureg.stateVec.real[q_index]+qureg.stateVec.imag[q_index]*qureg.stateVec.imag[q_index];
		alpha_sum += amp;
		if(i >= ref_hamil_energies.size()){
			loge("NOT ENOUGH REF ENERGIES STORED IN MEMORY.");
		}
		else{
			//loge(std::to_string(alpha));
			//logw(std::to_string(ref_hamil_energies[i].first) + " * " + std::to_string((amp / alpha)));
			energy += ref_hamil_energies[i].first * (amp / alpha);
			//loge(std::to_string(amp / alpha));
		}
		i++;
	}//std::cerr<<".";
	energy -= (alpha_sum - alpha) * ref_hamil_energies[i-1].first * (amp / alpha);
	//loge(std::to_string(-((alpha_sum - alpha))*(amp / alpha)));

	/*bool first=true;
	while(i--){
		std::cerr<<ref_hamil_energies[i].second<<" "<<ref_hamil_energies[i].first<<" p: "<<
		(first?(1-(alpha_sum - alpha)):1) * ((zero_state_amp_per_elem + qureg.stateVec.real[ref_hamil_energies[i].second]*qureg.stateVec.real[ref_hamil_energies[i].second]+qureg.stateVec.imag[ref_hamil_energies[i].second]*qureg.stateVec.imag[ref_hamil_energies[i].second])/alpha)<<"\n";
		first=false;
	}
	std::cerr<<"\n";*/

	/*if(energy == 0){
		std::cerr<<"." << i-1 << "\n";
		std::cerr<<ref_hamil_energies[i-1].first << " " << ref_hamil_energies[i-1].second << "\n";
		throw;
	}*/

	return energy/* / alpha*/ ; //normalize

	//QUEST CODE
		/*Complex localExpec = statevec_calcExpecDiagonalOpLocal(qureg, op);


	    if (qureg.numChunks == 1)
	        return localExpec;

	    qreal localReal = localExpec.real;
	    qreal localImag = localExpec.imag;
	    qreal globalReal, globalImag;
	    MPI_Allreduce(&localReal, &globalReal, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(&localImag, &globalImag, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);

	    Complex globalExpec;
	    globalExpec.real = globalReal;
	    globalExpec.imag = globalImag;
	    return globalExpec;*/
	  //QUEST CODE END

	//return calcExpecDiagonalOp(qureg, hamDiag).real;

}

void Accelerator::initialize(Hamiltonian* hamIn){

	int num_qubits = hamIn->nbQubits;

	qubits_list = new int[num_qubits]();
	all_x_list = new pauliOpType[num_qubits]();

	for(int i = 0; i < num_qubits; ++i){
		qubits_list[i]=i;
		all_x_list[i]=PAULI_X;
	}

	logd("Initializing " + std::to_string(num_qubits) + " qubits");

	for(int i = 1; i < env.numRanks; ++i){
		MPI_Send(&control_setting_new_seed, 1, MPI_INT, i, control_tag, MPI_COMM_WORLD);
	}

	unsigned long int keys[1];
	keys[0] = 1997;
	seedQuEST(&env, keys, 1);
	logd("Setting seed to " + std::to_string(keys[0]));

	for(int i = 1; i < env.numRanks; ++i){
		MPI_Send(&num_qubits, 1, MPI_INT, i, new_qureg_tag , MPI_COMM_WORLD);
	}

	qureg = createQureg(num_qubits, env);

	int coeffsSize = hamIn->coeffs.size();
	for(int i = 1; i < env.numRanks; ++i){
		MPI_Send(&coeffsSize, 1, MPI_INT, i, new_hamiltonian_tag , MPI_COMM_WORLD);
	}

	hamiltonian = createPauliHamil(num_qubits, coeffsSize);

	MPI_Bcast(&hamIn->coeffs[0], coeffsSize, MPI_QuEST_REAL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&hamIn->pauliOpts[0], coeffsSize*num_qubits, MPI_INT, 0, MPI_COMM_WORLD);

	hamiltonian.termCoeffs = &hamIn->coeffs[0]; //conversion to c array
	hamiltonian.pauliCodes = (enum pauliOpType*)(&hamIn->pauliOpts[0]);

	hamDiag = createDiagonalOp(num_qubits, env, 1);
	initDiagonalOpFromPauliHamil(hamDiag, hamiltonian);

	if(env.numRanks>1){
		loge("UNIMPLEMENTED");
		throw;
	}

	ref_hamil_energies.clear();

	std::vector<int> indexes(qureg.numAmpsPerChunk);
	std::iota(indexes.begin(), indexes.end(), 0); //zip with indices
	std::sort(indexes.begin(), indexes.end(), [&](int i, int j){return hamDiag.real[i] < hamDiag.real[j];}); //non-descending

	int counter = 0;
	for(auto &index : indexes){
		if(index == options.zero_reference_state){
			//logw("Zero excluded with counter " + std::to_string(counter) + " where E(0) = " + std::to_string(hamDiag.real[index]));
			if(hamDiag.real[index]!=0){
				loge("Tried to exclude something else than zero ground state! Did you run qaoa?");
				std::cerr<<"Tried to exclude index: " << index << " with energy " << hamDiag.real[index] << "\n";
			}
			else{
				continue;
			}
		}counter++;

		//logw(std::to_string(index)+"       " + std::to_string(hamDiag.real[index]));

		if(hamDiag.real[index] == 0){
			loge("Here we have not exluded zero with index " + std::to_string(index));
		}else
			ref_hamil_energies.push_back(RefEnergy(hamDiag.real[index], index));
		//if( double(counter++)/indexes.size() > options.samples_cut_ratio)
		//	break;
	}
}

void Accelerator::run_vqe_slave_process(){

	int value;
	MPI_Status status;

	bool quregInitialized = false;

	while(1){
		MPI_Recv(&value, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		if (status.MPI_TAG == control_tag) {
			if(value == control_val_exit){
					logd("Process " + std::to_string(env.rank) + " exiting work loop");
					break;
			}else if(value == control_setting_new_seed){
				unsigned long int keys[1];
				seedQuEST(&env, keys, 1);
				logd("Process " + std::to_string(env.rank) + " set seed to " + std::to_string(keys[0]));
			}else{
				loge("Invalid control value");
				throw;
			}

		}else if(status.MPI_TAG == new_circuit_tag){

			std::vector<double> ansatz_codes;

			ansatz_codes.resize(value*3);
  			MPI_Bcast(&ansatz_codes[0], value*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			logd("Process " + std::to_string(env.rank) + " received new ansatz codes");

			ansatz.circuit.gates.clear();

			for(int i = 0; i < ansatz_codes.size(); i+=3){
				ansatz.circuit.addGate(ansatz_codes[i],ansatz_codes[i+1],ansatz_codes[i+2]);
			}

		}else if(status.MPI_TAG == new_params_tag){

		}else if(status.MPI_TAG == new_qureg_tag){

			if(quregInitialized)
				destroyQureg(qureg, env);

			qureg = createQureg(value, env);
			quregInitialized = true;

			logd("Process " + std::to_string(env.rank) + " created new qureg");

		}else if(status.MPI_TAG == new_hamiltonian_tag){
			std::vector<double> coeffs;
			std::vector<int> pauliCodes;

			coeffs.resize(value);
			pauliCodes.resize(value*qureg.numQubitsInStateVec);

			MPI_Bcast(&coeffs[0], value, MPI_QuEST_REAL, 0, MPI_COMM_WORLD);
			MPI_Bcast(&pauliCodes[0], value*qureg.numQubitsInStateVec, MPI_INT, 0, MPI_COMM_WORLD);

			hamiltonian = createPauliHamil(qureg.numQubitsInStateVec, value);
			hamiltonian.termCoeffs = &coeffs[0]; //conversion to c array
			hamiltonian.pauliCodes = (enum pauliOpType*)(&pauliCodes[0]); //conversion to c array

			hamDiag = createDiagonalOp(qureg.numQubitsInStateVec, env, 1);
			initDiagonalOpFromPauliHamil(hamDiag, hamiltonian);

			logd("Process " + std::to_string(env.rank) + " initialized pauliHamil");

		}else if(status.MPI_TAG == calc_exp_val_tag){

			logd("Process " + std::to_string(env.rank) + " calcExpVal");

			std::vector<double> x(value);
			MPI_Bcast(&x[0], value, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			/*if(ansatz.circuit.qaoa_ansatz){

				int p=1;

				if(x.size() != p*2){
					loge("Wrong number of parameters");
				}



				initPlusState(qureg);
				loge("Plus state initialized");

				for(int i = 0; i < p; ++i){

				}

			}else*/
			loge("MPI not supported. throw");throw;

			initZeroState(qureg);
			run(ansatz.circuit, x);
			calcExpecDiagonalOp(qureg, hamDiag).real;

			//calcExpecDiagonalOp(qureg, hamDiag);
		}else if(status.MPI_TAG == measure_with_cache_tag){

		}else{
			loge("Unknown control tag");
			throw;
		}
	}

}

void Accelerator::finalize(){
	logd("Destroying qureg");
	destroyQureg(qureg, env);
}

Accelerator::Accelerator(AcceleratorOptions options){

	if(options.accelerator_type != "quest"){
		loge("No other accelerator than QuEST implemented");
		throw;
	}
	if(options.alpha_f == "constant"){
		this->alpha_f= alpha_constant_f;
	}else if(options.alpha_f == "linear"){
		this->alpha_f= alpha_linear_f;
	}else{
		loge("Unknown alpha function name");
		throw;
	}

	this->options = options;

}
