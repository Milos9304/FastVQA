#include "accelerator.h"
#include "logger.h"
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cfloat>

//#include <chrono>

# if QuEST_PREC==1
	#define QREAL_MAX FLOAT_MAX
# elif QuEST_PREC==2
    #define QREAL_MAX DOUBLE_MAX
# elif QuEST_PREC==4
   #define QREAL_MAX LONG_DOUBLE_MAX
# endif

namespace FastVQA{

const int max_qubits=40;

AlphaFunction Accelerator::alpha_constant_f = [](double init_val, double final_val, int iter_i, int max_iters){
	return init_val;
};

AlphaFunction Accelerator::alpha_linear_f = [](double init_val, double final_val, int iter_i, int max_iters){

	int steps=35;

	if(iter_i < max_iters*steps){
		int j = iter_i/steps;
		return init_val+((double)j/max_iters)*(final_val-init_val);
	}
	/*if(iter_i < max_iters)
		return init_val+((double)iter_i/max_iters)*(final_val-init_val);*/

	return final_val;
};

RefEnergies Accelerator::getEigenspace(){
	if(!ref_hamil_energies_set)
		throw_runtime_error("You need to run accelerator->initialize first!");
	return ref_hamil_energies;
}

RefEnergies Accelerator::getSolutions(){
	if(!ref_hamil_energies_set)
		throw_runtime_error("You need to run accelerator->initialize first!");

	RefEnergies res;

	//int min_sol=ref_hamil_energies[0].first;
	//int i = 0;
	//while(ref_hamil_energies[i].first == min_sol){
	//	res.push_back(ref_hamil_energies[i]);
	//	i++;
	//}

	for(auto &e: ref_hamil_energies){
		if(e.isConsideredSolution)
			res.push_back(e);
	}

	return res;
}

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

		//int num_qubits = measurement.size();
		 //CRITICAL PART; STRING QUBIT ORDERING IS REVERSED!!
		multiplier *= (measurement[/*num_qubits-1-*/z_index[0]] == '0') ? 1 : -1;
		if(num_zs == 2)
			multiplier *= (measurement[/*num_qubits-1-*/z_index[1]] == '0') ? 1 : -1;

		result += multiplier * coeff;
	}

	return result;
}

void Accelerator::finalConfigEvaluator(ExperimentBuffer* buffer, std::vector<double> final_params, int nbSamples){

	//if(!this->options.exclude_zero_state && ref_hamil_energies[0].value == 0)
	//	loge("Zero not excluded properly!");

	RefEnergies ground_states = this->getSolutions();

	//if(this->options.exclude_zero_state){
		/*
		qreal min_gs = DBL_MAX;
		for(auto &hamil_energy : ref_hamil_energies){
			if(hamil_energy.first != 0 && hamil_energy.first < min_gs)
				min_gs = hamil_energy.first;
		}
		for(auto &hamil_energy : ref_hamil_energies){
			if(hamil_energy.first == min_gs)
				ground_states.push_back(hamil_energy);
		}*/
	//}
	/*else{
		else{loge("here");
			std::sort(ref_hamil_energies.begin(), ref_hamil_energies.end(),
					[](const RefEnergy& a, const RefEnergy& b) {
						return a.first < b.first;
			});

			if(ref_hamil_energies[0].first == 0)
				loge("Zero not excluded");

			int j = 0;
			while(ref_hamil_energies[j++].first == ref_hamil_energies[0].first){
				ground_states.push_back(ref_hamil_energies[j-1]);
			}
		}
	}*/

	if(qureg.numQubitsInStateVec > max_qubits){
		loge("Not enough binary places to hold final opt_config. You are simulating more than" + std::to_string(max_qubits) + " qubits. Increase the number in the code.");
	}

	long long int i = 0;
	for(auto &ground_state: ground_states){
		i = ground_state.index;
		std::string opt_config = std::bitset<max_qubits>(i).to_string();
		opt_config=opt_config.substr(max_qubits-qureg.numQubitsInStateVec, qureg.numQubitsInStateVec);
		std::reverse(opt_config.begin(), opt_config.end());

		ExperimentBuffer::ExperimentBufferSolution solution(opt_config, ground_state.value, qureg.stateVec.real[i]*qureg.stateVec.real[i]+qureg.stateVec.imag[i]*qureg.stateVec.imag[i]);
		buffer->final_solutions.push_back(solution);
	}

	buffer->opt_val=ground_states[0].value;
	//buffer->opt_config=opt_config;
	//logw("Hits " + std::to_string(hits) + " / " + std::to_string(nbSamples));
	//buffer->hit_rate= qureg.stateVec.real[i]*qureg.stateVec.real[i]+qureg.stateVec.imag[i]*qureg.stateVec.imag[i];

	//std::sort(amps.begin(), amps.end();
}

void Accelerator::run_with_new_params(Circuit circuit, const std::vector<double> &x){

	int i = 0;
	double param_val=0;

	for(auto &gate : circuit.gates){

		if(gate.param->name != ""){
			param_val = x[i++];
		}

		apply_gate(gate, param_val);
	}

}

double Accelerator::_energy_evaluation(double* ground_state_overlap_out, int iteration_i){

	//if(!this->options.exclude_zero_state && ref_hamil_energies[0].value == 0)
	//		logd("Zero was not excluded properly!", this->log_level);

	//OLD VQE WAY

	/*long long int ground_index = ref_hamil_energies[0].index;
	*ground_state_overlap_out = qureg.stateVec.real[ground_index]*qureg.stateVec.real[ground_index]+qureg.stateVec.imag[ground_index]*qureg.stateVec.imag[ground_index];

	double energy=0;
	double alpha_sum=0;
	double amp=0;

	long long unsigned int i = 0;

	std::vector<long long unsigned int> zero_state_indices = options.zero_reference_states;
	long long unsigned int zero_state_index;
	if(zero_state_indices.size() > 1){
		loge("TODO: Multiple zero_state_indices not implemented yet");
	}zero_state_index = zero_state_indices.size() > 0 ? zero_state_indices[0] : -1;

	double zero_state_amp;
	if(ansatz.circuit.qaoa_ansatz || zero_state_indices.size() == 0)
		zero_state_amp = 0;
	else
		zero_state_amp = qureg.stateVec.real[zero_state_index]*qureg.stateVec.real[zero_state_index]+qureg.stateVec.imag[zero_state_index]*qureg.stateVec.imag[zero_state_index];

	double zero_state_amp_per_elem = zero_state_amp / (qureg.numAmpsPerChunk-1);

	double alpha = 0.2;//alpha_f(options.samples_cut_ratio, options.final_alpha, iteration_i, options.max_alpha_iters);

	while(alpha_sum < alpha){

		 if(i >= ref_hamil_energies.size()){
			 break;
		 }

		long long int q_index = ref_hamil_energies[i].index;
		amp = zero_state_amp_per_elem + qureg.stateVec.real[q_index]*qureg.stateVec.real[q_index]+qureg.stateVec.imag[q_index]*qureg.stateVec.imag[q_index];
		alpha_sum += amp;
		if(i >= ref_hamil_energies.size()){
			loge("NOT ENOUGH REF ENERGIES STORED IN MEMORY.");
		}
		else{
			energy += ref_hamil_energies[i].value * (amp / alpha);
		}
		i++;
	}
	energy -= (alpha_sum - alpha) * ref_hamil_energies[i-1].value * (amp / alpha);*/

	double energy = 0;
	for(auto &e: ref_hamil_energies){
		long long int q_index = e.index;
		energy += e.value * (qureg.stateVec.real[q_index]*qureg.stateVec.real[q_index]+qureg.stateVec.imag[q_index]*qureg.stateVec.imag[q_index]);
	}
	return energy;
}

double Accelerator::calc_expectation(ExperimentBuffer* buffer){
	double overlap;
	if(ansatz.circuit.qaoa_ansatz || ansatz.circuit.cm_qaoa_ansatz){
		loge("UNIMPLEMENTED");
		throw;
	}
	run_circuit(ansatz.circuit);
	return _energy_evaluation(&overlap, 0);
}

double Accelerator::calc_expectation(ExperimentBuffer* buffer, const std::vector<double> &x, int iteration_i, double* ground_state_overlap_out){
	//int x_size = x.size();
	if(env.numRanks > 1)
		throw_runtime_error("UNIMPLEMENTED");

	if(ansatz.circuit.qaoa_ansatz && ansatz.circuit.cm_qaoa_ansatz)
		throw_runtime_error("ansatz.circuit.qaoa_ansatz && ansatz.circuit.cm_qaoa_ansatz cannot be both true");

	if(ansatz.circuit.qaoa_ansatz || ansatz.circuit.cm_qaoa_ansatz){

		if(x.size() != (unsigned)ansatz.num_params)
			throw_runtime_error("Wrong number of parameters");

		int p = ansatz.num_params/2;

		//logd("Initializing plus state", this->log_level); //too verbosive
		initPlusState(qureg);

		//PRINT
		/*for(int j = 0; j < qureg.numAmpsTotal; ++j)
			std::cerr<<qureg.stateVec.real[j]<<"+"<<qureg.stateVec.imag[j]<<"i ";
		std::cerr<<std::endl;*/
		//PRINT

		for(int i = 0; i < p; ++i){

			//applyTrotterCircuit(qureg, hamiltonian,	x[2*i], 1, 1);
			qreal gamma = x[2*i];

			//std::cerr<<std::setprecision (15)<<"gamma: "<<gamma<<std::endl;
			for(long long i = 0; i < qureg.numAmpsTotal; ++i){
				qreal h = hamDiag.real[i]; //we know hamDiag is real

				qreal a = cos(gamma*h);
				qreal b = -sin(gamma*h);
				qreal c = qureg.stateVec.real[i];
				qreal d = qureg.stateVec.imag[i];

				//std::cerr<<std::setprecision (15)<<"h: "<<h<<" a: "<<a<<" b: "<<b<<" c: "<<c<<" d: "<<d<<std::endl;

				qureg.stateVec.real[i] = a*c-b*d;
				qureg.stateVec.imag[i] = b*c+a*d;

			}

			//PRINT
			/*std::cerr<<"after Hc: ";
			for(int j = 0; j < qureg.numAmpsTotal; ++j)
				std::cerr<<qureg.stateVec.real[j]<<"+"<<qureg.stateVec.imag[j]<<"i ";
			std::cerr<<std::endl;*/
			//PRINT

			if(ansatz.circuit.qaoa_ansatz){
				for(int j = 0; j < qureg.numQubitsInStateVec; ++j)
					rotateX(qureg, j, 2*x[2*i+1]);
				//multiRotatePauli(qureg, qubits_list, all_x_list, qureg.numQubitsInStateVec, -2*x[2*i+1]);
			}else{ //ansatz.circuit.cm_qaoa_ansatz

				std::string constraint_bits = std::bitset<max_qubits>(this->zero_index).to_string();
				constraint_bits=constraint_bits.substr(max_qubits-qureg.numQubitsInStateVec, qureg.numQubitsInStateVec);
				//std::reverse(constraint_bits.begin(), constraint_bits.end());

				for(int j = 0; j < qureg.numQubitsInStateVec; ++j){

					if(constraint_bits[constraint_bits.size()-1-j] == '1')
						pauliX(qureg, j);

					controlledRotateX(qureg, j, (j+1) % qureg.numQubitsInStateVec,2*x[2*i+1]);

					if(constraint_bits[constraint_bits.size()-1-j] == '1')
						pauliX(qureg, j);
				}
				//PRINT
				/*std::cerr<<"after Hm: ";
				for(int j = 0; j < qureg.numAmpsTotal; ++j)
					std::cerr<<qureg.stateVec.real[j]<<"+"<<qureg.stateVec.imag[j]<<"i ";
				std::cerr<<std::endl;;*/
				//PRINT
			}
		}


	}else{

		initZeroState(qureg);
		run_with_new_params(ansatz.circuit, x);
	}
	return _energy_evaluation(ground_state_overlap_out, iteration_i);
}

void Accelerator::__initialize(int num_qubits){

	if(options.createQuregAtEachInilization){
		logd("Initializing " + std::to_string(num_qubits) + " qubits", this->log_level);
		this->env = createQuESTEnv();
		unsigned long int keys[1];
		keys[0] = 1997;
		seedQuEST(&env, keys, 1);
		logd("Setting seed to " + std::to_string(keys[0]), this->log_level);
		this->qureg = createQureg(num_qubits, env);
	}else{
		if(this->qureg.numQubitsInStateVec != num_qubits)
			throw_runtime_error("Different number of qubits requested that those represented by qureg!");

		logd("Skipping qureg initialization. Be sure you know what you're doing!", this->log_level);
	}

}

void Accelerator::initialize(CostFunction cost_function, int num_qubits){

	logd("Calculating hamiltonian terms explicitly.", this->log_level);

	this->hamiltonian_specified = false;
	this->__initialize(num_qubits);

	ref_hamil_energies.clear();
	std::vector<long long unsigned int> indexes(qureg.numAmpsPerChunk);
	std::iota(indexes.begin(), indexes.end(), 0); //zip with indices
	std::sort(indexes.begin(), indexes.end(), [&](int i, int j){return cost_function(i) < cost_function(j);}); //non-descending

	throw_runtime_error("Unimplemented! RefEnergy(cost_function(index), index, false) should replace false with a false/true marking a solution!");

	for(auto &index : indexes){
		ref_hamil_energies.push_back(RefEnergy(cost_function(index), index, false));
		//std::cerr<<"."<<cost_function(index)<<" "<<index<<std::endl;
	}


}

void Accelerator::initialize(PauliHamiltonian* hamIn, bool use_external_hamDiag, DiagonalOpDuplicate *diagonalOpDuplicate){

	logd("Accelerator initialize", this->log_level);

	this->hamiltonian_specified = true;

	int num_qubits = hamIn->nbQubits;
	this->__initialize(num_qubits);

	int coeffsSize = hamIn->coeffs.size();

	/*if(pauliHamilInitialized){
		destroyPauliHamil(pauliHamiltonian);
		pauliHamilInitialized = false;
	}*/

	/*if(diag && hamDiagInitialized){
		destroyDiagonalOp(hamDiag, env);
		hamDiagInitialized = false;
	}*/


	PauliHamil pauliHamiltonian;
	pauliHamiltonian.numQubits = num_qubits;
	pauliHamiltonian.numSumTerms = coeffsSize;
	//h.termCoeffs = malloc(coeffsSize * sizeof *h.termCoeffs);
	//h.pauliCodes = malloc(num_qubits*coeffsSize * sizeof *h.pauliCodes);
    //for(int i=0; i<num_qubits*coeffsSize; i++)
    //	h.pauliCodes[i] = PAULI_I;

	pauliHamiltonian.termCoeffs = &hamIn->coeffs[0]; //conversion to c array
	pauliHamiltonian.pauliCodes = (enum pauliOpType*)(&hamIn->pauliOpts[0]);

	for(auto &c:hamIn->pauliOpts){ //check for diagonal hamiltonian
		if(c == 1 || c == 2){
			throw_runtime_error("TODO: NON DIAGONAL HAMILTONIAN NOT IMPLEMENTED");
		}
	}

	if(hamDiag.numQubits != num_qubits){
		hamDiag.numElemsPerChunk = (1LL << num_qubits) / env.numRanks;
		hamDiag.chunkId = env.rank;
		hamDiag.numChunks = env.numRanks;
		hamDiag.hermitian = 1;
		if(hamDiag.numQubits > 0)
			free(hamDiag.real);
		hamDiag.real = (qreal*) calloc(hamDiag.numElemsPerChunk, sizeof(qreal));
		hamDiag.numQubits = num_qubits;
	}

	if(use_external_hamDiag && diagonalOpDuplicate != nullptr){

		if(diagonalOpDuplicate->numQubits == -1){
			initDiagonalOpFromPauliHamil(hamDiag, pauliHamiltonian);
			diagonalOpDuplicate->numQubits = num_qubits;
			diagonalOpDuplicate->real = std::vector<qreal>(hamDiag.real, hamDiag.real + hamDiag.numElemsPerChunk);
		}else{
			std::copy(diagonalOpDuplicate->real.begin(), diagonalOpDuplicate->real.end(), hamDiag.real);
		}

	}else{

		//if(diagonalOpDuplicate != nullptr){
			//diagonalOpDuplicate->numQubits = num_qubits;
			//diagonalOpDuplicate->real = std::vector<qreal>(hamDiag.real, hamDiag.real + hamDiag.numElemsPerChunk);
		//}else{
			initDiagonalOpFromPauliHamil(hamDiag, pauliHamiltonian);
		//}

	}

	logd("PauliHamiltonian initialized", this->log_level);

	if(env.numRanks>1){
		throw_runtime_error("TODO: DISTRIBUTED COMPUTATION UNIMPLEMENTED");
	}

	ref_hamil_energies.clear();

	std::vector<long long unsigned int> indexes(qureg.numAmpsPerChunk);
	std::iota(indexes.begin(), indexes.end(), 0); //zip with indices
	std::sort(indexes.begin(), indexes.end(), [&](int i, int j){return hamDiag.real[i] < hamDiag.real[j];}); //non-descending

	int counter = 0;
	if(options.zero_reference_states.size() > 1)
		loge("TODO: Add more zero reference states compatibility");

	if(hamIn -> custom_solutions.size() == 0)
		logw("Unimplemented! TODO: If no custom solutions are provided, choose the lowest energy one.", this->log_level);

	for(auto &index : indexes){

		//TODO: add more zero_reference states
		if(options.zero_reference_states.size() > 1 && index == options.zero_reference_states[0]){
			//logw("Zero excluded with counter " + std::to_string(counter) + " where E(0) = " + std::to_string(hamDiag.real[index]));
			if(hamDiag.real[index]!=0){
				loge("Tried to exclude something else than zero ground state! Did you run qaoa?");
			}
			else{
				continue;
			}
		}counter++;


		//logw(std::to_string(index)+"       " + std::to_string(hamDiag.real[index]));

		if(hamDiag.real[index] == 0){
			logd("Here we have not exluded zero at index " + std::to_string(index), this->log_level);
		}
		//std::cerr<<index<<" "<<hamDiag.real[index]<<"\n";

		bool isConsideredSolution=false;

		for(const auto& s: hamIn->custom_solutions){
			std::string index_b = std::bitset<max_qubits>(index).to_string();
			index_b=index_b.substr(max_qubits-num_qubits, num_qubits);

			if(index_b.size() != s.size())
				throw_runtime_error("index_b.size != s.size()");

			bool matches=true;
			for(unsigned int i = 0; i < index_b.size(); ++i){
				if(s[i] == '.')
					continue;
				if(s[i] != index_b[i]){
					matches = false;
					break;
				}
			}

			if(matches){
				isConsideredSolution = true;
				break;
			}

		}
		ref_hamil_energies.push_back(RefEnergy(hamDiag.real[index], index, isConsideredSolution));

		//if( double(counter++)/indexes.size() > options.samples_cut_ratio)
		//	break;
	}
	ref_hamil_energies_set = true;

}

void Accelerator::initialize(int num_qubits){
	logd("Initializing qureg with " + std::to_string(num_qubits) + " qubits", this->log_level);
	this->__initialize(num_qubits);
}

void Accelerator::run_vqe_slave_process(){
	logw("Slave process not implemented in non-distributed version!", this->log_level);
}

void Accelerator::finalize(){
	logd("Destroying qureg", this->log_level);
	if(qureg.numQubitsInStateVec != -1)
		destroyQureg(qureg, env);
	qureg.numQubitsInStateVec = -1;
}

Accelerator::Accelerator(AcceleratorOptions options){

	this->log_level = options.log_level;

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
	if(!options.createQuregAtEachInilization){

		if(options.createQuregAtEachInilization_num_qubits < 0){
			loge("If createQuregAtEachInilization=false, need to specify number of qubits by setting createQuregAtEachInilization_num_qubits");
			throw;
		}
		logd("Creating QuEST environment", this->log_level);
		this->env = createQuESTEnv();

		logd("Initializing " + std::to_string(options.createQuregAtEachInilization_num_qubits) + " qubits", this->log_level);
		this->qureg = createQureg(options.createQuregAtEachInilization_num_qubits, env);
	}

	hamDiag.numQubits = 0;

}
}
