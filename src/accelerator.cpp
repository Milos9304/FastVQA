#include "accelerator.h"
#include "logger.h"
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cfloat>

# if QuEST_PREC==1
	#define QREAL_MAX FLOAT_MAX
# elif QuEST_PREC==2
    #define QREAL_MAX DOUBLE_MAX
# elif QuEST_PREC==4
   #define QREAL_MAX LONG_DOUBLE_MAX
# endif

namespace FastVQA{

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

	int min_sol=std::get<0>(ref_hamil_energies[0]);
	int i = 0;
	while(ref_hamil_energies[i].first == min_sol){
		res.push_back(ref_hamil_energies[i]);
		i++;
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

	//bool second_eigenergy = true;

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

	if(!this->options.exclude_zero_state && ref_hamil_energies[0].first == 0)
		loge("Zero not excluded properly!");

	std::vector<RefEnergy> ground_states;

	//if(this->options.exclude_zero_state){
		long long int i = 0;
		qreal min_gs = DBL_MAX;
		for(auto &hamil_energy : ref_hamil_energies){
			if(hamil_energy.first != 0 && hamil_energy.first < min_gs)
				min_gs = hamil_energy.first;
		}
		for(auto &hamil_energy : ref_hamil_energies){
			if(hamil_energy.first == min_gs)
				ground_states.push_back(hamil_energy);
		}
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

	const int max_qubits=40;
	if(qureg.numQubitsInStateVec > max_qubits){
		loge("Not enough binary places to hold final opt_config. You are simulating more than" + std::to_string(max_qubits) + " qubits. Increase the number in the code.");
	}


	for(auto &ground_state: ground_states){
		i = ground_state.second;
		std::string opt_config = std::bitset<max_qubits>(i).to_string();
		opt_config=opt_config.substr(max_qubits-qureg.numQubitsInStateVec,qureg.numQubitsInStateVec);
		std::reverse(opt_config.begin(), opt_config.end());

		ExperimentBuffer::ExperimentBufferSolution solution(opt_config, ground_state.first, qureg.stateVec.real[i]*qureg.stateVec.real[i]+qureg.stateVec.imag[i]*qureg.stateVec.imag[i]);
		buffer->final_solutions.push_back(solution);
	}

	buffer->opt_val=ground_states[0].first;
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

	if(!this->options.exclude_zero_state && ref_hamil_energies[0].first == 0)
			loge("Zero was not excluded properly!");

	long long int ground_index = ref_hamil_energies[0].second;
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

	double alpha = alpha_f(options.samples_cut_ratio, options.final_alpha, iteration_i, options.max_alpha_iters);

	while(alpha_sum < alpha){

		 if(i >= ref_hamil_energies.size()){
			 break;
		 }

		long long int q_index = ref_hamil_energies[i].second;
		amp = zero_state_amp_per_elem + qureg.stateVec.real[q_index]*qureg.stateVec.real[q_index]+qureg.stateVec.imag[q_index]*qureg.stateVec.imag[q_index];
		alpha_sum += amp;
		if(i >= ref_hamil_energies.size()){
			loge("NOT ENOUGH REF ENERGIES STORED IN MEMORY.");
		}
		else{
			energy += ref_hamil_energies[i].first * (amp / alpha);
		}
		i++;
	}
	energy -= (alpha_sum - alpha) * ref_hamil_energies[i-1].first * (amp / alpha);
	return energy;
}

double Accelerator::calc_expectation(ExperimentBuffer* buffer){
	double overlap;
	if(ansatz.circuit.qaoa_ansatz){
		loge("UNIMPLEMENTED");
		throw;
	}
	run_circuit(ansatz.circuit);
	return _energy_evaluation(&overlap, 0);
}

double Accelerator::calc_expectation(ExperimentBuffer* buffer, const std::vector<double> &x, int iteration_i, double* ground_state_overlap_out){
	//int x_size = x.size();
	if(env.numRanks > 1){
		loge("UNIMPLEMENTED");
		throw;
	}

	if(ansatz.circuit.qaoa_ansatz){

		if(x.size() != (unsigned)ansatz.num_params){
			loge("Wrong number of parameters");
		}

		int p = ansatz.num_params/2;

		//logd("Initializing plus state", this->log_level); //too verbosive
		initPlusState(qureg);

		for(int i = 0; i < p; ++i){

			//applyTrotterCircuit(qureg, hamiltonian,	x[2*i], 1, 1);
			qreal gamma = x[2*i];
			for(long long i = 0; i < qureg.numAmpsTotal; ++i){
				qreal h = hamDiag.real[i]; //we know hamDiag is real

				qreal a = cos(gamma*h);
				qreal b = -sin(gamma*h);
				qreal c = qureg.stateVec.real[i];
				qreal d = qureg.stateVec.imag[i];

				qureg.stateVec.real[i] = a*c-b*d;
				qureg.stateVec.imag[i] = b*c+a*d;

			}

			for(int j = 0; j < qureg.numQubitsInStateVec; ++j)
				rotateX(qureg, j, 2*x[2*i+1]);
			//multiRotatePauli(qureg, qubits_list, all_x_list, qureg.numQubitsInStateVec, -2*x[2*i+1]);
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

	for(auto &index : indexes){
		ref_hamil_energies.push_back(RefEnergy(cost_function(index), index));
		//std::cerr<<"."<<cost_function(index)<<" "<<index<<std::endl;
	}


}

void Accelerator::initialize(PauliHamiltonian* hamIn){

	bool diag=true;

	logd("Calculating Hamiltonian terms explicitly.", this->log_level);

	this->hamiltonian_specified = true;

	int num_qubits = hamIn->nbQubits;
	this->__initialize(num_qubits);

	qubits_list = new int[num_qubits]();
	all_x_list = new pauliOpType[num_qubits]();

	for(int i = 0; i < num_qubits; ++i){
		qubits_list[i]=i;
		all_x_list[i]=PAULI_X;
	}

	int coeffsSize = hamIn->coeffs.size();

	if(pauliHamilInitialized){
		destroyPauliHamil(pauliHamiltonian);
		if(diag){
			destroyDiagonalOp(hamDiag, env);
		}
	}

	pauliHamiltonian = createPauliHamil(num_qubits, coeffsSize);

	pauliHamiltonian.termCoeffs = &hamIn->coeffs[0]; //conversion to c array
	pauliHamiltonian.pauliCodes = (enum pauliOpType*)(&hamIn->pauliOpts[0]);


	for(auto &c:hamIn->pauliOpts){ //check for diagonal hamiltonian
		if(c == 1 || c == 2){
			diag=false;
			break;
		}
	}

	if(diag){
		hamDiag = createDiagonalOp(num_qubits, env, 1);
		initDiagonalOpFromPauliHamil(hamDiag, pauliHamiltonian);
	}else{ //assuming PauliHamil
		throw_runtime_error("TODO: NON DIAGONAL HAMILTONIAN NOT IMPLEMENTED");
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
			logw("Here we have not exluded zero at index " + std::to_string(index), this->log_level);
		}
		//std::cerr<<index<<" "<<hamDiag.real[index]<<"\n";
		ref_hamil_energies.push_back(RefEnergy(hamDiag.real[index], index));
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
	destroyQureg(qureg, env);
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


}
}
