#include "accelerator.h"
#include "logger.h"
#include <iomanip>
#include <fstream>
#include <algorithm>

namespace fastVQA{

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

	RefEnergy ground_state;

	long long int i = 0;
	if(this->options.exclude_zero_state){
		while(ref_hamil_energies[i++].first == 0);
		ground_state = ref_hamil_energies[i-1];
	}else{

		if(this->options.choose_ground_state_with_smallest_index){


			/*for(int j = 0; j < ref_hamil_energies.size(); ++j){
				std::cout << ref_hamil_energies[j].first << " " << ref_hamil_energies[j].second << std::endl;
			}

			double p=0;
			for(long long j = 0; j < qureg.numAmpsTotal; ++j){
				p+=qureg.stateVec.real[j]*qureg.stateVec.real[j]+qureg.stateVec.imag[j]*qureg.stateVec.imag[j];
			}

			logw("Qreg prob=" + std::to_string(p));

			p=0;*/

			//double p = 0;
			int j = 0;
			while(ref_hamil_energies[j++].first == ref_hamil_energies[0].first){
				//long long index = ref_hamil_energies[j-1].second;
				//p+=qureg.stateVec.real[index]*qureg.stateVec.real[index]+qureg.stateVec.imag[index]*qureg.stateVec.imag[index];
			}
			//loge("Total prob=" + std::to_string(p));

			//std::pair<qreal, long long int>
			std::sort(ref_hamil_energies.begin(), ref_hamil_energies.begin() + (j-1),
					[](const RefEnergy& a, const RefEnergy& b) {
						return a.second < b.second;
			});

			ground_state = ref_hamil_energies[0];
			std::cerr<<ground_state.first<<" "<<ground_state.second<<"\n";
			//loge("I choose: " + std::to_string(ground_state.second));

		}else
			ground_state = ref_hamil_energies[i];

	}
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

	//int circuit_size = this->ansatz.circuit.gates.size();

	std::vector<double> gateCodes;
	for(auto &gate : this->ansatz.circuit.gates){
		gateCodes.push_back(gate.code);
		gateCodes.push_back(gate.qubit1);
		gateCodes.push_back(gate.qubit2);
	}
}

double Accelerator::calc_expectation(ExperimentBuffer* buffer, const std::vector<double> &x, int iteration_i, double* ground_state_overlap_out){
	//int x_size = x.size();
	if(env.numRanks > 1){
		loge("UNIMPLEMENTED");
		throw;
	}

	std::vector<double> x_copy(x);

	if(ansatz.circuit.qaoa_ansatz){

		if(x.size() != (unsigned)ansatz.num_params){
			loge("Wrong number of parameters");
		}

		int p = ansatz.num_params/2;

		initPlusState(qureg);

		for(int i = 0; i < p; ++i){

			//applyTrotterCircuit(qureg, hamiltonian,	x[2*i], 1, 1);
			qreal gamma = x[2*i];
			for(long long i = 0; i < qureg.numAmpsTotal; ++i){
				qreal h = hamDiag.real[i];

				qreal a = cos(gamma*h);
				qreal b = -sin(gamma*h);
				qreal c = qureg.stateVec.real[i];
				qreal d = qureg.stateVec.imag[i];

				qureg.stateVec.real[i] = a*c-b*d;
				qureg.stateVec.imag[i] = b*c+a*d;

			}
			multiRotatePauli(qureg, qubits_list, all_x_list, qureg.numQubitsInStateVec, x[2*i+1]);
		}

	}else{
		initZeroState(qureg);
		run(ansatz.circuit, x);
	}

	if(this->options.exclude_zero_state && ref_hamil_energies[0].first == 0)
		loge("Zero not excluded properly!");

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

void Accelerator::__initialize(int num_qubits){

	unsigned long int keys[1];
	keys[0] = 1997;
	seedQuEST(&env, keys, 1);
	logd("Setting seed to " + std::to_string(keys[0]), options.log_level);

	if(options.createQuregAtEachInilization){
		logd("Initializing " + std::to_string(num_qubits) + " qubits", options.log_level);
		this->qureg = createQureg(num_qubits, env);
	}else{
		logw("Skipping qureg initialization. Be sure you know what you're doing!", options.log_level);
	}

}

void Accelerator::initialize(CostFunction cost_function, int num_qubits){

	logd("Calculating hamiltonian terms explicitly.", options.log_level);

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

void Accelerator::initialize(Hamiltonian* hamIn){

	logd("Calculating hamiltonian terms explicitly.", options.log_level);

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

	hamiltonian = createPauliHamil(num_qubits, coeffsSize);

	hamiltonian.termCoeffs = &hamIn->coeffs[0]; //conversion to c array
	hamiltonian.pauliCodes = (enum pauliOpType*)(&hamIn->pauliOpts[0]);

	hamDiag = createDiagonalOp(num_qubits, env, 1);
	initDiagonalOpFromPauliHamil(hamDiag, hamiltonian);

	if(env.numRanks>1){
		throw_runtime_error("TODO: UNIMPLEMENTED");
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
		if(index == options.zero_reference_states[0]){
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
			loge("Here we have not exluded zero");
		}else
			ref_hamil_energies.push_back(RefEnergy(hamDiag.real[index], index));
		//if( double(counter++)/indexes.size() > options.samples_cut_ratio)
		//	break;
	}
}

void Accelerator::run_vqe_slave_process(){
	logw("Slave process not implemented in non-distributed version!");
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
	if(!options.createQuregAtEachInilization){

		if(options.createQuregAtEachInilization_num_qubits < 0){
			loge("If createQuregAtEachInilization=false, need to specify number of qubits by setting createQuregAtEachInilization_num_qubits");
			throw;
		}
		logd("Creating QuEST environment");
		this->env = createQuESTEnv();

		logd("Initializing " + std::to_string(options.createQuregAtEachInilization_num_qubits) + " qubits", options.log_level);
		this->qureg = createQureg(options.createQuregAtEachInilization_num_qubits, env);
	}


}
}
