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



	long long int i = 0;
	while(ref_hamil_energies[i++].first == 0);

	RefEnergy ground_state = ref_hamil_energies[i-1];
	i = ground_state.second;

	/*
		 *
		 */


		std::ofstream output_file("../experiment_files/kokotina.txt");
		output_file << std::fixed << std::showpoint;
		output_file << std::setprecision(10);
		for(long long int i = 0; i < qureg.numAmpsTotal; ++i){
			output_file << qureg.stateVec.real[i]*qureg.stateVec.real[i]+qureg.stateVec.imag[i]*qureg.stateVec.imag[i] << "\n";
		}

		output_file << i << "\n";

		output_file.close();
		loge("kokotina hotova");

		//throw;

		/*
		 *
		 */

		//loge("Delete fstream");




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



	long long int ground_index = ref_hamil_energies[0].second;
	if(ref_hamil_energies[0].first == 0)
		loge("Zero not excluded properly!");

	*ground_state_overlap_out = qureg.stateVec.real[ground_index]*qureg.stateVec.real[ground_index]+qureg.stateVec.imag[ground_index]*qureg.stateVec.imag[ground_index];



	//loge("Modified calc expec");
	return calcExpecDiagonalOp(qureg, hamDiag).real;

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

	long long unsigned int i = 0;

	std::vector<long long unsigned int> zero_state_indices = options.zero_reference_states;

	if(zero_state_indices.size() > 1){
		loge("TODO: Multiple zero_state_indices not implemented yet");
	}long long unsigned int zero_state_index = zero_state_indices[0];

	double zero_state_amp;
	if(ansatz.circuit.qaoa_ansatz)
		zero_state_amp = 0;
	else
		zero_state_amp = qureg.stateVec.real[zero_state_index]*qureg.stateVec.real[zero_state_index]+qureg.stateVec.imag[zero_state_index]*qureg.stateVec.imag[zero_state_index];

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

	int log_level = options.log_level;
	int num_qubits = hamIn->nbQubits;

	qubits_list = new int[num_qubits]();
	all_x_list = new pauliOpType[num_qubits]();

	for(int i = 0; i < num_qubits; ++i){
		qubits_list[i]=i;
		all_x_list[i]=PAULI_X;
	}

	logd("Initializing " + std::to_string(num_qubits) + " qubits", log_level);

	unsigned long int keys[1];
	keys[0] = 1997;
	seedQuEST(&env, keys, 1);
	logd("Setting seed to " + std::to_string(keys[0]));

	qureg = createQureg(num_qubits, env);

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
}
}
