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

	typedef std::pair<int, double> meas_freq_eval;
	typedef std::pair<std::string, meas_freq_eval> measurement;

	std::vector<measurement> measurements;

	Qureg qureg_cache = createQureg(qureg.numQubitsInStateVec, env);

	std::string classicalRefState_str(qureg.numQubitsRepresented, '0');
	classicalRefState_str[qureg.numQubitsRepresented-1]='1';

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
		if(/*!overlapPenalization || */measurementStr != classicalRefState_str){
			if(!found){
			  measurements.push_back(measurement(measurementStr, // bit string
									  meas_freq_eval(1, //frequency
											  evaluate_assignment(hamiltonian, measurementStr))));
			}else{
				measurements[index].second.first++;
			}
		}
	}
	if(/*m_maximize*/false)
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

	measurement optimalMeasurement;
	if(second_eigenergy && measurements.size() > 1 && measurements[0].second.second == 0.){
		optimalMeasurement = measurements[1];
		logd("Returning second highest energy level");
	}else
		optimalMeasurement = measurements[0];

	int i = 0;
	int hits = 0;
	while(i < measurements.size() && abs(measurements[i++].second.second - optimalMeasurement.second.second)<10e-1){
	  hits += measurements[i-1].second.first;			}

	loge("Measurements size = " + std::to_string(measurements.size()));

	buffer->opt_config=optimalMeasurement.first;
	buffer->opt_val=optimalMeasurement.second.second;
	buffer->hit_rate=hits / double(nbSamples);

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

double Accelerator::calc_expectation(ExperimentBuffer* buffer, const std::vector<double> &x){

	int x_size = x.size();

	for(int i = 1; i < env.numRanks; ++i){
		MPI_Send(&x_size, 1, MPI_INT, i, calc_exp_val_tag , MPI_COMM_WORLD);
	}

	std::vector<double> x_copy(x);

	MPI_Bcast(&x_copy[0], x_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	initZeroState(qureg);

	run(ansatz.circuit, x);
	return calcExpecDiagonalOp(qureg, hamDiag).real;

}

void Accelerator::initialize(Hamiltonian* hamIn){

	int num_qubits = hamIn->nbQubits;

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

Accelerator::Accelerator(std::string accelerator){

	if(accelerator != "quest"){
		loge("No other accelerator than QuEST implemented");
		throw;
	}

}
