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


double Accelerator::calc_expectation(ExperimentBuffer* buffer, Ansatz* ansatz, const std::vector<double> &x){

	initZeroState(qureg);
	run(ansatz->circuit, x);
	return calcExpecDiagonalOp(qureg, hamDiag).real;

}

void Accelerator::initialize(Hamiltonian* hamIn){

	int num_qubits = hamIn->nbQubits;

	logd("Initializing " + std::to_string(num_qubits) + " qubits");

	unsigned long int keys[1];
	keys[0] = 1997;
	seedQuEST(&env, keys, 1);
	logd("Setting seed to " + std::to_string(keys[0]));

	qureg = createQureg(num_qubits, env);

	hamiltonian = createPauliHamil(num_qubits, hamIn->coeffs.size());
	hamiltonian.termCoeffs = &hamIn->coeffs[0]; //conversion to c array
	hamiltonian.pauliCodes = (enum pauliOpType*)(&hamIn->pauliOpts[0]);

	hamDiag = createDiagonalOp(num_qubits, env, 1);
	initDiagonalOpFromPauliHamil(hamDiag, hamiltonian);

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
