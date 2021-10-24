/*
 * gate_encodings.h
 *
 *  Created on: Apr 26, 2021
 *      Author: Milos Prokop
 */

#ifndef QUACC_GATE_ENCODINGS_H_
#define QUACC_GATE_ENCODINGS_H_

typedef unsigned char GateCode;

class Gate{

private:
	static const uint8_t two_qubit=128;
	static const uint8_t one_qubit=0;

public:

	Gate(GateCode code, int qubit1, double param=.0){
		if(code & two_qubit){
		  loge(std::to_string(code) + "is a two qubit gate");
		  throw;
		}

		this->code=code;
		this->qubit1=qubit1;
		this->param=param;
	}

	Gate(GateCode code, int qubit1, int qubit2, double param=.0){
		if(!(code & two_qubit)){
		  loge(std::to_string(code) + "is a one qubit gate");
		  throw;
		}

		this->code=code;
		this->qubit1=qubit1;
		this->qubit2=qubit2;
		this->param=param;
	}

	GateCode code;
	bool is_twoQubit(){return code & two_qubit;}
	int qubit1, qubit2;
	double param;

	//Gate codes
	static const GateCode g_I = 1 | one_qubit;
	static const GateCode g_H = 2 | one_qubit;
	static const GateCode g_CNOT = 3 | two_qubit;
	static const GateCode g_Rz = 4 | one_qubit;
	static const GateCode g_Ry = 5 | one_qubit;
	static const GateCode g_Rx = 6 | one_qubit;
	static const GateCode g_X = 7 | one_qubit;
	static const GateCode g_Y = 8 | one_qubit;
	static const GateCode g_Z = 9 | one_qubit;
	static const GateCode g_CPhase = 10 | two_qubit;
	static const GateCode g_Swap = 11 | two_qubit;
	static const GateCode g_iSwap = 12 | two_qubit;
	//static const GateCode g_fSim = 13;
	static const GateCode g_Measure = 14 | one_qubit;
	static const GateCode g_CZ = 15 | two_qubit;
	static const GateCode g_CY = 16 | two_qubit;
	static const GateCode g_CRZ = 17 | two_qubit;
	static const GateCode g_CH = 18 | two_qubit;
	//static const GateCode g_S = 19;
	//static const GateCode g_Sdg = 20;
	static const GateCode g_T = 21 | two_qubit;
	//static const GateCode g_Tdg = 22;
	static const GateCode g_U = 23 | one_qubit;
	//static const GateCode g_U1 = 24;
	//static const GateCode g_IfStmt = 25;
	//static const GateCode g_XY = 26;
	//static const GateCode g_Reset = 27;
	static const GateCode g_initPlusState = 100;
	static const GateCode g_overlap_instruction = 101;

};

#endif /* QUACC_GATE_ENCODINGS_H_ */
