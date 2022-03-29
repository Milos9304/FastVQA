/*
 * gate.h
 *
 *  Created on: Apr 26, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_GATE_H_
#define FASTVQA_GATE_H_

#include <memory>
#include <string>
namespace fastVQA{

typedef unsigned char GateCode;

class Parameter{
public:
	static const std::string blank_param_name;

	std::string name;
	double value;
	Parameter(std::string name, double value){
		this->name = name;
		this->value = value;
	}
};

class Gate{

private:
	static const unsigned char two_qubit=128;
	static const unsigned char one_qubit=0;

public:

	static const std::shared_ptr<Parameter> blankParam;

	Gate(GateCode code, int qubit1, std::shared_ptr<Parameter> param=blankParam);

	Gate(GateCode code, int qubit1, int qubit2, std::shared_ptr<Parameter> param=blankParam);

	GateCode code;
	bool is_twoQubit(){return code & two_qubit;}
	int qubit1=-1, qubit2=-1;
	std::shared_ptr<Parameter> param;

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
}
#endif /* FASTVQA_GATE_H_ */
