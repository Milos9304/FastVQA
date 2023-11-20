/*
 * qaoaOptions.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_QAOA_QAOAOPTIONS_H_
#define FASTVQA_QAOA_QAOAOPTIONS_H_

#include "vqaOptions.h"


namespace FastVQA{

class QAOAOptions : public VQAOptions{

	public:
		int p=1;
		bool extendedParametrizedMode = false;

		std::vector<double> initial_params; //if empty, they are generated at random

		//TODO: explain
		bool simplifiedSimulation = true;

		std::string getParameterScheme(){
			return extendedParametrizedMode ? "Extended" : "Standard";
		}

};
}

#endif /* FASTVQA_QAOA_QAOAOPTIONS_H_ */
