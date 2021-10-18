/*
 * qaoaOptions.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_QAOA_QAOAOPTIONS_H_
#define SRC_QAOA_QAOAOPTIONS_H_

#include "../vqaOptions.h"

class QAOAOptions : public VQAOptions{

	public:
		int p=1;
		bool extendedParametrizedMode = false;

		std::string getParameterScheme(){
			return extendedParametrizedMode ? "Extended" : "Standard";
		}

};

#endif /* SRC_QAOA_QAOAOPTIONS_H_ */
