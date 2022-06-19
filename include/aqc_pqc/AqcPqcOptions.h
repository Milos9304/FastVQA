/*
 * AqcPqcOptions.h
 *
 *  Created on: June 17, 2022
 *      Author: Milos Prokop
 */

#ifndef AQC_PQC_OPTIONS_H_
#define AQC_PQC_OPTIONS_H_

#include "../ansatz.h"
#include "AqcPqcOptions.h"

namespace FastVQA{

class AqcPqcOptions{

	public:

		// Name of the initial hamiltonian as defined in FastVQA/ansatz.h
		std::string ansatz_name = "Ry_CNOT_all2all_Ry";

};

}

#endif /* AQC_PQC_OPTIONS_H_ */



