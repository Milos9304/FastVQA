/*
 * AqcPqc.h
 *
 *  Created on: June 17, 2022
 *      Author: Milos Prokop
 */

#ifndef AQC_PQC_H_
#define AQC_PQC_H_

#include "../ansatz.h"
#include "AqcPqcOptions.h"

namespace FastVQA{

class AqcPqc{

	public:

		void run_aqpc(Hamiltonian* h_target, AqcPqcOptions* options);

};

}

#endif /* AQC_PQC_H_ */
