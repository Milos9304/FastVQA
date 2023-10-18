/*
 * vqeOptions.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_VQE_VQEOPTIONS_H_
#define SRC_VQE_VQEOPTIONS_H_

#include "vqaOptions.h"
namespace FastVQA{
class VQEOptions : public VQAOptions{

	public:
		// Name of the ansatz as defined in FastVQA/ansatz.h
		std::string ansatz_name = "Ry_CNOT_all2all_Ry";

};
}
#endif /* SRC_VQE_VQEOPTIONS_H_ */
