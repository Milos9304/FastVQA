/*
 * qaoa.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_QAOA_QAOA_H_
#define SRC_QAOA_QAOA_H_

#include "PauliOperator.hpp"
#include "../indicators/progress_bar.hpp"


//void run_qaoa(xacc::quantum::PauliOperator, bool verbose);

class QAOAOptions{

	public:
		int max_iters = 0;
		bool verbose = true;

};

void run_qaoa(std::string hamiltonian, std::string name, indicators::ProgressBar* bar, QAOAOptions options);/*{
	run_qaoa(xacc::quantum::PauliOperator(hamiltonian), verbose);
}*/

#endif /* SRC_QAOA_QAOA_H_ */
