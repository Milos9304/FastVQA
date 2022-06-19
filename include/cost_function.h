/*
 * cost_function.h
 *
 *  Created on: Mar 30, 2022
 *      Author: Milos Prokop
 */

#ifndef SRC_COST_FUNCTION_H_
#define SRC_COST_FUNCTION_H_

#include <functional>

namespace FastVQA{
	/*Input:
		decimal encoding of binary string
		bit_size

	  Output:
		cost value
	*/
	typedef std::function<double(long long int i)> CostFunction;

}

#endif /* SRC_COST_FUNCTION_H_ */
