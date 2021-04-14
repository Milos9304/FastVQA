/*
 * utils.hpp
 *
 *  Created on: Apr 14, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_UTILS_HPP_
#define SRC_UTILS_HPP_

class OnlineMeanVarianceCalculator{

	public:

		double mean = 0;
		double variance = 0;

		void update(double sample);

		void reset(){
			mean = variance = 0;
		}

	private:
		int n = 0;
};



#endif /* SRC_UTILS_HPP_ */
