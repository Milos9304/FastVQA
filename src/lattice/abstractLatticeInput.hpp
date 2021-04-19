/*
 * abstractLatticeInput.hpp
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTICE_ABSTRACTLATTICEINPUT_HPP_
#define SRC_LATTICE_ABSTRACTLATTICEINPUT_HPP_

#include <string>

class AbstractLatticeInput{

	public:
		std::string name;
		int n;

		virtual std::string toHamiltonianString() = 0;
		virtual ~AbstractLatticeInput(){}

};



#endif /* SRC_LATTICE_ABSTRACTLATTICEINPUT_HPP_ */
