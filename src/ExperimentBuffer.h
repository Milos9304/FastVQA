/*
 * ExperimentBuffer.h
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_EXPERIMENT_BUFFER_H_
#define SRC_EXPERIMENT_BUFFER_H_

class ExperimentBuffer{
public:
   double opt_val;
   std::vector<double> opt_params;
   std::string opt_config;
   double expected_energy;
   double hit_rate;
};

#endif /* SRC_EXPERIMENT_BUFFER_H_ */
