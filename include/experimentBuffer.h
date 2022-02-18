/*
 * experimentBuffer.h
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#ifndef FASTVQA_EXPERIMENT_BUFFER_H_
#define FASTVQA_EXPERIMENT_BUFFER_H_

namespace fastVQA{

class ExperimentBuffer{
public:
   double opt_val;
   std::vector<double> opt_params;
   std::string opt_config;
   double expected_energy;
   double hit_rate;

   std::vector<double> intermediateEnergies;
   std::vector<double> intermediateGroundStateOverlaps;
};
}
#endif /* FASTVQA_EXPERIMENT_BUFFER_H_ */
