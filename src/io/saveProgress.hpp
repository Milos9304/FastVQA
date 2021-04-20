/*
 * saveProgress.hpp
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_IO_SAVEPROGRESS_HPP_
#define SRC_IO_SAVEPROGRESS_HPP_

#include <vector>
#include <string>

/*
#include <string>
#include <fstream>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
*/
class save_instance{

	private:
		/*friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version){
		        ar & coefficients;
		        ar & expected_energy;
		        ar & sv;
		        ar & sv_energy;
		        ar & hit_rate;
		}
*/
	public:
		save_instance();
		save_instance(std::vector<double> coefficients, double expected_energy, double sv_energy, double hit_rate){
			this->coefficients=coefficients;
			this->expected_energy=expected_energy;
			this->sv_energy=sv_energy;
			this->hit_rate=hit_rate;
		}

		save_instance(std::vector<double> coefficients, double expected_energy, std::string sv, double sv_energy, double hit_rate){
					this->coefficients=coefficients;
					this->expected_energy=expected_energy;
					this->sv = sv;
					this->sv_energy=sv_energy;
					this->hit_rate=hit_rate;
		}


		std::vector<double> coefficients;
		std::string sv=""; //found shortest vector
		double expected_energy;
		double sv_energy;
		double hit_rate;

};
void saveProgress(std::string filename, std::vector<double> coefficients, double expected_energy, double sv_energy, double hit_rate);
void saveProgress(std::string filename, std::vector<double> coefficients, double expected_energy, std::string sv, double sv_energy, double hit_rate);
bool loadProgress(std::string filename, std::vector<double>* coefficients, double* expected_energy, double* sv_energy, double* hit_rate);



#endif /* SRC_IO_SAVEPROGRESS_HPP_ */
