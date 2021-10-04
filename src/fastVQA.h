/*
 * fastVQA.h
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LATTIQ_H_
#define SRC_LATTIQ_H_

using namespace indicators;

Color colors[9] = {Color::red, Color::green, Color::yellow, Color::blue, Color::magenta, Color::cyan};

#define bar_opts(counter, num_lattices, lattice_name)   option::BarWidth{50},\
														option::Start{"["},\
														option::Fill{"="},\
														option::Lead{">"},\
														option::Remainder{" "},\
														option::End{"]"},\
														option::PrefixText{std::to_string(counter+1) + "/" + std::to_string(num_lattices) + " " + lattice_name},\
														option::ForegroundColor{colors[counter % 7]},\
														option::ShowElapsedTime{true},\
														option::ShowRemainingTime{true},\
														option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},\
														option::MaxProgress{qaoaOptions.max_iters}


std::string process(std::string const& s)
{
    std::string::size_type pos = s.find('_');
    if (pos != std::string::npos)
    {
        return s.substr(0, pos);
    }
    else
    {
        return s;
    }
}


#endif /* SRC_LATTIQ_H_ */
