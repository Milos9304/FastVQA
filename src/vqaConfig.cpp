/*
 * config_io.cpp
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#include "logger.h"
#include "vqaConfig.h"

MatrixInt VqaConfig::loadLatticeFromFile(std::string filename, bool *success){

	std::ifstream file(filename);
	if(!file.is_open()){
		*success = false;
		loge("Unable to open hamiltonian file " + filename);
		throw std::runtime_error("Could not open file");
	}

	std::vector<std::string> lines;

	std::string line;
	while (getline(file, line)){
		if (!line.empty())
			lines.push_back(line);
	}

	file.close();

	std::string line0 = lines[0];

	int rows = lines.size();
	int cols = 1;
	for (unsigned int i = 0; i < line0.length(); i++) {
		if (isspace(line0.at(i))){
		  cols++;
		}
	}

	MatrixInt lattice;
	lattice.resize(rows, cols);

	int row_index = 0;
	int col_index = 0;

	for(auto &line: lines){

		 if(line[1] == '[')
			 line = line.substr(2, line.length()-1);
		 else if(line[line.length()-2] == ']')
			 line = line.substr(1, line.length()-2);
		 else
			 line = line.substr(1, line.length()-1);

		 std::istringstream in(line, std::istringstream::in);

		 int n;

		 col_index = 0;
		 while (in >> n)
			 lattice(row_index, col_index++) = n;

		 row_index += 1;
	}

	*success = true;
	return lattice;
}

VqaConfig::VqaConfig(std::string pathname){

	std::string lattice_files;

	std::ifstream ifs(pathname);
	std::string line;
	std::istringstream line_stream;

	if(!ifs.is_open()){
		throw std::runtime_error("Unable to open config file " + pathname);
	}

	bool lookForType = true;

	//
	// LOAD CONFIGURATION FILE
	//

	while (std::getline(ifs, line)) {

		if(line[0] == '#')
				continue;

		line_stream.str(line.substr(line.find("=")+1));

		if(lookForType and line.find("type") == std::string::npos)
			throw std::runtime_error("Invalid config file type.");

		lookForType = false;

		if (line.find("type") != std::string::npos) {
			if(line_stream.str() != "qaoa")
				throw std::runtime_error("Invalid config file type.");
		}

		else if (line.find("lattice_files") != std::string::npos) {

			line_stream >> lattice_files;

			std::istringstream ss(lattice_files);
			std::string lattice_file;

			while(std::getline(ss, lattice_file, ',')) {
				bool success;
				MatrixInt lattice = loadLatticeFromFile(lattice_file, &success);

				if(success){

					std::string name = lattice_file;

					const size_t last_slash_idx = name.find_last_of("\\/");
					if (std::string::npos != last_slash_idx)
						name.erase(0, last_slash_idx + 1);

					const size_t period_idx = name.rfind('.');
					if (std::string::npos != period_idx)
						name.erase(period_idx);

					lattices.push_back(Lattice(lattice, name));
				}

				logi(lattice_file + " loaded");
			}

		}

		else if (line.find("lattice_dir") != std::string::npos) {

					line_stream >> lattice_files;

					std::istringstream ss(lattice_files);
					std::string lattice_dir;

					while(std::getline(ss, lattice_dir, ',')) {
						bool success;

						for (const auto &lattice_path : std::filesystem::directory_iterator(lattice_dir)){

							std::string lattice_file = lattice_path.path().filename();
							MatrixInt lattice = loadLatticeFromFile(lattice_path.path().relative_path(), &success);

							if(success){

								std::string name = lattice_file;

								const size_t last_slash_idx = name.find_last_of("\\/");
								if (std::string::npos != last_slash_idx)
									name.erase(0, last_slash_idx + 1);

								const size_t period_idx = name.rfind('.');
								if (std::string::npos != period_idx)
									name.erase(period_idx);

								lattices.push_back(Lattice(lattice, name));
							}
							logi(lattice_file + " loaded");
						}
					}
				}

		else if (line.find("verbose") != std::string::npos) {
			std::string temp;
			line_stream >> temp;
			if(temp == "true")
				verbose = true;
		}

		line_stream.clear();
	}

	ifs.close();

}

/*
std::ifstream file(filename);
	std::vector<std::string> lines;

	std::string line;
	while (getline(file, line)){
	    if (!line.empty())
	        lines.push_back(line);
	}

	file.close();

	std::string line0 = lines[0];

	int rows = lines.size();
	int cols = 1;
	for (int i = 0; i < line0.length(); i++) {
	    if (isspace(line0.at(i))){
	      cols++;
	    }
	}

	matrix res("B", rows, cols, true);

	int row_index = 0;
	int col_index = 0;

	for(auto &line: lines){

         if(line[1] == '[')
        	 line = line.substr(2, line.length()-1);
         else if(line[line.length()-2] == ']')
        	 line = line.substr(1, line.length()-2);
         else
        	 line = line.substr(1, line.length()-1);

		 std::istringstream in(line, std::istringstream::in);

		 double n;

		 col_index = 0;
		 while (in >> n)
		    res.data[col_index++ + row_index * res.cols] = n;

		 row_index += 1;
	}

	logd(filename + " loaded");
*/
