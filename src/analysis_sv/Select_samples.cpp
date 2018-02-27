/*
 * Select_samples.cpp
 *
 *  Created on: Feb 27, 2018
 *      Author: sedlazec
 */

#include "Select_samples.h"

bool genotype_parse(char * buffer) {
	//cout << "buffer: " << buffer[0] << buffer[1] << buffer[2] << endl;

	if ((buffer[0] == '0' && buffer[2] == '1') || (buffer[0] == '1' && buffer[2] == '1')) {
		return true;
	}
	if (strncmp(buffer, "./.:0:0,0:--:NaN:NaN", 20) != 0) {
		return true;
	}
	//0/0 ./.
	return false;
}

std::vector<std::vector<int> > parase_matrix(std::string vcf_file, std::vector<std::string> & names, std::map<int, bool> taken_ids) {
	std::string buffer;
	std::ifstream myfile;
	std::vector<std::vector<int> > matrix;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}

	getline(myfile, buffer);
	while (!myfile.eof()) {
		if (names.empty() && (buffer[0] == '#' && buffer[1] == 'C')) { //parse names
			int count = 0;
			std::string id = "";
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count >= 9 && buffer[i] != '\t') {
					id += buffer[i];
				}
				if (buffer[i] == '\t') {
					if (!id.empty()) {
						names.push_back(id);
						id = "";
					}
					count++;
				}
			}
			if (!id.empty()) {
				names.push_back(id);
			}

		} else if (buffer[0] != '#') { //parse svs;

			bool discard = false;
			int count = 0;
			std::string entry;
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count >= 9 && buffer[i - 1] == '\t') {
					if (genotype_parse(&buffer[i])) {
						if(taken_ids.find((int)entry.size())!=taken_ids.end()){
							discard=true;
							break;
						}
						entry += '1';

					} else {
						entry += '0';
					}
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			if (matrix.empty()) { //init pairwise matrix;
				std::vector<int> tmp;
				tmp.resize(names.size(), 0);
				matrix.resize(names.size(), tmp);
			}

			if (!discard) {
				for (size_t i = 0; i < entry.size(); i++) {
					for (size_t j = 0; j < entry.size(); j++) {
						if (entry[i] == '1' && entry[j] == '1') {
							matrix[i][j]++;
						}
					}
				}
			}
		}

		getline(myfile, buffer);
	}
	myfile.close();
	return matrix;
}
void print_mat(std::vector<std::vector<int> > svs_count_mat) {
	for (size_t i = 0; i < svs_count_mat.size(); i++) {
		for (size_t j = 0; j < svs_count_mat[i].size(); j++) {
			std::cout << svs_count_mat[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;

}
void select_greedy(std::string vcf_file, std::string output) {
	std::map<int, bool> taken_ids;
	std::vector<std::string> sample_names;

	std::vector<std::vector<int> > svs_count_mat = parase_matrix(vcf_file, sample_names, taken_ids); //span a  NxN matrix and stores the shared SVs
	//print_mat(svs_count_mat);

	int total_svs = 0;
	for (size_t i = 0; i < sample_names.size(); i++) {
		//select max on main diag
		int max = 0;
		int max_id = -1;

		for (size_t j = 0; j < sample_names.size(); j++) {
		//	cout << svs_count_mat[j][j] << "\t";
			if (max < svs_count_mat[j][j]) {
				max = svs_count_mat[j][j];
				max_id = j;
			}
		}
		std::cout << std::endl;

		total_svs += max;
		//print max_id and max
		//erase joined svs over matrix given the pairwise matrix.
		taken_ids[max_id] = true;
		svs_count_mat = parase_matrix(vcf_file, sample_names, taken_ids); //span a  NxN matrix and stores the shared SVs
	}
}

