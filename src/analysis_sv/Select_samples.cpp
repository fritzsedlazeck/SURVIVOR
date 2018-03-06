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
	} else if (buffer[0] == '0' && buffer[2] == '0') {
		return false;
	}
	if (strncmp(buffer, "./.:0:0,0:--:NaN:NaN", 20) != 0) {
		return false;
	}
	//0/0 ./.
	return false;
}

std::vector<int>  parase_matrix(std::string vcf_file, std::vector<std::string> & names, std::map<int, bool> taken_ids, int &num) {
	std::cout<<"Parsing... "<<endl;
	std::string buffer;
	std::ifstream myfile;
	std::vector<int>  matrix;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}

	getline(myfile, buffer);
	int line=0;
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
			if (matrix.empty()) { //init pairwise matrix;
				std::vector<int> tmp;
				matrix.resize(names.size(), 0);
			}
			line++;
			num++;
			//bool discard = false;
			int count = 0;
			bool include = false;
			int num=0;

			std::string entries;
			entries.resize(names.size(),'0');
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count >= 9 && buffer[i - 1] == '\t') {
					if (genotype_parse(&buffer[i])) {
						if (taken_ids.find(num) != taken_ids.end()) {
							include=false;
							break;
						}
						include=true;
						entries[num]='1';
					}
					num++;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			if(include){
				for(size_t j=0;j<entries.size();j++){
					if(entries[j]=='1'){
						matrix[j]++;
					}
				}
			}
			if(line%10000==0){
				std::cout <<"lines: "<<line<<std::endl;
			}
		}

		getline(myfile, buffer);
	}

	cout<<"Fin"<<endl;
	myfile.close();
	return matrix;
}
void print_mat(std::vector<int>  svs_count_mat) {
	for (size_t i = 0; i < svs_count_mat.size(); i++) {
			std::cout << svs_count_mat[i] << "\t";
	}
	std::cout << std::endl;
	std::cout << std::endl;

}
void select_greedy(std::string vcf_file, std::string output) {
	std::map<int, bool> taken_ids;
	std::vector<std::string> sample_names;
	int total_svs = 0;

	//we can actually just use a vector instead!
	std::vector<int> svs_count_mat = parase_matrix(vcf_file, sample_names, taken_ids, total_svs); //span a  NxN matrix and stores the shared SVs
	//print_mat(svs_count_mat);

	FILE *file;
	file = fopen(output.c_str(), "w");

	fprintf(file, "%s", "Sample\t#SVs\t#_SVs_captured\t%_SVs_captured\n");
	std::cout << "Parsed vcf file with " << sample_names.size() << " samples" << endl;
	int captured_svs = 0;
	for (size_t i = 0; i < sample_names.size(); i++) {
		//select max on main diag
		int max = 0;
		int max_id = -1;

		for (size_t j = 0; j < sample_names.size(); j++) {
			//	cout << svs_count_mat[j][j] << "\t";
			if (max < svs_count_mat[j]) {
				max = svs_count_mat[j];
				max_id = j;
			}
		}
		captured_svs += max;
		std::cout <<"RANK:\t"<<i<<"\t"<<sample_names[max_id]<<"\t"<<max<<"\t"<<captured_svs<<"\t"<< (double) captured_svs / (double) total_svs<< std::endl;
		fprintf(file, "%s", sample_names[max_id].c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", max);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", captured_svs);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", (double) captured_svs / (double) total_svs);
		fprintf(file, "%c", '\n');

		if (max == 0) {
			break;
		}
	//	print_mat(svs_count_mat);
		//print max_id and max
		//erase joined svs over matrix given the pairwise matrix.
		taken_ids[max_id] = true;
		total_svs = 0;
		svs_count_mat = parase_matrix(vcf_file, sample_names, taken_ids, total_svs); //span a  NxN matrix and stores the shared SVs
	}
	fclose(file);
}

