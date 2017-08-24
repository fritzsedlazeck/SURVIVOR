/*
 * Summ_mat.cpp
 *
 *  Created on: Jul 5, 2017
 *      Author: sedlazec
 */

#include "Summ_mat.h"

void process_patterns(std::vector<std::string> mat, FILE *& file) {

	std::map<std::string, int> patterns;
	std::vector<int> vec;
	for (size_t j = 0; j < mat[0].size(); j++) {
		std::string patt;
		for (size_t i = 0; i < mat.size(); i++) {
			patt += mat[i][j];
		}
		if (patterns.find(patt) == patterns.end()) {
			patterns[patt] = 1;
		} else {
			patterns[patt]++;
		}
	}

	for (std::map<std::string, int>::iterator i = patterns.begin(); i != patterns.end(); i++) {
		while ((*i).second >= (int) vec.size()) {
			vec.push_back(0);
		}
		vec[(*i).second]++;
	}
	std::stringstream ss;
	for (size_t i = 1; i < vec.size(); i++) {
		fprintf(file, "%i", (int) vec[i]);
		fprintf(file, "%c", ';');
	}
	fprintf(file, "%c", '\n');
}
char parse_inf(char * buffer) {

	size_t i = 0;
	while (buffer[i] != '\t' && buffer[i] != '\n') {
		if (strncmp("NaN", &buffer[i], 3) == 0) {
			return '0';
		}
		i++;
	}
	return '1';
}
void summarize_svs_table_window(std::string venn_file, int window, std::string output) {

	size_t buffer_size = 200000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(venn_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: " << venn_file.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);

	FILE *file;
	file = fopen(output.c_str(), "w");

	int last_pos = 0;

	std::vector<std::string> mat;
	std::string last_chr = "";
	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int pos = 0;
			std::string chr;
			int count = 0;
			std::string pattern;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 0 && buffer[i] != '\t') {
					chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					pos = atoi(&buffer[i]);

					if (pos - last_pos > window || (strcmp(chr.c_str(), last_chr.c_str()) != 0)) {
						//process entries;
						if (mat.size() > 0) {
							fprintf(file, "%s", last_chr.c_str());
							fprintf(file, "%c", ':');
							fprintf(file, "%i", (int) last_pos);
							fprintf(file, "%c", ':');
							//	fprintf(file, "%s", (*patterns.begin()).first.c_str());
							//	fprintf(file, "%c", ':');
							process_patterns(mat, file);
							mat.clear();
						}
						last_pos = pos;
						last_chr = chr;

					}
				}
				if (count > 9 && buffer[i - 1] == '\t') {
					pattern += parse_inf(&buffer[i]);

				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			mat.push_back(pattern);
		}
		myfile.getline(buffer, buffer_size);
	}
	fclose(file);
}

