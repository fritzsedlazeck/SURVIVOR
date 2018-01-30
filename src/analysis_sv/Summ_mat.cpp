/*
 * Summ_mat.cpp
 *
 *  Created on: Jul 5, 2017
 *      Author: sedlazec
 */

#include "Summ_mat.h"

void process_patterns(std::vector<std::string> mat, FILE *& file, FILE *& file2) {

	std::map<std::string, int> patterns;
	std::map<std::string, std::vector<int> > patterns2;
	std::vector<int> vec;
	for (size_t j = 0; j < mat[0].size(); j++) {
		std::string patt;
		bool useful = false;
		for (size_t i = 0; i < mat.size(); i++) {
			patt += mat[i][j];
			if (mat[i][j] == '1') {
				useful = true;
			}
		}
		if (useful) {

			patterns2[patt].push_back(j);
		}
		if (patterns.find(patt) == patterns.end()) {
			patterns[patt] = 1;
		} else {
			patterns[patt]++;
		}

	}

	for (std::map<std::string, std::vector<int> >::iterator i = patterns2.begin(); i != patterns2.end(); i++) {
		fprintf(file2, "%c", '\t');
		fprintf(file2, "%i", (*i).second[0]);
		for (size_t j = 1; j < (*i).second.size(); j++) {
			fprintf(file2, "%c", ',');
			fprintf(file2, "%i", (*i).second[j]);
		}
	}
	fprintf(file2, "%c", '\n');

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

	FILE* file2;
	std::string out = output;
	out += "patient_hist";
	file2 = fopen(out.c_str(), "w");

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
							fprintf(file2, "%s", last_chr.c_str());
							fprintf(file2, "%c", ':');
							fprintf(file2, "%i", (int) last_pos);
							process_patterns(mat, file, file2);

							vector<short> patients;
							patients.assign(mat[0].size(), 0);
							for (size_t i = 0; i < mat.size(); i++) {
								for (size_t j = 0; j < mat[i].size(); j++) {
									if (mat[i][j] == '1') {
										patients[j]++;
									}
								}
							}

							for (size_t i = 0; i < patients.size(); i++) {
								fprintf(file2, "%c", '\t');
								fprintf(file2, "%i", (int) patients[i]);

							}
							fprintf(file2, "%c", '\n');
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

void summarize_svs_table_window_stream(int window, std::string output) {

	FILE *file;
	file = fopen(output.c_str(), "w");

	FILE* file2;
	std::string out = output;
	out += "perpatient";
	file2 = fopen(out.c_str(), "w");

	int last_pos = 0;

	std::vector<std::string> mat;
	std::string last_chr = "";
	while (!cin.eof()) {
		std::string buffer;
		getline(cin, buffer);
		if (!cin.fail()) {
			if (buffer[0] != '#') {
				int pos = 0;
				std::string chr;
				int count = 0;
				std::string pattern;
				for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
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

								fprintf(file2, "%s", last_chr.c_str());
								fprintf(file2, "%c", ':');
								fprintf(file2, "%i", (int) last_pos);

								process_patterns(mat, file, file2);

								vector<short> patients;
								patients.assign(mat[0].size(), 0);
								for (size_t i = 0; i < mat.size(); i++) {
									for (size_t j = 0; j < mat[i].size(); j++) {
										if (mat[i][j] == '1') {
											patients[j]++;
										}
									}
								}
								mat.clear();

							}
							last_pos = pos;
							last_chr = chr;

						}
					}
					if (count == 7 && strncmp(&buffer[i], "SUPP_VEC=", 9) == 0) {
						std::string tmp = buffer.substr(i + 9);
						std::size_t found = tmp.find_first_of(";");
						pattern = tmp.substr(0, found);
					}
					if (buffer[i] == '\t') {
						count++;
					}
				}
				mat.push_back(pattern);
			}
		}
	}
	fclose(file);
}

