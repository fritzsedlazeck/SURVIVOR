/*
 * CorrectAllele.cpp
 *
 *  Created on: Jun 18, 2015
 *      Author: fsedlaze
 */

#include "CorrectAllele.h"

bool eval_sv(bool flag, double p22, double p32, double ratiop22,
		double ratiop23) {
	if (flag) {
		//DEL: (P22 and P32 < 1e-10) and (alt.JB22.ratio and alt.JB32.ratio < 0.2)
		return ((p22 < 1e-10 && p32 < 1e-10)
				&& (ratiop22 < 0.2 && ratiop23 < 0.2));
	}
	//DUP: (P22 and P32 < 1e-10) and (alt.JB22.ratio and alt.JB32.ratio > 1.8)
	return ((p22 < 1e-10 && p32 < 1e-10) && (ratiop22 < 1.8 && ratiop23 < 1.8));
}

std::map<std::string, std::map<std::string, bool> > get_entries(
		std::string table) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(table.c_str(), std::ifstream::in);
	myfile.getline(buffer, buffer_size);
	std::map<std::string, std::map<std::string, bool> > entries;
	std::string prevname = "";
	bool flag = (bool) (buffer[2] == 'E');
	while (!myfile.eof()) {
		int count = 0;
		std::string name = "";
		double p22 = 100;
		double p32 = 100;
		double ratiop22 = 0;
		double ratiop32 = 0;
		std::string id;
		for (size_t i = 0; buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				id += buffer[i];
			}
			if (count == 2 && buffer[i] != '\t') {
				name += buffer[i];
			}
			if (count == 4 && buffer[i - 1] == '\t') {
				p22 = atof(&buffer[i]);
			}
			if (count == 5 && buffer[i - 1] == '\t') {
				p32 = atof(&buffer[i]);
			}
			if (count == 6 && buffer[i - 1] == '\t') {
				ratiop22 = atof(&buffer[i]);
			}
			if (count == 7 && buffer[i - 1] == '\t') {
				ratiop32 = atof(&buffer[i]);
				bool test=eval_sv(flag, p22, p32, ratiop22, ratiop32);
				if(test){
					std::cout<<name<<std::endl;
				}else{
					//std::cout<<"No HIT"<<std::endl;
				}
				entries[id][name] = eval_sv(flag, p22, p32, ratiop22, ratiop32);

			}
			if (buffer[i] == '\t') {
				count++;
			}
		}

		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	return entries;
}

void correct_alleles(std::string vcf_file, std::string table,
		std::string output) {

	std::map<std::string, std::map<std::string, bool> > entries = get_entries(
			table);
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	myfile.getline(buffer, buffer_size);
	FILE *file;
	file = fopen(output.c_str(), "w");
	while (!myfile.eof() && buffer[0] == '#' && buffer[1] == '#') {
		fprintf(file, "%s", buffer);
		fprintf(file, "%c", '\n');
		myfile.getline(buffer, buffer_size);
	}
	fprintf(file, "%s", buffer);
	fprintf(file, "%c", '\n');
	std::vector<std::string> names;
	std::string name;
	int count = 0;
	for (size_t i = 0;
			i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
		if (count > 8 && buffer[i] != '\t') {
			name += buffer[i];
		}
		if (buffer[i] == '\t') {
			count++;
			if (!name.empty()) {
				std::cout<<name<<std::endl;
				names.push_back(name);
				name.clear();
			}
		}
	}

	myfile.getline(buffer, buffer_size);

	//parse vcf:
	while (!myfile.eof()) {
		std::string chr;
		std::string header;
		std::string pos;
		std::string type;
		int count = 0;
		for (size_t i = 0;
				i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n';
				i++) {
			if (count == 0 && buffer[i] != '\t') {
				chr += buffer[i];
			}
			if (count == 1 && buffer[i] != '\t') {
				pos += buffer[i];
			}
			if ((count == 4 && buffer[i] != '<')
					&& (buffer[i] != '>' && buffer[i] != '\t')) {
				type += buffer[i];
			}
			if (count == 8 && buffer[i - 1] == '\t') {
				std::string name = type;
				name += '.';
				name += chr;
				name += '.';
				name += pos;
			//	std::cout<<name<<std::endl;
				if (entries.find(name) == entries.end()) {
					//print full;
					fprintf(file, "%s", buffer);
					fprintf(file, "%c", '\n');
					break;
				} else {
					int count2 = 0;
					int id=0;
					size_t j=0;
					while((j) < buffer_size&& buffer[(j)] != '\0' && buffer[(j)] != '\n'){
					//for (size_t j = 0; j < buffer_size; j++) { ///print header up to entries:
						if (count2 < 9) {
							fprintf(file, "%c", buffer[(j)]);
						} else if (buffer[(j) - 1] == '\t') {
						 	if (entries[name][names[id]]) {
								fprintf(file, "%s", "1/1");
							} else {
								fprintf(file, "%s", "0/0");
							}
							j+=2;
							id++;
						}else{
							fprintf(file, "%c", buffer[(j)]);
						}
						if (buffer[(j)] == '\t') {
							count2++;
						}
						j++;
					}
					fprintf(file, "%c", '\n');
				}
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}

		myfile.getline(buffer, buffer_size);
	}
	myfile.close();

}
