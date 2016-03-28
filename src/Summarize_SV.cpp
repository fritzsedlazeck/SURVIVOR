/*
 * Summarize_SV.cpp
 *
 *  Created on: Nov 18, 2015
 *      Author: fsedlaze
 */

#include "Summarize_SV.h"

void adjust(std::vector<int> & vec, int dist) {
	while (vec.size() < dist+1) {
		vec.push_back(0);
	}
}
void summary_SV(std::string filename, std::string output) {
	std::vector<int> len_Del;
	std::vector<int> len_Dup;
	std::vector<int> len_Inv;
	int TRA = 0;

	std::map<std::string, std::map<std::string, int> > SV_chrs;

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			std::string chr;
			int start = 0;
			int stop = 0;
			std::string type;
			//cout<<buffer<<endl;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 0 && buffer[i] != '\t') {
					chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					start = atoi(&buffer[i]);
				}
				if (count == 4 && (buffer[i] != '\t' && (buffer[i] != '>' && buffer[i] != '<'))) {
					type += buffer[i];
				}
				if (count == 7 && strncmp("END=", &buffer[i], 4) == 0) {
					stop = atoi(&buffer[i + 4]);
					break;
				}

				if (buffer[i] == '\t') {
					count++;
				}
			}

			int step = 1000;
			int dist = abs(stop - start) / step;
			//cout<<dist<<endl;
			if (SV_chrs.find(chr) != SV_chrs.end() || SV_chrs[chr].find(type) != SV_chrs[chr].end()) {
				SV_chrs[chr][type]++;
			} else {
				SV_chrs[chr][type] = 1;
			}

			if (strcmp(type.c_str(), "DEL") == 0) {
				if (dist > len_Del.size()) {
				//	std::cout << type.c_str()<<" " << dist << endl;
				}
				adjust(len_Del, dist);
				len_Del[dist]++;
			} else if (strcmp(type.c_str(), "INV") == 0) {
				if (dist > len_Inv.size()) {
					//	std::cout << type.c_str()<<" " << dist << endl;
					//	std::cout<<buffer<<std::endl;
				}
				adjust(len_Inv, dist);
				len_Inv[dist]++;
			} else if (strcmp(type.c_str(), "DUP") == 0) {
				if (dist > len_Dup.size()) {
				//			std::cout << type.c_str()<<" " << dist << endl;
				}
				adjust(len_Dup, dist);
				len_Dup[dist]++;
			} else if (strcmp(type.c_str(), "TRA") == 0) {
				TRA++;
			}
		}
		myfile.getline(buffer, buffer_size);
	}

	std::cout<<"Parsing done: "<<TRA<<endl;
	myfile.close();

	int maxim = std::max(std::max((int) len_Del.size(), (int) len_Dup.size()), (int) len_Inv.size());
	FILE * file;
	file = fopen(output.c_str(), "w");
	fprintf(file, "%s", "Del(1kb)\tDup(1kb)\tInv(1kb)\tTRA(1kb)\n");
	for (size_t i = 0; i < maxim + 1; i++) {
		if (i < len_Del.size()) {
			fprintf(file, "%i", len_Del[i]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if (i < len_Dup.size()) {
			fprintf(file, "%i", len_Dup[i]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if (i < len_Inv.size()) {
			fprintf(file, "%i", len_Inv[i]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if (i == 0) {
			fprintf(file, "%i", TRA);
		} else {
			fprintf(file, "%i", 0);
		}

		fprintf(file, "%c", '\n');
	}
	fclose(file);

	std::string out = output;
	out += "_CHR";
	file = fopen(out.c_str(), "w");
	bool flag = true;
	for (std::map<std::string, std::map<std::string, int> >::iterator i = SV_chrs.begin(); i != SV_chrs.end(); i++) {

		if (flag) { //print the header:
			fprintf(file, "%s", "Chr");
			for (std::map<std::string, int>::iterator j = (*i).second.begin(); j != (*i).second.end(); j++) {
				fprintf(file, "%c", '\t');
				fprintf(file, "%s", (*j).first.c_str());
			}
			fprintf(file, "%c", '\n');
			flag = false;
		}

		fprintf(file, "%s", (*i).first.c_str());
		fprintf(file, "%c", '\t');
		if ((*i).second.find("DEL") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["DEL"]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if ((*i).second.find("DUP") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["DUP"]);
		} else {
			fprintf(file, "%i", 0);
		}

		fprintf(file, "%c", '\t');
		if ((*i).second.find("INV") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["INV"]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if ((*i).second.find("INS") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["INS"]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if ((*i).second.find("TRA") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["TRA"]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}
