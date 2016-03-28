/*
 * DetectDif.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: fsedlaze
 */

#include "DetectDif.h"

int get_stop(char * buffer, int & type) {
	size_t i = 0;
	int stop=0;
	if (strncmp(&buffer[i], "TRA",3) == 0) {
		type = 3;
	} else if (strncmp(&buffer[i], "INV",3) == 0) {
		type = 2;
	} else if (strncmp(&buffer[i], "DUP",3) == 0) {
		type = 1;
	} else if (strncmp(&buffer[i], "DEL",3) == 0) {
		type = 0;
	}
	i+=3;
	while (buffer[i] != '\t') {
		if (buffer[i] == '.' && buffer[i + 1] == '.') {

			stop = atoi(&buffer[i + 2]);

		}
		i++;
	}
	return stop;
}

std::vector<svs_str> parse_strains(std::string file, int call) {
	std::ifstream myfile;

	myfile.open(file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Could not open file: " << file.c_str() << std::endl;
		exit(0);
	}

	size_t buffer_size = 2000000;
	char * buffer = new char[buffer_size];
	myfile.getline(buffer, buffer_size);
	std::vector<svs_str> strain;
	while (!myfile.eof()) {
		int count = 0;
		svs_str svs;
		svs.joined = false;
		bool called=false;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				svs.chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				svs.start = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				svs.stop = get_stop(&buffer[i], svs.type);

			}
			if (count == 3 + call && buffer[i - 1] == '\t') {
				if (buffer[i + 2] == '1') {
					called = true;
				}
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		if (called) {
			strain.push_back(svs);
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	return strain;
}

void detect_divergence(std::string file, float precent_overlap, std::string output) {
	//parse info from vcf file.
	std::vector<svs_str> p1, p2;
	p1 = parse_strains(file, 0);
	p2 = parse_strains(file, 1);



	for (size_t i = 0; i < p1.size(); i++) {
		for (size_t j = 0; j < p2.size(); j++) {
			if (p1[i].type < 2 && p1[i].type == p2[j].type) { //only del and dups:
				if (strcmp(p1[i].chr.c_str(), p2[i].chr.c_str()) == 0) {
					double dist = 0;
					if ((p1[i].start < p2[j].start && p2[j].start < p1[i].stop) || (p1[i].start < p2[j].stop && p2[j].stop < p1[i].stop)) {
						int start = max(p1[i].start, p2[j].start);
						int stop = min(p1[i].stop, p2[j].stop);
						dist = stop - start;
						double len = max((p1[i].stop - p1[i].start), (p2[j].stop - p2[j].start));
						if (dist / len > precent_overlap) {
							//mark them as joined!
							p1[i].joined = true;
							p2[i].joined = true;
						}
					}
				}
			}
		}
	}
	cout << "P1:" << endl;
	for (size_t i = 0; i < p1.size(); i++) {
		cout << "\t";
		if (p1[i].joined) {
			cout << "Joined ";
		} else {
			cout << "NOT    ";
		}
		cout << p1[i].start << " " << p1[i].stop <<" "<< p1[i].type << endl;
	}

	cout << "P2:" << endl;
	for (size_t i = 0; i < p2.size(); i++) {
		cout << "\t";
		if (p2[i].joined) {
			cout << "Joined ";
		} else {
			cout << "NOT    ";
		}
		cout << p2[i].start << " " << p2[i].stop <<" "<< p2[i].type << endl;
	}

}
