/*
 * MUMmer_overlap.cpp
 *
 *  Created on: Dec 27, 2017
 *      Author: sedlazec
 */

#include "MUMmer_overlap.h"

vector<string> parse_filenames(std::string filename) {
	std::vector<std::string> names;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "File Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		names.push_back(std::string(buffer));
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();

	return names;
}

void comp_entries_mummer(std::vector<strvcfentry> & entries, std::string filename, int max_dist) {

	std::size_t found = filename.find("short_");
	std::string id_denovo = "";
	if (found != std::string::npos) {
		id_denovo = filename.substr(found + 6);
	} else {
		cout << "not found" << endl;
		id_denovo = filename;
	}

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "File Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof() && buffer[0] != '[') { //avoid headers!
		myfile.getline(buffer, buffer_size);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		int count = 0;
		std::string chr = "";
		int start = 0;
		int stop = 0;
		int len = 0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				chr += buffer[i];
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				start = atoi(&buffer[i]);
			}
			if (count == 3 && buffer[i - 1] == '\t') {
				stop = atoi(&buffer[i]);

			}
			if (count == 4 && buffer[i - 1] == '\t') {
				len = atoi(&buffer[i]);
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		for (size_t i = 0; i < entries.size(); i++) {
			if (entries[i].num_reads.second == 0) {
				if (strcmp(entries[i].start.chr.c_str(), chr.c_str()) == 0) {
					if (abs(entries[i].start.pos - start) < max_dist) {
						//match!

						//cout<<"HIT1 "<<entries[i].start.chr <<" "<<entries[i].start.pos <<" "<<entries[i].stop.pos <<" "<<start<<" "<<stop <<endl;
						entries[i].num_reads.first++;
						entries[i].num_reads.second++;
						entries[i].calls[id_denovo] = 1;
					}
				} else if (strcmp(entries[i].stop.chr.c_str(), chr.c_str()) == 0) {
					if (abs(entries[i].stop.pos - start) < max_dist) {
						//match!

						//	cout<<"HIT2 "<<entries[i].start.chr <<" "<<entries[i].start.pos <<" "<<entries[i].stop.pos <<" "<<start<<" "<<stop <<endl;
						entries[i].num_reads.first++;
						entries[i].num_reads.second++;
						entries[i].calls[id_denovo] = 1;
					}
				}
			}
		}

		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	for (size_t i = 0; i < entries.size(); i++) {
		entries[i].num_reads.second = 0; //reset flag!
	}
}

void overlapp_mummer(std::string vcf_SVs_file, std::string mummer_files, int max_dist, std::string output) {
//parse VCF file:
	std::vector<strvcfentry> entries = parse_vcf(vcf_SVs_file, 0);

	for (size_t i = 0; i < entries.size(); i++) {
		//init to use it later!
		entries[i].num_reads.first = 0; //total counts
		entries[i].num_reads.second = 0; //flag to no count twice!
	}
//compare to MUMMer files:
	vector<string> filenames = parse_filenames(mummer_files);
	for (size_t i = 0; i < filenames.size(); i++) {
		comp_entries_mummer(entries, filenames[i], max_dist);
	}

	//combine info:
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(vcf_SVs_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "File Parser: could not open file: " << vcf_SVs_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);

	FILE *file;
	file = fopen(output.c_str(), "w");
	int line = 0;
	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				fprintf(file, "%c", buffer[i]);
				if (count == 7 && buffer[i + 1] == '\t') {
					if (entries[line].calls.size() > 1) {
						fprintf(file, "%s", ";DeNovo_mummer=");
						bool draw_comma=false;
						for (std::map<std::string, std::string>::iterator j = entries[line].calls.begin(); j != entries[line].calls.end(); j++) {
							if (strcmp((*j).first.c_str(), vcf_SVs_file.c_str()) != 0) {
								if (draw_comma) {
									fprintf(file, "%c", ',');
								}
								fprintf(file, "%s", (*j).first.c_str());
								draw_comma=true;
							}
						}
					}

				}
				if (buffer[i] == '\t') {
					count++;
				}

			}
			fprintf(file, "%c", '\n');
			line++;
		} else {
			fprintf(file, "%s", buffer);
			fprintf(file, "%c", '\n');
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	fclose(file);
}

