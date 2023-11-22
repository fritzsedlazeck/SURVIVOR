/*
 * Density_VCF.cpp
 *
 *  Created on: May 15, 2020
 *      Author: sedlazec
 */

#include "Density_VCF.h"

void density_VCF(std::string vcf_file, int window, std::string output) {

	std::string buffer;
	std::ifstream myfile;

	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}

	std::vector<strvcfentry> calls;
	getline(myfile, buffer);

	int num_samples = 0;
	while (!myfile.eof()) {
		if (buffer[0] == '#' && buffer[1] == 'C') { //find out how many samples:
			int count = 0;

			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count >= 9 && buffer[i - 1] == '\t') {
					num_samples++;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
		}
		getline(myfile, buffer);
	}
	myfile.close();

	std::cerr << "Found: " << num_samples << " Samples" << std::endl;
	std::vector<strvcfentry> entries = parse_vcf(vcf_file, -1);

	std::vector<double> dups;
	std::vector<double> dels;
	std::vector<double> inv;
	std::vector<double> ins;

	for (size_t i = 0; i < entries.size(); i++) {
		int start = entries[i].start.pos / window;
		int stop = entries[i].stop.pos / window;

		while (stop + 1 > dups.size()) {
			dups.push_back(0);
			dels.push_back(0);
			inv.push_back(0);
			ins.push_back(0);
		}

		double ratio = (double) entries[i].supp / (double) num_samples;
		int pos = start;
		//	for (int pos = start; pos < stop; pos++) {
		if (entries[i].type == 0) { //DEL
			dels[pos] += ratio;
		} else if (entries[i].type == 1) { //DUP
			dups[pos] += ratio;
		} else if (entries[i].type == 2) { //INV
			inv[pos] += ratio;
		} else if (entries[i].type == 4) { //INS
			ins[pos] += ratio;
		}
		//	}
	}
	FILE * file;
	file = fopen(output.c_str(), "w");
	fprintf(file, "%s", "Start\tStop\tDEL\tDUP\tINV\tINS\n");
	for (size_t i = 0; i < dups.size(); i++) {
		//for (size_t pos = (int) i * window; pos != (int) (i + 1) * window; pos++) {
		fprintf(file, "%i", (int) i * window);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", (int) (i + 1) * window);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", dels[i]);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", dups[i]);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", inv[i]);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", ins[i]);
		fprintf(file, "%c", '\n');
		//}
	}
	fclose(file);
}

