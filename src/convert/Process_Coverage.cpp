/*
 * Process_Coverage.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: fsedlaze
 */

#include "Process_Coverage.h"

void summarize_badcoverage(std::string filename, int win_size, int min_cov, std::string output) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Lumpy Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}

	FILE *file;
	file = fopen(output.c_str(), "w");

	myfile.getline(buffer, buffer_size);
	myfile.getline(buffer, buffer_size); //avoiding header
	int start = win_size * -1;
	int stop = 0;
	std::string chr = "";
	int pos = 0;
	std::string chr_prev = "";
	while (!myfile.eof()) {
		chr.clear();
		int count = 0;
		int cov = -1;

		//REF     POS     COV
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				if (!chr_prev.empty() && strcmp(chr_prev.c_str(), chr.c_str()) != 0) {

					if (start != win_size * -1) {
						fprintf(file, "%s", chr_prev.c_str());
						fprintf(file, "%c", '\t');
						fprintf(file, "%i", start);
						fprintf(file, "%c", '\t');
						fprintf(file, "%i", pos);
						fprintf(file, "%c", '\n');
						start = win_size * -1;
						stop = 1;
					}
					chr_prev = chr;
				}
				pos = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				cov = atoi(&buffer[i]);
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}

		if (cov <= min_cov) {
			if (start == win_size * -1) {
				if (pos - win_size > 0) {
					start = pos - win_size;
				} else {
					start = 0;
				}
				stop = pos;
			}
			if (stop - pos < win_size) { // extend current window.
				stop = pos;
			}
		} else if (start != win_size * -1 && pos - stop > win_size) { //report:
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", start);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", stop);
			fprintf(file, "%c", '\n');
			start = win_size * -1;
			stop = pos;
		}
		chr_prev = chr;
		myfile.getline(buffer, buffer_size);
	}
	if (start != win_size * -1) {
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", start);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", pos);
		fprintf(file, "%c", '\n');
	}
	fclose(file);
	myfile.close();

}

