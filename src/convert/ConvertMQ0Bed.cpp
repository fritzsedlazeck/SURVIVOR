/*
 * ConvertMQ0Bed.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: fsedlaze
 */

#include "ConvertMQ0Bed.h"

void comp_mq0bed(std::string cov_file, int border, int cov_tresh) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(cov_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Cov Parser: could not open file: " << cov_file.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	int start = 0;
	int current = -1;
	int prev = current;
	std::string start_chr;
	std::string chr;
	while (!myfile.eof()) {
		int count = 0;
		int cov = 0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {

			if (count == 0 && buffer[i] != '\t') {
				chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				current = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				cov = atoi(&buffer[i]);
				break;
			}

			if (buffer[i] == '\t') {
				count++;
			}
		}
		//std::cout<<current<<" "<<prev<<std::endl;
		if (cov > cov_tresh) {
			if (prev != -1 && abs(current - prev) > border) {
				//print
				if (start - border > 1) {
					std::cout << start_chr << "\t" << start - border << "\t" << prev + border << std::endl;
				} else {
					std::cout << start_chr << "\t" << 1 << "\t" << prev + border << std::endl;
				}
				start = current;
				start_chr = chr;
			} else if (prev == -1) {
				start = current;
				start_chr = chr;
			}
			prev = current;
		}
		chr.clear();
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
}
