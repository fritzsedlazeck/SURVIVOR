/*
 * Convert_MUMmer.cpp
 *
 *  Created on: Jul 31, 2017
 *      Author: sedlazec
 */

#include "Convert_MUMmer.h"

std::string print_entry_mummer(std::string chr, std::string type, int start, int stop, int len) {

	std::ostringstream convert;   // stream used for the conversion
	convert << chr;
	convert << "\t";
	convert << start;      // insert the textual representation of 'Number' in the characters in the stream
	convert << "\t";
	convert << type;
	convert << "00";
	convert << "MUMmer\tN\t<";
	convert << type;
	convert << ">\t.\tLowQual\tIMPRECISE;SVTYPE=";
	convert << type;
	convert << ";SVMETHOD=MUMmer;CHR2=";
	convert << chr;
	convert << ";END=";
	convert << stop;
	convert << ";SVLEN=";
	convert << len;
	convert << ";PE=";
	convert << 1;
	convert << "\tGT:GL:GQ:FT:RC:DR:DV:RR:RV\t";
	std::stringstream s;
	s << "1/1:0,0,0:0:PASS:0:0:";
	s << 1;
	s << ":0:0";

	return convert.str();
}

void convert_mummer_svs(std::string mummer, int min_len, std::string output) {

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(mummer.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "MUMmer Parser: could not open file: " << mummer.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size); //avoid header
	while (!myfile.eof() && buffer[0] != '[') {
		myfile.getline(buffer, buffer_size);
	}
	myfile.getline(buffer, buffer_size);
	FILE *file;
	file = fopen(output.c_str(), "w");

	while (!myfile.eof()) {
		int count = 0;
		std::string chr = "";
		std::string type = "";
		int start = 0;
		int stop = 0;
		int len = 0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				chr += buffer[i];
			}
			if (count == 1 && buffer[i] != '\t') {
				type += buffer[i];
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

		if (len > min_len) {
			fprintf(file, "%s",  print_entry_mummer(chr, type, start, stop, len).c_str());
			fprintf(file, "%c", '\n');
		}

		myfile.getline(buffer, buffer_size);
	}
	fclose(file);
	myfile.close();
}

