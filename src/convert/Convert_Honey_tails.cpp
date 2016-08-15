/*
 * Convert_Honey_tails.cpp
 *
 *  Created on: Jun 6, 2016
 *      Author: fsedlaze
 */
#include "Convert_Honey_tails.h"

short get_type_honey(std::string type){

	if (strncmp(type.c_str(), "DEL", 3) == 0 ) {
		return 0;
	} else if (strncmp(type.c_str(), "DUP", 3) == 0) {
		return 1;
	} else if (strncmp(type.c_str(), "INV", 3) == 0) {
		return 2;
	} else if (strncmp(type.c_str(), "TLOC", 4) == 0) {
		return 3;
	} else if (strncmp(type.c_str(), "INS", 3) == 0) {
		return 4;
	} else {
		std::cerr << "Unknown type! "<<type << std::endl;
	}
	return -1;

}

void parse_honey_tails(std::string assemblytics,int minlen, std::vector<strvcfentry> & entries) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(assemblytics.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Pindel Parser: could not open file: " << assemblytics.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size); //avoid header
	myfile.getline(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		//	std::cout<<buffer<<std::endl;
		int count = 0;
		strvcfentry tmp;
		std::string type;
		int dist=0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 2 && buffer[i] != '\t') {
				tmp.start.chr += buffer[i];
			}

			if (count == 5 && buffer[i] != '\t') {
				tmp.stop.chr += buffer[i];
			}
			if (count == 3 && buffer[i - 1] == '\t') {
				tmp.start.pos = atoi(&buffer[i]);
			}
			if (count == 6 && buffer[i - 1] == '\t') {
				tmp.stop.pos = atoi(&buffer[i]);
			}
			if (count == 9 && buffer[i] != '\t') {
				type += buffer[i];
			}
			//if (count == 7 && buffer[i - 1] == '\t') {
			//	dist= abs(atoi(&buffer[i]));
			//}
			if (count == 4 && buffer[i - 1] == '\t') {
				dist= abs(atoi(&buffer[i]));
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		tmp.type = get_type_honey(type);
		if(tmp.stop.pos-tmp.start.pos > minlen){
			entries.push_back(tmp);
		}
		myfile.getline(buffer, buffer_size);
	}
}

std::string print_entry_honey(strvcfentry & region) {

//	III     5104    DEL00000002     N       <DEL>   .       LowQual IMPRECISE;CIEND=-305,305;CIPOS=-305,305;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.5.9;CHR2=III;END=15991;SVLEN=10887;CT=3to5;PE=2;MAPQ=60        GT:GL:GQ:FT:RC:DR:DV:RR:RV      1/1:-12,-0.602059,0:6:LowQual:816:0:2:0:0

	std::ostringstream convert;   // stream used for the conversion
	convert << region.start.chr;
	convert << "\t";
	convert << region.start.pos;      // insert the textual representation of 'Number' in the characters in the stream
	convert << "\t";
	convert << trans_type(region.type);
	convert << "00";
	convert << "Honey\tN\t<";
	convert << trans_type(region.type);
	convert << ">\t.\tLowQual\tIMPRECISE;SVTYPE=";
	convert << trans_type(region.type);
	convert << ";SVMETHOD=Honey_tails;CHR2=";
	convert << region.stop.chr;
	convert << ";END=";
	convert << region.stop.pos;
	convert << ";SVLEN=";
	convert << region.stop.pos - region.start.pos;
	convert << ";PE=";
	convert << 1;
	convert << "\tGT:GL:GQ:FT:RC:DR:DV:RR:RV\t";
	std::stringstream s;
	s << "1/1:0,0,0:0:PASS:0:0:";
	s << 1;
	s << ":0:0";
	//std::cout<<convert.str()<<std::endl;
	return convert.str();
}

void process_Honey(std::string assemblytics, int minlen, std::string output) {
	std::vector<strvcfentry> entries;
	parse_honey_tails(assemblytics, minlen, entries);
	FILE *file;
	file = fopen(output.c_str(), "w");
	for (size_t i = 0; i < entries.size(); i++) {
		fprintf(file, "%s", print_entry_honey(entries[i]).c_str());
		fprintf(file, "%c", '\n');
	}

	fclose(file);
}
