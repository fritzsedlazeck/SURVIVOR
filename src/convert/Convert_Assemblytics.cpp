/*
 * Convert_Assemblytics.cpp
 *
 *  Created on: May 26, 2016
 *      Author: fsedlaze
 */

#include "Convert_Assemblytics.h"

short get_type_assemblytics(std::string type) {
	if(strcmp(type.c_str(), "Tandem_expansion") == 0){
		return 1;
	}else if (strcmp(type.c_str(), "Deletion") == 0 || (strcmp(type.c_str(), "Repeat_contraction") == 0 || strcmp(type.c_str(), "Tandem_contraction") == 0)) {
		return 0;
	} else if (strcmp(type.c_str(), "Insertion") == 0 || strcmp(type.c_str(), "Repeat_expansion") == 0 ) {
		return 4;
	} else {
		std::cerr << "Unknown type! " << type << std::endl;
	}
	return -1;
}
void parse_assemblytics(std::string assemblytics,int minlen, std::vector<strvcfentry> & entries) {
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
	while (!myfile.eof()) {
		//	std::cout<<buffer<<std::endl;
		int count = 0;
		strvcfentry tmp;
		std::string type;
		int dist=0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				tmp.start.chr += buffer[i];
				tmp.stop.chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				tmp.start.pos = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				tmp.stop.pos = atoi(&buffer[i]);
			}
			if (count == 6 && buffer[i] != '\t') {
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
		tmp.type = get_type_assemblytics(type);
		if(tmp.type!=0 && tmp.type!=1){
			tmp.stop.pos+=dist;
		}
		if(tmp.stop.pos-tmp.start.pos > minlen){
			entries.push_back(tmp);
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
}

std::string print_entry(strvcfentry & region) {

//	III     5104    DEL00000002     N       <DEL>   .       LowQual IMPRECISE;CIEND=-305,305;CIPOS=-305,305;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.5.9;CHR2=III;END=15991;SVLEN=10887;CT=3to5;PE=2;MAPQ=60        GT:GL:GQ:FT:RC:DR:DV:RR:RV      1/1:-12,-0.602059,0:6:LowQual:816:0:2:0:0

	std::ostringstream convert;   // stream used for the conversion
	convert << region.start.chr;
	convert << "\t";
	convert << region.start.pos;      // insert the textual representation of 'Number' in the characters in the stream
	convert << "\t";
	convert << trans_type(region.type);
	convert << "00";
	convert << "Assym\tN\t<";
	convert << trans_type(region.type);
	convert << ">\t.\tLowQual\tIMPRECISE;SVTYPE=";
	convert << trans_type(region.type);
	convert << ";SVMETHOD=PINDELv0.2.5a8;CHR2=";
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

void process_Assemblytics(std::string assemblytics,int minlen, std::string output) {

	std::vector<strvcfentry> entries;
	parse_assemblytics(assemblytics,minlen, entries);
	FILE *file;
	file = fopen(output.c_str(), "w");
	for (size_t i = 0; i < entries.size(); i++) {
		fprintf(file, "%s", print_entry(entries[i]).c_str());
		fprintf(file, "%c", '\n');
	}

	fclose(file);

}
