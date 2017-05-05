/*
 * Convert_Bionano.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: fsedlaze
 */
//83      228     1       1       835435.0        864300.7        230129707.0     230159791.0     -1.00   deletion        181     182     -1      101     107     23768   23776
#include "Convert_Bionano.h"

short trans_type_bio(std::string type) {
	if (strncmp(type.c_str(), "deletion", 8) == 0) {
		return 0;
		//} else if (strncmp(type.c_str(), "DUP", 3) == 0) {
		//	return 1;
	} else if (strncmp(type.c_str(), "inversion", 9) == 0) { //TODO: inversion paired??
		return 2;
	} else if (strcmp(type.c_str(), "translocation_interchr") == 0 || strcmp(type.c_str(), "translocation_intrachr") == 0) { //TODO: intra != inter!!
		return 3;
	} else if (strncmp(type.c_str(), "insertion", 9) == 0) {
		return 4;
	} else {
		std::cerr << "Unknown type! " << type << std::endl;
	}
	return -1;
}

std::string print_entry_bio(strvcfentry & region) {

//	III     5104    DEL00000002     N       <DEL>   .       LowQual IMPRECISE;CIEND=-305,305;CIPOS=-305,305;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.5.9;CHR2=III;END=15991;SVLEN=10887;CT=3to5;PE=2;MAPQ=60        GT:GL:GQ:FT:RC:DR:DV:RR:RV      1/1:-12,-0.602059,0:6:LowQual:816:0:2:0:0

	std::ostringstream convert;   // stream used for the conversion
	convert << region.start.chr;
	convert << "\t";
	convert << region.start.pos;      // insert the textual representation of 'Number' in the characters in the stream
	convert << "\t";
	convert << trans_type(region.type);
	convert << "00";
	convert << "Bionanom\tN\t<";
	convert << trans_type(region.type);
	convert << ">\t.\tLowQual\tIMPRECISE;SVTYPE=";
	convert << trans_type(region.type);
	convert << ";CHR2=";
	convert << region.stop.chr;
	convert << ";END=";
	convert << region.stop.pos;
	convert << ";SVLEN=";
	convert << region.sv_len;
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

void parse_Bionano(std::string bionano, std::vector<strvcfentry>& entries) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(bionano.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Bionano Parser: could not open file: " << bionano.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size); //avoid header
	myfile.getline(buffer, buffer_size);
	while (buffer[0] == '#' && !myfile.eof()) {
		myfile.getline(buffer, buffer_size);
	}

	while (!myfile.eof()) {
		int count = 0;
		strvcfentry tmp;
		std::string type;
		int query_start = 0;
		int query_stop = 0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 2 && buffer[i] != '\t') {
				tmp.start.chr += buffer[i];
			}
			if (count == 3 && buffer[i] != '\t') {
				tmp.stop.chr += buffer[i];
			}
			if (count == 4 && buffer[i - 1] == '\t') {
				query_start = atoi(&buffer[i]);
			}
			if (count == 5 && buffer[i - 1] == '\t') {
				query_stop = atoi(&buffer[i]);
			}
			if (count == 6 && buffer[i - 1] == '\t') {
				tmp.start.pos = atoi(&buffer[i]);
			}
			if (count == 7 && buffer[i - 1] == '\t') {
				tmp.stop.pos = atoi(&buffer[i]);
			}
			if (count == 9 && buffer[i] != '\t') {
				type += buffer[i];
			}
			if (count == 10 && buffer[i - 1] == '\t') {
				tmp.type = trans_type_bio(type);
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		double factor = (query_stop - query_start) / 2;
		if (tmp.type == 0) {
			tmp.start.pos = tmp.start.pos + factor;
			tmp.stop.pos = tmp.stop.pos - factor;
			tmp.sv_len=tmp.stop.pos-tmp.start.pos;
		} else if (tmp.type == 4) {
			tmp.sv_len=query_start - query_stop + tmp.stop.pos - tmp.start.pos;

			//tmp.sv_len=abs((query_start - query_stop) - (tmp.start.pos - tmp.stop.pos));
			tmp.start.pos=(query_stop +  query_start)/2;
			tmp.stop.pos = tmp.start.pos+1;
		}
		entries.push_back(tmp);
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();

}

void parse_GC(std::string bionano, std::vector<strvcfentry>& entries) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(bionano.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Bionano Parser: could not open file: " << bionano.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size); //avoid header
	myfile.getline(buffer, buffer_size);
	while (buffer[0] == '#' && !myfile.eof()) {
		myfile.getline(buffer, buffer_size);
	}

	while (!myfile.eof()) {
		if (buffer[0] != '#' && buffer[0] != '>') {
			int count = 0;
			strvcfentry tmp;
			std::string type = "";
			tmp.type = -1;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 1 && buffer[i] != '\t') {
					type += buffer[i];
				}
				if (count == 2 && buffer[i - 1] == '\t') {
					//	std::cout<<type<<std::endl;
					if (strcmp(type.c_str(), "inversion") == 0) {
						//	std::cout<<"HIT"<<std::endl;
						tmp.type = 2;
					} else if (strcmp(type.c_str(), "deletion") == 0) {
						std::cout << "HIT" << std::endl;
						tmp.type = 0;
					} else {
						type.clear(); //flag for ignoring!
						break;
					}
				}
				if (count == 5 && buffer[i] != '\t') { //all events are on the same chrs!
					tmp.start.chr += buffer[i];
					tmp.stop.chr += buffer[i];
				}
				if (count == 6 && buffer[i - 1] == '\t') {
					tmp.start.pos = atoi(&buffer[i]);
				}
				if (count == 7 && buffer[i - 1] == '\t') {
					tmp.stop.pos = atoi(&buffer[i]);
				}
				if (count == 8 && buffer[i - 1] == '\t') {
					tmp.sv_len = atoi(&buffer[i]);
				}

				if (buffer[i] == '\t') {
					count++;
				}
			}
			if (tmp.type != -1) {
				entries.push_back(tmp);
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();

}

void process_Bionano(std::string bionano, std::string output) {
	std::vector<strvcfentry> entries;
	parse_Bionano(bionano, entries);
	FILE *file;
	file = fopen(output.c_str(), "w");
	for (size_t i = 0; i < entries.size(); i++) {
		fprintf(file, "%s", print_entry_bio(entries[i]).c_str());
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}

void process_CG(std::string gc_file, std::string output) {
	std::vector<strvcfentry> entries;
	parse_GC(gc_file, entries);
	FILE *file;
	file = fopen(output.c_str(), "w");
	for (size_t i = 0; i < entries.size(); i++) {
		fprintf(file, "%s", print_entry_bio(entries[i]).c_str());
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}
