/*
 * Merge_VCF.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */

#include "Merge_VCF.h"

//read in all the vcf filenames:
std::vector<std::string> parse_filename(std::string filename) {
	std::vector<std::string> names;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: "
				<< filename.c_str() << std::endl;
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

strcoordinate parse_stop(const char * buffer) {
	size_t i = 0;
	bool chr_flag = false;
	strcoordinate pos;
	pos.chr="";
	while (buffer[i] != '\t') {
		if (strncmp(&buffer[i],";END=",5)==0 ){
			pos.pos = atoi(&buffer[i+5]);
		}
		if((strncmp(&buffer[i],"END=",4)==0&&i==0)){
			pos.pos = atoi(&buffer[i+4]);
		}
		if (strncmp(&buffer[i],"CHR2=",5)==0){
			i=i+5;
			chr_flag=true;
		}
		if(buffer[i]==';'){
			chr_flag=false;
		}
		if(chr_flag){
			pos.chr += buffer[i];
		}

		i++;
	}
	//std::cout<<"end: "<<pos.chr<<" "<<pos.pos<<std::endl;
	return pos;
}
std::pair <bool,bool>parse_strands(const char * buffer) {
	std::pair<bool,bool> strands;
	size_t i = 0;
	while (buffer[i] != '\t') {
		if (strncmp(&buffer[i],"3to5",4)==0){
			strands.first=false;
			strands.second=true;
			return strands;
		}
		if (strncmp(&buffer[i],"3to3",4)==0){
			strands.first=false;
			strands.second=false;
			return strands;
		}
		if (strncmp(&buffer[i],"5to3",4)==0){
			strands.first=true;
			strands.second=false;
			return strands;
		}
		if (strncmp(&buffer[i],"5to5",4)==0){
			strands.first=true;
			strands.second=true;
			return strands;
		}
		i++;
	}
	return strands;
}

short get_type(std::string type) {
	if (strncmp(type.c_str(), "DEL", 3) == 0) {
		return 0;
	} else if (strncmp(type.c_str(), "DUP", 3) == 0) {
		return 1;
	} else if (strncmp(type.c_str(), "INV", 3) == 0) {
		return 2;
	} else if (strncmp(type.c_str(), "TRA", 3) == 0) {
		return 3;
	} else if (strncmp(type.c_str(), "INS", 3) == 0) {
		return 4;
	} else {
		std::cerr << "Unknown type!" <<type<< std::endl;

	}
	return -1;
}
//for each file parse the entries
std::vector<strvcfentry> parse_vcf(std::string filename) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: "
				<< filename.c_str() << std::endl;
		exit(0);
	}
	std::vector<strvcfentry> calls;
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			strvcfentry tmp;
			tmp.sup_lumpy=0;
			//std::cout<<buffer<<std::endl;
			for (size_t i = 0;
					i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n';
					i++) {

				if (count == 0 && buffer[i] != '\t') {
					tmp.start.chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					tmp.start.pos = atoi(&buffer[i]);
					//std::cout<<tmp.start.pos<<std::endl;
				}
				if (count == 7 && buffer[i - 1] == '\t') {
					tmp.stop = parse_stop(&buffer[i]);
					tmp.strands=parse_strands(&buffer[i]);
					//std::cout<<tmp.stop.chr<<std::endl;
				}
				if (count == 4 && buffer[i - 1] == '<') {

					tmp.type = get_type(std::string(&buffer[i]));
				}
				if (count == 9 && buffer[i - 1] == '\t') {
					tmp.calls[filename] = std::string(&buffer[i]);
					//std::cout<<std::string(&buffer[i])<<std::endl;
					//std::cout<<tmp.calls[filename]<<std::endl;
					break;
				}

				if (count >= 0 && count < 9) {
					tmp.header += buffer[i];
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			calls.push_back(tmp);
			tmp.calls.clear();
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	std::cout<<calls.size()<<std::endl;
	return calls;
}

int overlap(strvcfentry tmp, std::vector<strvcfentry> & final_vcf,
		int max_dist) {
	for (size_t i = 0; i < final_vcf.size(); i++) {
		//check type:
		if (final_vcf[i].type == tmp.type) {
			//check chrs:
			if (strcmp(final_vcf[i].stop.chr.c_str(), tmp.stop.chr.c_str()) == 0
					&& strcmp(final_vcf[i].start.chr.c_str(),
							tmp.start.chr.c_str()) == 0) {
				//check coordinates:
				if (abs(final_vcf[i].stop.pos - tmp.stop.pos) < max_dist
						&& abs(final_vcf[i].start.pos - tmp.start.pos)
								< max_dist) {
					return i;
				}
			}
		}
	}
	return -1;
}

//detect overlap and merge:
void merge_entries(std::string filename, int max_dist,
		std::vector<strvcfentry> & final_vcf) {
	//get new entries
	std::vector<strvcfentry> new_entries = parse_vcf(filename);
	//merge entires:
	for (size_t i = 0; i < new_entries.size(); i++) {
		int id = overlap(new_entries[i], final_vcf, max_dist);
		if (id > -1) {
			//std::cout<<"add " <<new_entries[i].calls[filename]<<std::endl;
			final_vcf[id].calls[filename] = new_entries[i].calls[filename]; //add call to entries;
		} else {
			//std::cout<<"push " <<new_entries[i].calls[filename]<<std::endl;
			final_vcf.push_back(new_entries[i]);
		}
	}
}

std::string get_header(std::vector<std::string> names) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(names[0].c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: "
				<< names[0].c_str() << std::endl;
		exit(0);
	}
	std::string header;
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if (buffer[0] == '#' && buffer[1] == '#') {
			header += std::string(buffer);
			header += '\n';
		} else if (buffer[0] == '#') {
			int count = 0;
			for (size_t i = 0;
					i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n';
					i++) {
				if (count < 9) {
					header += buffer[i];
				}
				if (count == 9) {
					break;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			break;
		}
		myfile.getline(buffer, buffer_size);
	}

	for (size_t i = 0; i < names.size(); i++) {
		header += '\t';
		header += names[i];
	}
	header += '\n';
	myfile.close();
	return header;
}

void print_merged_vcf(std::string outputfile, std::string header,
		std::vector<strvcfentry> &final_vcf, std::vector<std::string> names) {
	FILE *file;
	file = fopen(outputfile.c_str(), "w");

	fprintf(file, "%s", header.c_str());
	header.clear();

	for (size_t i = 0; i < final_vcf.size(); i++) {
		fprintf(file, "%s", final_vcf[i].header.c_str());
		for (size_t t = 0; t < names.size(); t++) {
			fprintf(file, "%c", '\t');
			if (final_vcf[i].calls.find(names[t]) != final_vcf[i].calls.end()) { //found an entry
				fprintf(file, "%s", final_vcf[i].calls[names[t]].c_str());
			} else {
				fprintf(file, "%s", "./.:0,0.0,0.0:0:NotDetected:0:0:0:0:0");
			}
		}
		fprintf(file, "%c", '\n');
	}
}

//main:
void merge_vcf(std::string filenames, int max_dist, std::string outputfile) {

	std::vector<std::string> names = parse_filename(filenames);
	std::cout << "found in file: " << names.size() << std::endl;
	std::vector<strvcfentry> final_vcf;

	for (size_t i = 0; i < names.size(); i++) {
		merge_entries(names[i], max_dist, final_vcf);
		std::cout << "merged: " << final_vcf.size() << std::endl;
	}
	std::cout << "get header:" << std::endl;
	std::string header = get_header(names);
	std::cout << "print:" << std::endl;
	print_merged_vcf(outputfile, header, final_vcf, names);
}
