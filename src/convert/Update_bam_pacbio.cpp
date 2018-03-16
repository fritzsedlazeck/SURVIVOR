/*
 * Update_bam_pacbio.cpp
 *
 *  Created on: Mar 15, 2018
 *      Author: sedlazec
 */

#include "Update_bam_pacbio.h"


std::vector<std::string> parse_header(std::string unmapped_sam){
	std::vector<std::string> header;


	return header;
}



void merge_header(std::string unmapped_sam,std::string mapped_sam,FILE *& file2) {
	std::vector<std::string> header_un=parse_header(unmapped_sam);


}

void update_entries(std::map<std::string, std::string> & entries, std::string unmapped_sam) {
	std::cout<<"check unmapped"<<std::endl;
	std::string buffer;
	std::ifstream myfile;
	myfile.open(unmapped_sam.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Sam Parser: could not open file: " << unmapped_sam.c_str() << std::endl;
		exit(0);
	}
	getline(myfile, buffer);
	while (!myfile.eof()) {
		if (buffer[0] != '@') {
			size_t found = buffer.find_first_of('\t');
			std::string id = buffer.substr(0, found);
			if (entries.find(id) != entries.end()) { //found!
				int count = 0;
				for (size_t i = 0; i < buffer.size(); i++) {
					if (count > 9) {
						entries[id] += buffer[i];
					}
					if (buffer[i] == '\t') {
						count++;
					}
				}
			}
		}
		getline(myfile, buffer);
	}

}

void process_sam_forpacbio(std::string unmapped_sam, std::string mapped_sam, std::string output_sam) {


	std::string buffer;
	std::ifstream myfile;
	myfile.open(mapped_sam.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Sam Parser: could not open file: " << mapped_sam.c_str() << std::endl;
		exit(0);
	}

	FILE *file2;
	file2 = fopen(output_sam.c_str(), "w");

	merge_header(unmapped_sam,mapped_sam,file2);

	std::map<std::string, std::string> entries;
	getline(myfile, buffer);
	while (!myfile.eof()) { //avoid header.
		if (buffer[0] != '@') {
			//parse part of the mapped entries into a map (e.g. step size =100000)
			size_t found = buffer.find_first_of('\t');
			std::string id = buffer.substr(0, found);
			entries[id] = buffer;
			if (entries.size() > 1000) {
				std::cout<<"check entries"<<std::endl;
				//check orig file and update them
				update_entries(entries, unmapped_sam);
				for (std::map<std::string, std::string>::iterator i = entries.begin(); i != entries.end(); i++) {
					fprintf(file2, "%s", (*i).second.c_str());
					fprintf(file2, "%c", '\n');
				}
				entries.clear();
			}
		}

		getline(myfile, buffer);
	}
	myfile.close();
	//check orig file and update them
	update_entries(entries, unmapped_sam);
	for (std::map<std::string, std::string>::iterator i = entries.begin(); i != entries.end(); i++) {
		fprintf(file2, "%s", (*i).second.c_str());
		fprintf(file2, "%c", '\n');
	}
	fclose(file2);


}

