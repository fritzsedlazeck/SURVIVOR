/*
 * MT_identifier.cpp
 *
 *  Created on: Aug 15, 2017
 *      Author: sedlazec
 */

#include "MT_identifier.h"
int parse_read_start(std::string cigar, bool strand) {
	std::string::size_type sz;   // alias of size_t
	int pos = 0;
	if (strand) {
		size_t found = cigar.find_first_of("S");
		pos = std::stoi(cigar.substr(0, found), &sz);
	} else {
		std::string tmp="";
		std::string::reverse_iterator i = cigar.rbegin();
		while (i != cigar.rend() && (*i) != 'M') {
			i++;
			tmp+=(*i);
		}
		cout<<"Check Cigar backend:" << tmp<<std::endl;
		pos = std::stoi(tmp, &sz);
	}
	return pos;
}

void insert_sort(segment_str tmp, std::vector<segment_str> & segments) {

	size_t i = 0;
	while (i < segments.size()) {
	//	if (tmp.pos<)
			i++;
	}
}
void parse_entries(const char * buffer, std::vector<segment_str> & segments) {
	size_t i = 5; //avoid SA:Z:
	int count = 0;
	segment_str tmp;
	std::string cigar = "";
	while (buffer[i] != '\t') {
		if (count == 0 && buffer[i] != ',') {
			tmp.chr += buffer[i];
		}
		if (count == 1 && buffer[i - 1] == ',') {
			tmp.pos = atoi(&buffer[i]);
		}
		if (count == 2 && buffer[i - 1] == ',') {
			tmp.strand = (buffer[i] == '+');
		}
		if (count == 3 && buffer[i] != ',') {
			cigar += buffer[i];
		}
		if (count == 4 && buffer[i - 1] == ',') {
			//parse cigar and set read start:
			tmp.read_start=parse_read_start(cigar, tmp.strand);
			tmp.MQ = atoi(&buffer[i]);
		}
		if (buffer[i] == ';') {
			//store;
			insert_sort(tmp, segments);
			count = 0;
		}
		if (buffer[i] == ',') {
			count++;
		}
	}
}
void detect_MT_copies(std::string chr_identifier) {
	//Q1: Start lockations always full lenght?
	//Q2: Avg vs. max copies per read

	chr_identifier = "MT";

	std::vector<int> start_pos; //we can do that interactively...

	int min_len = 0;
	while (!cin.eof()) {
		string line;
		getline(cin, line);
		if (!cin.fail()) {
			size_t found = 0;
			found = line.find_first_of(chr_identifier);
			if (found != std::string::npos) { //only if the line includes MT!
				if (line[0] == '@') {
					//get the length of the chr!
					//found: @SQ     SN:MT   LN:19431
					found = line.find_first_of("LN:");
					found += 3;
					std::string::size_type sz;   // alias of size_t
					int len = std::stoi(line.substr(found), &sz);   //get chr size
					min_len = len * 2;
					std::cout << "LEN: " << min_len << std::endl;
				} else {
					int count = 0;
					int sequence = 0;
					std::string cigar = "";
					std::vector<segment_str> segments;
					segment_str tmp;
					for (size_t i = 0; i < line.size(); i++) {
						if (count == 2 && line[i] != '\t') {
							tmp.chr += line[i];
						}
						if (count == 1 && line[i - 1] == '\t') {
							tmp.strand = (line[i] == '0'); //should be 0 (+) or 16 (-)
						}

						//parse cigar and read start ;
						if (count == 5 && line[i] != '\t') {
							cigar += line[i];
						}
						if (count == 10 && line[i] != '\t') {
							sequence++;
						}
						if (count == 11 && line[i - 1] == '\t') {

							if (sequence > min_len || strncmp(chr_identifier.c_str(), tmp.chr.c_str(), chr_identifier.size()) != 0) {
								sequence = 0;
								break; //early terminate we dont need to parse the rest.
							}
							//parse cigar and read start ;
							tmp.read_start = parse_read_start(cigar, tmp.strand);
							segments.push_back(tmp);
						}
						if (count == 20 && line[i] == '\t') {
							parse_entries(line.substr(i).c_str(), segments);
						}

						if (line[i] == '\t') {
							count++;
						}
					}
					if (sequence > min_len) {

					}
				}
			}
		} else {
			break;
		}
	}

}
