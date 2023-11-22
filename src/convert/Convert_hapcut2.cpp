/*
 * Convert_hapcut2.cpp
 *
 *  Created on: Mar 14, 2018
 *      Author: sedlazec
 */

#include "Convert_hapcut2.h"

map<std::string, int> parse_hapcut(std::string hapcut2, std::string target_chr, int & phaseblock_id) {
	std::string buffer;
	std::ifstream myfile;

	map<std::string, int> hapcut;

	myfile.open(hapcut2.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Hapcut2 Parser: could not open file: " << hapcut2.c_str() << std::endl;
		exit(0);
	}

	getline(myfile, buffer);
	while (!myfile.eof()) {
		if (buffer[0] != 'B' && buffer[0] != '*') { //avoid headers!
			int count = 0;
			std::string chr = "";
			int pos = -1;
			bool first_gt = true;
			bool second_gt = true; //not needed but good to check!
			for (size_t i = 0; i < buffer.size(); i++) {
				if (count == 1 && buffer[i - 1] == '\t') {
					first_gt = (bool) (buffer[i] != '0');
				}
				if (count == 2 && buffer[i - 1] == '\t') {
					second_gt = (bool) (buffer[i] != '0');
				}

				if (count == 3 && buffer[i] != '\t') {
					chr += buffer[i];
				}
				if (count == 4 && buffer[i - 1] == '\t') {
					pos = atoi(&buffer[i]);
					break;
				}

				if (buffer[i] == '\t') {
					count++;
				}
			}
			if (strcmp(chr.c_str(), target_chr.c_str()) == 0) {
				if ((first_gt && !second_gt) || (!first_gt && second_gt)) {
					std::stringstream ss;
					ss << chr;
					ss << "_";
					ss << pos;

					if (hapcut.find(ss.str()) == hapcut.end()) {
						hapcut[ss.str()] = phaseblock_id;
						if (!first_gt) { //negative ID for none
							hapcut[ss.str()] = hapcut[ss.str()] * -1;
						}
					} else {
						cerr << "A position was found twice: " << ss.str() << endl;
					}
				}
			}
		} else if(buffer[0] != 'B') {
			phaseblock_id++;
		}
		getline(myfile, buffer);
	}
	myfile.close();
	return hapcut;
}

void process_hapcut(std::string orig_snp, std::string hapcut2, std::string output) {

	//parse VCF file. just if we dected a 0/1 we check the hapcut results.
	std::string buffer;
	std::ifstream myfile;

	myfile.open(orig_snp.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SNP Parser: could not open file: " << orig_snp.c_str() << std::endl;
		exit(0);
	}
	FILE * file = fopen(output.c_str(), "w");
	map<std::string, int> hapcut_res;
	std::string old_chr;
	getline(myfile, buffer);
	int phaseblock = 1;
	//int num = 0;
	while (!myfile.eof()) {
		if (buffer[0] == '#') {
			if (buffer[1] == 'C') {
				fprintf(file, "%s", "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n");
			}
			fprintf(file, "%s", buffer.c_str());
			fprintf(file, "%c", '\n');
		} else {
			//	num++;

			std::size_t found = buffer.find_last_of('\t');
			//check found + 1 if '1' && if found+3 !='0' ->
			if ((buffer[found + 1] == '0' && buffer[found + 3] == '1') || (buffer[found + 1] == '1' && buffer[found + 3] == '0')) {
				//search pos and replace 3 chars.
				int count = 0;
				std::string chr = "";
				int pos = 0;
				for (size_t i = 0; i < buffer.size(); i++) {
					if (count == 0 && buffer[i] != '\t') {
						chr += buffer[i];
					}
					if (count == 1 && buffer[i - 1] == '\t') {
						pos = atoi(&buffer[i]);
						break;
					}
					if (buffer[i] == '\t') {
						count++;
					}
				}
				if (strcmp(chr.c_str(), old_chr.c_str()) != 0) {
					//load new chr set:
					cout << "Parsing hapcut2 output for " << chr;
					hapcut_res = parse_hapcut(hapcut2, chr, phaseblock);
					cout << " SNPs parsed " << hapcut_res.size() << endl;
					old_chr = chr;
				}

				if (!chr.empty()) {

					std::stringstream ss;
					ss << chr;
					ss << "_";
					ss << pos;
					if (hapcut_res.find(ss.str()) != hapcut_res.end()) {

						buffer.insert(found, ":PS");
						found += 3;
						//			cout << "MATCH: " << ss.str() << endl;
						if (hapcut_res[ss.str()] > 0) {
							buffer[found + 1] = '1';
							buffer[found + 2] = '|';
							buffer[found + 3] = '0';
						} else {
							buffer[found + 1] = '0';
							buffer[found + 2] = '|';
							buffer[found + 3] = '1';
						}
						std::stringstream id;
						id << ":";
						id << abs(hapcut_res[ss.str()]);

						buffer.append(id.str());

					}
				}

			}
			fprintf(file, "%s", buffer.c_str());
			fprintf(file, "%c", '\n');

		}

		getline(myfile, buffer);
	}
	myfile.close();
	fclose(file);
}
