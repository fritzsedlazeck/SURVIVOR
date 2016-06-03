/*
 * Compoverla_VCF.cpp
 *
 *  Created on: Feb 27, 2015
 *      Author: fsedlaze
 */

#include "Compoverlap_VCF.h"

int overlap_bothdir(strvcfentry vcf1, std::vector<strvcfentry> vcf2, int max_dist) {
	for (size_t i = 0; i < vcf2.size(); i++) {
		//check type:
		if (vcf2[i].type == vcf1.type) {
			//check chrs:
			if (strcmp(vcf2[i].stop.chr.c_str(), vcf1.stop.chr.c_str()) == 0 && strcmp(vcf2[i].start.chr.c_str(), vcf1.start.chr.c_str()) == 0) {
				//check coordinates:
				if (abs(vcf2[i].stop.pos - vcf1.stop.pos) < max_dist && abs(vcf2[i].start.pos - vcf1.start.pos) < max_dist) {
					return i;
				}
			}

			if (strcmp(vcf2[i].start.chr.c_str(), vcf1.stop.chr.c_str()) == 0 && strcmp(vcf2[i].stop.chr.c_str(), vcf1.start.chr.c_str()) == 0) {
				//check coordinates:
				if (abs(vcf2[i].start.pos - vcf1.stop.pos) < max_dist && abs(vcf2[i].stop.pos - vcf1.start.pos) < max_dist) {
					return i;
				}
			}
		}
	}
	return -1;
}

void print_header(std::string vcf_file, FILE *& out) {
	std::vector<strsimul> simulated;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "BED Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if (buffer[0] == '#') {
			fprintf(out, "%s", buffer);
			fprintf(out, "%c", '\n');
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
}
void print_entry(strvcfentry entry, FILE *& out) {
	std::string tmp = entry.header;
	int count = 0;
	for (size_t i = 0; i < tmp.size(); i++) {
		if (count == 7&&tmp[i-1]=='\t') {
			if(!entry.caller_supports.empty()){
				fprintf(out, "%s", "SUP=");
				fprintf(out, "%i", entry.caller_supports[0]);
				for (size_t j = 1; j < entry.caller_supports.size(); j++) {
					fprintf(out, "%c", ',');
					fprintf(out, "%i", entry.caller_supports[j]);
				}
				fprintf(out, "%c", ';');
			}
		}
		fprintf(out, "%c", tmp[i]);
		if (tmp[i] == '\t') {
			count++;
		}
	}
	//fprintf(out, "%s", entry.header.c_str());
	for (std::map<std::string, std::string>::iterator tz = entry.calls.begin(); tz != entry.calls.end(); tz++) {
		fprintf(out, "%s", (*tz).second.c_str());
	}
	fprintf(out, "%c", '\n');
}

void comp_overlap_vcf(std::string vcf1_file, std::string vcf2_file, int max_dist, std::string output) {
	std::vector<strvcfentry> vcf1 = parse_vcf(vcf1_file);
	std::vector<strvcfentry> vcf2 = parse_vcf(vcf2_file);

	std::cout << vcf1.size() << " " << vcf2.size() << std::endl;

	FILE * combined;
	FILE * unique_lumpy;
	FILE * unique_delly;
	std::string out = output;
	out += "_overlap.vcf";
	combined = fopen(out.c_str(), "w");

	out = output;
	out += "_uniq_delly.vcf";
	unique_delly = fopen(out.c_str(), "w");

	out = output;
	out += "_uniq_lumpy.vcf";
	unique_lumpy = fopen(out.c_str(), "w");

	print_header(vcf1_file, combined);
	print_header(vcf1_file, unique_delly);
	print_header(vcf2_file, unique_lumpy);

	for (size_t i = 0; i < vcf1.size(); i++) {
		vcf1[i].sup_lumpy = 0;
	}
	for (size_t i = 0; i < vcf2.size(); i++) {
		vcf2[i].sup_lumpy = 0;
	}
	for (size_t i = 0; i < vcf1.size(); i++) {
		int id = overlap_bothdir(vcf1[i], vcf2, max_dist);
		if (id > -1) {
			vcf1[i].sup_lumpy = 1;
			vcf2[id].sup_lumpy = 1;
		}
	}

	std::vector<int> tmp;
	tmp.resize(5, 0);
	std::vector<int> overlap_vcf1 = tmp;
	std::vector<int> overlap_vcf2 = tmp;
	std::vector<int> unique_vcf1 = tmp;
	std::vector<int> unique_vcf2 = tmp;

	std::cout << "vcf1.size " << vcf1.size() << std::endl;
	for (size_t i = 0; i < vcf1.size(); i++) {
		if (vcf1[i].sup_lumpy == 1) {
			overlap_vcf1[vcf1[i].type]++;
			overlap_vcf2[vcf2[i].type]++;
			//	print_entry( vcf1[i],combined);
		} else {
			unique_vcf1[vcf1[i].type]++;
			//	print_entry( vcf1[i],unique_delly);
		}
	}
	for (size_t i = 0; i < vcf2.size(); i++) {
		if (vcf2[i].sup_lumpy == 0) {
			unique_vcf2[vcf2[i].type]++;
			//std::cout << "Not " << trans_type(vcf2[i].type) << " "
			//		<< vcf2[i].start.chr << " " << vcf2[i].start.pos << " "
			//		<< vcf2[i].stop.chr << " " << vcf2[i].stop.pos << std::endl;
			//	print_entry( vcf2[i],unique_lumpy);
		}
	}
	fclose(unique_delly);
	fclose(unique_lumpy);
	fclose(combined);

//0=DEL,1=DUP,2=INV,3=TRA
	std::cout << "Overlap VCF1:" << " DEL " << overlap_vcf1[0] << " DUP " << overlap_vcf1[1] << " INV " << overlap_vcf1[2] << " TRA " << overlap_vcf1[3] << std::endl;
	std::cout << "Overlap VCF2:" << " DEL " << overlap_vcf2[0] << " DUP " << overlap_vcf2[1] << " INV " << overlap_vcf2[2] << " TRA " << overlap_vcf2[3] << std::endl;
	std::cout << "Uniqe VCF1:" << " DEL " << unique_vcf1[0] << " DUP " << unique_vcf1[1] << " INV " << unique_vcf1[2] << " TRA " << unique_vcf1[3] << std::endl;
	std::cout << "Uniqe VCF2:" << " DEL " << unique_vcf2[0] << " DUP " << unique_vcf2[1] << " INV " << unique_vcf2[2] << " TRA " << unique_vcf2[3] << std::endl;

}
