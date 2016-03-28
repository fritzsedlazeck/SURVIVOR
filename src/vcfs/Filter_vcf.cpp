/*
 * Filter_vcf.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */
#include "Filter_vcf.h"
//parses all genotype calls and summarizes them.
void translate(const char * entry, int min_alternative_pairs,
		float min_alt_ref_ratio, int max_genotype, strentry & result) {
//	GT:GL:GQ:FT:RC:DR:DV:RR:RV
	size_t i = 0;
	int count = 0;
	float dv = 0; //"# high-quality variant pairs"
	float dr = 0; //"# high-quality reference pairs"
	int gq = 0; //"Genotype Quality"
	int rr = 0;
	int rv = 0;
	int lumpy = 0;
	//std::cout<<entry<<std::endl;
	while (entry[i] != '\t' && entry[i] != '\n' && entry[i] != '\0') {
		if (count == 2 && entry[i - 1] == ':') {
			gq = atoi(&entry[i]);
		}
		if (count == 5 && entry[i - 1] == ':') {
			dr = atoi(&entry[i]);
		}
		if (count == 6 && entry[i - 1] == ':') {
			dv = atoi(&entry[i]);
		}
		if (count == 7 && entry[i - 1] == ':') {
			rr = atoi(&entry[i]);
		}
		if (count == 8 && entry[i - 1] == ':') {
			rv = atoi(&entry[i]);
		}
		if (count == 9 && entry[i - 1] == ':') {
			lumpy = atoi(&entry[i]);
		}
		if (entry[i] == ':') {
			count++;
		}
		i++;
	}
//	std::cout<<dv<<" "<<dv/dr<<" "<<gq<<std::endl;
	if (dv == 0 && dr == 0) {
		result.not_covered++;
	} else if (dv + rv < min_alternative_pairs) {
		result.not_valid++;
	} else if (dr > 0 && (dv + rv) / (dr+rr) < min_alt_ref_ratio) {
		result.not_valid++;
	} else if (gq > max_genotype) {
		result.not_valid++;
	} else {
		//std::cout<<"valid"<<std::endl;
		result.valid++;
	}
}

//check if SV should be kept:
bool pass_filter(strentry entry, strregion region,
		std::vector<strregion> ignore_regions) {
	//TODO: compare to bed file to filter out regions!
	for (size_t i = 0; i < ignore_regions.size(); i++) {
		if (strcmp(region.start.chr.c_str(),
				ignore_regions[i].start.chr.c_str()) == 0) {
			if (region.start.pos > ignore_regions[i].start.pos
					&& region.start.pos < ignore_regions[i].stop.pos) {
				return false;
			}
		}

		if (strcmp(region.stop.chr.c_str(), ignore_regions[i].start.chr.c_str())
				== 0) {
			if (region.stop.pos > ignore_regions[i].start.pos
					&& region.stop.pos < ignore_regions[i].stop.pos) {
				return false;
			}
		}

	}
	if (entry.valid > 0) {
		return true;
	}
	return false;
}

//parse the bed file that defines regions that should be ignored:
std::vector<strregion> parse_bed(std::string filename) {
	std::vector<strregion> ignore_regions;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "BED Parser: could not open file: " << filename.c_str()
				<< std::endl;
		return ignore_regions;
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		int count = 0;
		strregion tmp;
		for (size_t i = 0;
				i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n';
				i++) {
			if (count == 0 && buffer[i] != '\t') {
				tmp.start.chr += buffer[i];
				tmp.stop.chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				tmp.start.pos = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				tmp.stop.pos = atoi(&buffer[i]);
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		ignore_regions.push_back(tmp);
		myfile.getline(buffer, buffer_size);
	}

	myfile.close();
	return ignore_regions;
}

void filter_vcf(std::string vcf_file, std::string genomic_regions,
		int min_alternative_pairs, float min_alt_ref_ratio, int max_genotype,
		std::string outputvcf) {

	std::vector<strregion> ignore_regions = parse_bed(genomic_regions);
	//parse vcf file and write vcf file in one go:

	std::vector<std::string> names;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: "
				<< vcf_file.c_str() << std::endl;
		exit(0);
	}
	FILE *file;
	file = fopen(outputvcf.c_str(), "w");

	myfile.getline(buffer, buffer_size);
	int deleted = 0;
	while (!myfile.eof()) {
		if (buffer[0] == '#') { //write header info:
			fprintf(file, "%s", buffer);
			fprintf(file, "%c", '\n');
		} else {
			int count = 0;
			strregion region;
			strentry entry;
			entry.not_covered = 0;
			entry.not_valid = 0;
			entry.valid = 0;
			for (size_t i = 0;
					i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n';
					i++) {
				if (count == 0 && buffer[i] != '\t') {
					region.start.chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					region.start.pos = atoi(&buffer[i]);
				}
				if (count == 7 && strncmp(&buffer[i], "CHR2=", 5) == 0) {
					region.stop = parse_stop(&buffer[i]);
				}
				if (count >= 9 && buffer[i - 1] == '\t') {
					translate(&buffer[i], min_alternative_pairs,
							min_alt_ref_ratio, max_genotype, entry);
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}

			if (pass_filter(entry, region, ignore_regions)) {
				std::string tmp = std::string(buffer);

				std::size_t found = tmp.find("\tLowQual\t");
				if (found != std::string::npos) {
					tmp.replace(tmp.find("\tLowQual\t"), 9, "\tPASS\t");
				}

				fprintf(file, "%s", tmp.c_str());
				fprintf(file, "%c", '\n');
			} else {
				deleted++;
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	std::cout << "WE deleted: " << deleted << std::endl;
	myfile.close();
	fclose(file);
}
