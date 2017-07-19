/*
 * Overlap_snps.cpp
 *
 *  Created on: Jul 18, 2017
 *      Author: sedlazec
 */

#include "Overlap_snps.h"
void parse_entry(char * buffer, std::string & name, int &len) {
	size_t i = 13;

	while (buffer[i] != ',') {
		name += buffer[i];
		i++;
	}
	while (buffer[i] != '=') {
		i++;
	}
	i++;
	len = atoi(&buffer[i]);
}
std::map<std::string, std::string> parse_vcf(std::string vcf_file) {
	std::map<std::string, std::string> genome;
	/*##contig=<ID=1,length=249250621>
	 ##contig=<ID=2,length=243199373>
	 ##contig=<ID=3,length=198022430>
	 ##contig=<ID=4,length=191154276>
	 ##contig=<ID=5,length=180915260>
	 ##contig=<ID=6,length=171115067>
	 ##contig=<ID=7,length=159138663>
	 ##contig=<ID=8,length=146364022>
	 ##contig=<ID=9,length=141213431>
	 ##contig=<ID=10,length=135534747>
	 ##contig=<ID=11,length=135006516>
	 ##contig=<ID=12,length=133851895>
	 ##contig=<ID=13,length=115169878>
	 ##contig=<ID=14,length=107349540>
	 ##contig=<ID=15,length=102531392>
	 ##contig=<ID=16,length=90354753>
	 ##contig=<ID=17,length=81195210>
	 ##contig=<ID=18,length=78077248>
	 ##contig=<ID=19,length=59128983>
	 ##contig=<ID=20,length=63025520>
	 ##contig=<ID=21,length=48129895>
	 ##contig=<ID=22,length=51304566>
	 ##contig=<ID=X,length=155270560>
	 ##contig=<ID=Y,length=59373566>
	 ##contig=<ID=MT,length=16569>
	 */
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SNP Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);

	while (!myfile.eof() && buffer[0] == '#') {
		if (strncmp("##contig", buffer, 8) == 0) {
			std::string name;
			int len;
			parse_entry(buffer, name, len);
			genome[name].assign(len, 'N');
		}
		myfile.getline(buffer, buffer_size);
	}
	return genome;
}
void overlap_snps(std::string svs_file, std::string snp_file, int max_dist, int min_svs, std::string output) {

	float min_eufreq = 0.05;

	Parameter::Instance()->min_freq = 0.05;   //TODO parameter?

	std::map<std::string, std::string> genome = parse_vcf(svs_file);

	std::vector<strvcfentry> entries = parse_vcf(svs_file, min_svs);
	std::cout << "merging entries: " << entries.size() << std::endl;
	for (size_t j = 0; j < entries.size(); j++) {
		int start = std::max(0, (int) entries[j].start.pos - max_dist);
		int stop = std::min((int) genome[entries[j].start.chr].size(), (int) entries[j].stop.pos + max_dist);

		if (entries[j].type == 0 || (entries[j].type == 1 || entries[j].type == 6)) {
			char sv = 'D';   //del
			if (entries[j].type == 1) {
				sv = 'P'; //DUP
			} else if (entries[j].type == 6) {
				sv = 'C'; //CNV
			}
			for (int pos = start; pos < stop; pos++) {
				genome[entries[j].start.chr][pos] = sv;
			}
		} else {
			char sv = 'V'; //INV
			if (entries[j].type == 3) {
				sv = 'T'; //TRA
			} else if (entries[j].type == 4) {
				sv = 'I'; //INS
			} else if (entries[j].type == 5) {
				sv = 'U';
			}

			for (int pos = start; pos < entries[j].start.pos + max_dist && pos < (int) genome[entries[j].start.chr].size(); pos++) {
				genome[entries[j].start.chr][pos] = sv;
			}

			for (int pos = stop; pos < entries[j].stop.pos + max_dist && pos < (int) genome[entries[j].stop.chr].size(); pos++) {
				genome[entries[j].stop.chr][pos] = sv;
			}
		}
	}

	std::cout << "fin filling" << std::endl;

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(snp_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SNP Parser: could not open file: " << snp_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);

	FILE *file;
	file = fopen(output.c_str(), "w");

	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			breakpoint_str snp;
			float freq = 0;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 0 && buffer[i] != '\t') {
					snp.chr += buffer[i];
					//std::cout << "to find" << snp.chr << " ";
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					snp.position = atoi(&buffer[i]);
					//std::cout << snp.position << std::endl;
				}
				if (count > 6 && strncmp("EUR_AF=", &buffer[i], 7) == 0) {
					freq = atof(&buffer[i + 7]);
					//	std::cout<<"fre: "<<freq<<std::endl;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}

			if (freq > min_eufreq) {
				fprintf(file, "%s", snp.chr.c_str());
				fprintf(file, "%c", '\t');
				fprintf(file, "%i", snp.position);
				fprintf(file, "%c", '\t');
				fprintf(file, "%f", freq);
				fprintf(file, "%c", '\t');

				if (genome[snp.chr][snp.position] == 'T') {
					fprintf(file, "%s", "TRA");
				} else if (genome[snp.chr][snp.position] == 'I') {
					fprintf(file, "%s", "INS");
				} else if (genome[snp.chr][snp.position] == 'U') {
					fprintf(file, "%s", "UNK");
				} else if (genome[snp.chr][snp.position] == 'D') {
					fprintf(file, "%s", "DEL");
				} else if (genome[snp.chr][snp.position] == 'P') {
					fprintf(file, "%s", "DUP");
				} else if (genome[snp.chr][snp.position] == 'C') {
					fprintf(file, "%s", "CNV");
				} else if (genome[snp.chr][snp.position] == 'V') {
					fprintf(file, "%s", "INV");
				} else {
					fprintf(file, "%s", "NA");
				}
				fprintf(file, "%c", '\n');
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	fclose(file);

}

std::vector<breakpoint_str> parse_chrs(char * buffer) {
	size_t i = 0;
	std::vector<breakpoint_str> snps;
	breakpoint_str tmp;
	tmp.position = -1;
	while (buffer[i] != '\t' && buffer[i] != ' ') {
		if (buffer[i] != ';') {
			tmp.chr += buffer[i];
		} else {
			snps.push_back(tmp);
			tmp.chr.clear();
		}
		i++;
	}
	snps.push_back(tmp);
	return snps;
}

void parse_pos(char *buffer, std::vector<breakpoint_str>&snps) {
	size_t i = 0;
	size_t id = 1;
	breakpoint_str tmp;
	snps[0].position = atoi(&buffer[i]);
	while (buffer[i] != '\t' && buffer[i] != ' ') {
		if (buffer[i - 1] == ';') {
			snps[id].position = atoi(&buffer[i]);
			id++;
		}
		i++;
	}
	snps[id].position = atoi(&buffer[i]);
}
void overlap_snps_gwas(std::string svs_file, std::string snp_file, int max_dist, int min_svs, std::string output) {

	Parameter::Instance()->min_freq = 0.05;   //TODO parameter?
	std::map<std::string, std::string> genome = parse_vcf(svs_file);

	std::vector<strvcfentry> entries = parse_vcf(svs_file, min_svs);
	std::cout << "merging entries: " << entries.size() << std::endl;
	for (size_t j = 0; j < entries.size(); j++) {
		int start = std::max(0, (int) entries[j].start.pos - max_dist);
		int stop = std::min((int) genome[entries[j].start.chr].size(), (int) entries[j].stop.pos + max_dist);

		if (entries[j].type == 0 || (entries[j].type == 1 || entries[j].type == 6)) {
			char sv = 'D';   //del
			if (entries[j].type == 1) {
				sv = 'P'; //DUP
			} else if (entries[j].type == 6) {
				sv = 'C'; //CNV
			}
			for (int pos = start; pos < stop; pos++) {
				genome[entries[j].start.chr][pos] = sv;
			}
		} else {
			char sv = 'V'; //INV
			if (entries[j].type == 3) {
				sv = 'T'; //TRA
			} else if (entries[j].type == 4) {
				sv = 'I'; //INS
			} else if (entries[j].type == 5) {
				sv = 'U';
			}

			for (int pos = start; pos < entries[j].start.pos + max_dist && pos < (int) genome[entries[j].start.chr].size(); pos++) {
				genome[entries[j].start.chr][pos] = sv;
			}

			for (int pos = stop; pos < entries[j].stop.pos + max_dist && pos < (int) genome[entries[j].stop.chr].size(); pos++) {
				genome[entries[j].stop.chr][pos] = sv;
			}
		}
	}

	std::cout << "fin filling" << std::endl;

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(snp_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SNP Parser: could not open file: " << snp_file.c_str() << std::endl;
		exit(0);
	}
	FILE *file;
	file = fopen(output.c_str(), "w");

	myfile.getline(buffer, buffer_size);
	fprintf(file, "%s", buffer);
	fprintf(file, "%c", '\t');
	fprintf(file, "%s","SVs");
	fprintf(file, "%c", '\n');
	myfile.getline(buffer, buffer_size);

	while (!myfile.eof()) {

		int count = 0;
		std::vector<breakpoint_str> snps;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 11 && buffer[i - 1] == '\t') {
				snps = parse_chrs(&buffer[i]);
			}
			if (count == 12 && buffer[i - 1] == '\t') {
				parse_pos(&buffer[i], snps);
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		for (size_t i = 0; i < snps.size(); i++) {
			if (snps[i].position > 0) {
				breakpoint_str snp = snps[i];
				fprintf(file, "%s", buffer);
			//	fprintf(file, "%c", '\t');
			//	fprintf(file, "%i", snp.position);
			//	fprintf(file, "%c", '\t');
			//	fprintf(file, "%i", 0);
				fprintf(file, "%c", '\t');

				if (genome[snp.chr][snp.position] == 'T') {
					fprintf(file, "%s", "TRA");
				} else if (genome[snp.chr][snp.position] == 'I') {
					fprintf(file, "%s", "INS");
				} else if (genome[snp.chr][snp.position] == 'U') {
					fprintf(file, "%s", "UNK");
				} else if (genome[snp.chr][snp.position] == 'D') {
					fprintf(file, "%s", "DEL");
				} else if (genome[snp.chr][snp.position] == 'P') {
					fprintf(file, "%s", "DUP");
				} else if (genome[snp.chr][snp.position] == 'C') {
					fprintf(file, "%s", "CNV");
				} else if (genome[snp.chr][snp.position] == 'V') {
					fprintf(file, "%s", "INV");
				} else if (genome[snp.chr][snp.position] == 'N') {
					fprintf(file, "%s", "NA");
				} else {
					fprintf(file, "%s", "REP");
				}
				genome[snp.chr][snp.position] = 'R';
				fprintf(file, "%c", '\n');
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	fclose(file);

}

