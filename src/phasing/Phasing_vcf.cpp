/*
 * Phasing_vcf.cpp
 *
 *  Created on: Sep 26, 2018
 *      Author: sedlazec
 */

#include "Phasing_vcf.h"
std::vector<snp_str> parse_hapcut2(std::string hapcut_output) {
	std::vector<snp_str> snps;

	std::string buffer;
	std::ifstream myfile;
	myfile.open(hapcut_output.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Hapcut Parser: could not open file: " << hapcut_output.c_str() << std::endl;
		exit(0);
	}
	getline(myfile, buffer);
	int phase_block_id = 0;
	while (!myfile.eof()) {
		if (buffer[0] != 'B') {

			if (buffer[0] == '*') {
				//store new
				phase_block_id++;
			} else {
				snp_str tmp;
				tmp.parental = 0;
				tmp.gatk = 0;
				tmp.phase_block = phase_block_id;
				tmp.ratio = -1;
				//parse the actual snps:
				int count = 0;
				for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
					if (count == 1 && buffer[i] == '1') {
						tmp.haplotype = true;
					}
					if (count == 2 && buffer[i] == '1') {
						tmp.haplotype = false;
					}
					if (count == 3 && buffer[i] != '\t') {
						tmp.chr += buffer[i];
					}
					if (count == 4 && buffer[i - 1] == '\t') {
						tmp.position = atoi(&buffer[i]);
					}
					if (count == 6 && buffer[i] != '\t') {
						tmp.alt_allele = buffer[i];
						break;
					}
					if (buffer[i] == '\t') {
						count++;
					}
				}
				if (strcmp(tmp.chr.c_str(), "17") == 0) {
					snps.push_back(tmp);

				}

			}

		}
		getline(myfile, buffer);
	}

	myfile.close();
	return snps;
}

void update_parents_xatlas(std::string parents_vcf, std::vector<snp_str> &snps) {
	std::string buffer;
	std::ifstream myfile;
	myfile.open(parents_vcf.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Hapcut Parser: could not open file: " << parents_vcf.c_str() << std::endl;
		//	exit(0);
	} else {
		getline(myfile, buffer);
		while (!myfile.eof()) {
			int pos = 0;
			std::string chr;
			char alt_allele = ' ';
			short parental = -1; //0=na ; 1=father; 2=mother;
			int count = 0;
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 0 && buffer[i] != '\t') {
					chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					pos = atoi(&buffer[i]);
				}
				if (count == 4 && buffer[i] != '\t') {
					alt_allele = buffer[i];
				}
				if (count == 9 && buffer[i - 1] == '\t') {
					if (buffer[i] != '.') {
						//	std::cout<<"SET mother: "<<buffer[i-1]<<buffer[i]<<buffer[i+1]<<std::endl;
						parental = 2;
					}
				}
				if (count == 10 && buffer[i - 1] == '\t') {
					if (buffer[i] != '.') {
						if (parental == 0) {
							//	std::cout<<"SET father "<<buffer[i-1]<<buffer[i]<<buffer[i+1]<<std::endl;
							parental = 1;
						} else {
							parental = 0;
						}
					}
					break;
				}

				if (buffer[i] == '\t') {
					count++;
				}
			}

			if (parental != -1 && strncmp(chr.c_str(), "17", 2) == 0) {
				for (size_t i = 0; i < snps.size(); i++) {
					if (snps[i].position == pos && strncmp(snps[i].chr.c_str(), chr.c_str(), chr.size()) == 0) {
						//found the snp in the offspring (hapcut2 output)
						if (snps[i].alt_allele == alt_allele) {
							snps[i].parental = parental;
						} else {
							std::cerr << "Warning alt allele is different: " << snps[i].alt_allele << " != " << alt_allele << std::endl;
						}
						break;
					}
				}
			}

			getline(myfile, buffer);
		}
		myfile.close();
	}

}

void update_parents_gatk(std::string gatk_output, std::vector<snp_str> &snps) {
	std::string buffer;
	std::ifstream myfile;
	myfile.open(gatk_output.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "GATK Parser: could not open file: " << gatk_output.c_str() << std::endl;
		//exit(0);
	} else {
		getline(myfile, buffer);
		while (!myfile.eof()) {
			int pos = 0;
			std::string chr;
			char alt_allele = ' ';
			short parental = -1; //0=na ; 1=father; 2=mother;
			int count = 0;
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 0 && buffer[i] != '\t') {
					chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					pos = atoi(&buffer[i]);
				}
				if (count == 4 && buffer[i] != '\t') {
					alt_allele = buffer[i];
				}
				if (count == 8 && buffer[i - 1] == '\t') {
					if (strncmp("P1orP2", &buffer[i], 6) == 0) {
						parental = 0;
					} else if (buffer[i + 1] == '2') {
						parental = 1;
					} else if (buffer[i + 1] == '1') {
						parental = 2;
					} else {
						std::cout << "Prob parsing: " << buffer << std::endl;
					}
					break;
				}

				if (buffer[i] == '\t') {
					count++;
				}
			}
			if (parental != -1 && strncmp(chr.c_str(), "17", 2) == 0) {
				for (size_t i = 0; i < snps.size(); i++) {
					if (snps[i].position == pos && strncmp(snps[i].chr.c_str(), chr.c_str(), chr.size()) == 0) {
						//found the snp in the offspring (hapcut2 output)
						if (snps[i].alt_allele == alt_allele) {
							snps[i].gatk = parental;
						} else {
							std::cerr << "Warning alt allele is different (GATK): " << snps[i].alt_allele << " != " << alt_allele << std::endl;
						}
						break;
					}
				}
			}

			getline(myfile, buffer);
		}
		myfile.close();
	}
}

int is_in(std::string chr, int pos, std::vector<snp_str> &snps) {
	for (size_t i = 0; i < snps.size() && snps[i].position <= pos; i++) {
		if (snps[i].position == pos) {
			return (int) i;
		}
	}
	return -1;
}
double parse_gt(std::string buffer) {
	//parse GT freq:
	//	 0/1:51:214:265:. 2nd, 4th number: Major alt alle, read depth
	int count = 0;
	double alt = 0;
	double cov = 0;
	//std::cout<<buffer[0]<<buffer[1]<<buffer[2]<<buffer[3]<<std::endl;
	for (size_t i = 0; i < buffer.size(); i++) {
		if (count == 1 && buffer[i - 1] == ':') {
			alt = atof(&buffer[i]);

		}
		if (count == 3 && buffer[i - 1] == ':') {
			cov = atof(&buffer[i]);
			//	std::cout<<"Alele: " << alt<<" "<<cov<<std::endl;
			return alt / cov;
		}
		if (buffer[i] == ':') {
			count++;
		}
	}
	return -1;
}
void update_snps_alles(std::string snp_output, std::vector<snp_str> &snps) {

	std::string buffer;
	std::ifstream myfile;
	myfile.open(snp_output.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SNP Parser: could not open file: " << snp_output.c_str() << std::endl;
		//exit(0);
	} else {
		getline(myfile, buffer);
		while (!myfile.eof()) {
			int count = 0;
			int pos = 0;
			int index = -1;
			std::string chr = "";

			for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {

				if (count == 0 && buffer[i] != '\t') {
					chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					pos = atoi(&buffer[i]);
					index = is_in(chr, pos, snps);

				}
				if (count == 6 && buffer[i] != '\t') {
					if (index != -1) {
						snps[index].qual += buffer[i];
					}
				}
				if (count == 9 && buffer[i - 1] == '\t') {
					if (index != -1) {
						snps[index].ratio = parse_gt(&buffer[i]);
					}
				}

				if (buffer[i] == '\t') {
					count++;
				}

			}
			getline(myfile, buffer);
		}

		myfile.close();
	}
}
void parental_phasing(std::string parents_vcf, std::string hapcut_output, std::string gatk_output, std::string snp_file, std::string output) {
	std::vector<snp_str> snps = parse_hapcut2(hapcut_output);

//parse with parental SNPs:
	update_parents_xatlas(parents_vcf, snps);

//GATK parsing:
	update_parents_gatk(gatk_output, snps);

//Add allel frequencies/ratios from propant
	update_snps_alles(snp_file, snps);

//print results:

	FILE * file;
	file = fopen(output.c_str(), "w");

	for (size_t i = 0; i < snps.size(); i++) {
		fprintf(file, "%s", snps[i].chr.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", snps[i].position);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", snps[i].position + 1);
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", snps[i].qual.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%c", snps[i].alt_allele);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", snps[i].ratio);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", snps[i].phase_block);
		fprintf(file, "%c", '_');
		if (snps[i].haplotype) {
			fprintf(file, "%i", 1);
		} else {
			fprintf(file, "%i", 2);
		}
		fprintf(file, "%c", '\t');
		switch (snps[i].parental) {
		case 0:
			fprintf(file, "%s", "NA");
			break;
		case 1:
			fprintf(file, "%s", "father");
			break;
		case 2:
			fprintf(file, "%s", "mother");
			break;
		default:
			std::cerr << "Warning undefined parental value. Check parser!" << std::endl;
			break;
		}
		fprintf(file, "%c", '\t');
		switch (snps[i].gatk) {
		case 0:
			fprintf(file, "%s", "NA");
			break;
		case 1:
			fprintf(file, "%s", "father");
			break;
		case 2:
			fprintf(file, "%s", "mother");
			break;
		default:
			std::cerr << "Warning undefined gatk value. Check parser!" << std::endl;
			break;
		}
		fprintf(file, "%c", '\n');
	}
}

