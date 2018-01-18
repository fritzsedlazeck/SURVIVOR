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

//check this out!

void overlap_snps(std::string svs_file, std::string snp_file, int max_dist, int min_svs, int allele, std::string output) {

	Parameter::Instance()->min_freq = (double) allele / 100; // 0.05;   //TODO parameter?
	//cout << "AF: " << Parameter::Instance()->min_freq << endl;

	std::map<std::string, std::string> genome = parse_vcf(svs_file);

	std::vector<strvcfentry> entries = parse_vcf(svs_file, min_svs);

	/*	for (size_t j = 0; j < entries.size(); j++) {
	 cout << entries[j].start.pos << "\t" << entries[j].stop.pos << "\t" << entries[j].type << endl;
	 }
	 exit(0);*/
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

	FILE *file;
	file = fopen(output.c_str(), "w");

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(snp_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SNP Parser: could not open file: " << snp_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size); //avoid header
	//fprintf(file, "%s", buffer);
//	fprintf(file, "%s", "\tSV\n");

	myfile.getline(buffer, buffer_size);

	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			breakpoint_str snp;
			double pval = 1000;
			int start = 0;
			int stop = 0;
			int pop = 0;
			bool flag = true;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 11 && buffer[i] != '\t') {
					if (buffer[i] == ';') {
						flag = false;
					}
					if (flag) {
						snp.chr += buffer[i];
					}
					//std::cout << "to find" << snp.chr << " ";
				}
				if (count == 12 && buffer[i - 1] == '\t') {
					start = atoi(&buffer[i]);
					stop = start;
					//std::cout << snp.chr <<" " <<start << std::endl;
				}
				if (count == 30 && buffer[i - 1] == '\t') {
					pval = atof(&buffer[i]);
					//if (start == 119380704) {
					//	cout << pval<<" "<<buffer[i]<<buffer[i+1]<<buffer[i+2]<< endl;
					//}
				}
				if (count == 37 && buffer[i - 1] == '\t') {
					pop = atoi(&buffer[i]);

					//	cout<<pop<<" "<<pval<<endl;
					//	exit(0);
					break;
				}
				/*if (count == 2 && buffer[i - 1] == '\t') {
				 stop = atoi(&buffer[i]);
				 }*/
				if (buffer[i] == '\t') {
					count++;
					flag = true;
				}
			}
			if (start == 119380704) {
				cout << "pop: " << pop << " " << pval << endl;
			}
			if (pop >= 1000 && pval > 0) {
				std::string SV = "NA";
				for (size_t i = start; i <= stop; i++) {
					if (genome.find(snp.chr) != genome.end() && genome[snp.chr].size() > i) {
						if (genome[snp.chr][i] == 'T') {
							SV = "TRA";
							break;
						} else if (genome[snp.chr][i] == 'I') {
							SV = "INS";
							break;
						} else if (genome[snp.chr][i] == 'U') {
							SV = "UNK";
							break;
						} else if (genome[snp.chr][i] == 'D') {
							SV = "DEL";
							break;
						} else if (genome[snp.chr][i] == 'P') {
							SV = "DUP";
							break;
						} else if (genome[snp.chr][i] == 'C') {
							SV = "CNV";
							break;
						} else if (genome[snp.chr][i] == 'V') {
							SV = "INV";
							break;
						} else {
							SV = "NA";
						}
					}
				}

				/*fprintf(file, "%s", snp.chr.c_str());
				 fprintf(file, "%c", '\t');
				 fprintf(file, "%i", start);
				 fprintf(file, "%c", '\t');*/
				fprintf(file, "%e", pval);
				//	fprintf(file, "%s", std::string(buffer).c_str());
				fprintf(file, "%c", '\t');
				/*	if (genome.find(snp.chr) != genome.end() && genome[snp.chr].size() > snp.position) {
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
				 } else {
				 std::cout << "error: " << snp.chr << " " << snp.position << std::endl;
				 fprintf(file, "%s", "ERR");
				 }*/
				fprintf(file, "%s", SV.c_str());
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
	size_t id = 0;
	breakpoint_str tmp;
	int pos = atoi(&buffer[i]);
	while (buffer[i] != '\t') {
		if (buffer[i - 1] == ';') {
			if (id < snps.size()) {
				snps[id].position = pos;					//
			}
			pos = atoi(&buffer[i]);
			id++;
		}
		i++;
	}
	if (id < snps.size()) {
		snps[id].position = pos;

	} else {
		cout << "Strange" << endl;
	}
}

void parse_pos(std::string buffer, int& pos, std::string& chr) {
	size_t i = 0;
	int count = 0;

	while (buffer[i] != '\t') {
		if (count == 0 && buffer[i] != ':') {
			chr += buffer[i];
		}
		if (count == 1 && buffer[i - 1] == ':') {
			pos = atoi(&buffer[i]);
			break;
		}
		if (buffer[i] == ':') {
			count++;
		}
		i++;
	}
}

std::vector<strvcfentry> parse_random(std::string svs_file) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];

	std::ifstream myfile;
	myfile.open(svs_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SNP Parser: could not open file: " << svs_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	std::vector<strvcfentry> entries;
	while (!myfile.eof()) {
		strvcfentry tmp;
		int count = 0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
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
			if (count == 3 && buffer[i - 1] == '\t') {
				tmp.type = atoi(&buffer[i]);
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		entries.push_back(tmp);
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	return entries;

}
void overlap_snps_gwas(std::string svs_file, std::string random_SV, int max_dist, int min_svs, std::string output) {

	Parameter::Instance()->min_freq = 0.05;   //TODO parameter?
	std::map<std::string, std::string> genome = parse_vcf(svs_file);
	std::vector<strvcfentry> entries;

	if (strcmp(random_SV.c_str(), "N") == 0) {
		entries = parse_vcf(svs_file, min_svs);
	} else {
		entries = parse_random(random_SV);
	}

	std::cout << "merging entries: " << entries.size() << std::endl;
	for (size_t j = 0; j < entries.size(); j++) {
		//	cout << "Prev\t" << entries[j].start.chr << "\t" << entries[j].start.pos << "\t" << entries[j].stop.pos << "\t" << entries[j].type << endl;

		int start = std::max(0, (int) entries[j].start.pos - max_dist);
		int stop = std::min((int) genome[entries[j].start.chr].size(), (int) entries[j].stop.pos + max_dist);

		if (entries[j].type == 0 || (entries[j].type == 1 || entries[j].type == 6)) { //entire region
			char sv = 'D';   //del
			if (entries[j].type == 1) {
				sv = 'P'; //DUP
			} else if (entries[j].type == 6) {
				sv = 'C'; //CNV
			}
			//	cout<<"cover: "<<entries[j].start.chr<<" "<<start<<" "<<stop-start<<endl;
			for (int pos = start; pos < stop; pos++) {
				genome[entries[j].start.chr][pos] = sv;
			}
		} else { // over breakpoints
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

	cout << "Start parsing" << endl;
	FILE *file;
	file = fopen(output.c_str(), "w");

	int num = 0;
	while (!cin.eof()) {
		std::string buffer;
		getline(cin, buffer);
		if (!cin.fail()) {
			if (num == 0) {
				fprintf(file, "%s", buffer.c_str());
				fprintf(file, "%c", '\t');
				fprintf(file, "%s", "SVs");
				fprintf(file, "%c", '\n');
			} else {
				//	int count = 0;
				std::string chr = "";
				int pos = 0;
				parse_pos(buffer, pos, chr);
				/*	for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {

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
				 }*/
				//cout<<chr<<" "<<pos<<endl;
				fprintf(file, "%s", buffer.c_str());
				fprintf(file, "%c", '\t');
				std::string SV = "NA";
				if (genome.find(chr) != genome.end() && genome[chr].size() > pos) {
					if (genome[chr][pos] == 'T') {
						SV = "TRA";
					} else if (genome[chr][pos] == 'I') {
						SV = "INS";
					} else if (genome[chr][pos] == 'U') {
						SV = "UNK";
					} else if (genome[chr][pos] == 'D') {
						SV = "DEL";
					} else if (genome[chr][pos] == 'P') {
						SV = "DUP";
					} else if (genome[chr][pos] == 'C') {
						SV = "CNV";
					} else if (genome[chr][pos] == 'V') {
						SV = "INV";
					}
				}
				fprintf(file, "%s", SV.c_str());
				fprintf(file, "%c", '\n');
			}
			num++;
		} else {
			break;
		}
	}

	/*	exit(0);
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
	 fprintf(file, "%s", "SVs");
	 fprintf(file, "%c", '\t');
	 fprintf(file, "%s", "Occurs");
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
	 fprintf(file, "%c", '\t');
	 char nuc = genome[snp.chr][snp.position];
	 if (genome[snp.chr][snp.position] == 'T' || genome[snp.chr][snp.position] == 'S') {
	 fprintf(file, "%s", "TRA");
	 genome[snp.chr][snp.position] = 'S';
	 } else if (genome[snp.chr][snp.position] == 'I' || genome[snp.chr][snp.position] == 'J') {
	 fprintf(file, "%s", "INS");
	 genome[snp.chr][snp.position] = 'J';
	 } else if (genome[snp.chr][snp.position] == 'U' || genome[snp.chr][snp.position] == 'W') {
	 fprintf(file, "%s", "UNK");
	 genome[snp.chr][snp.position] = 'W';
	 } else if (genome[snp.chr][snp.position] == 'D' || genome[snp.chr][snp.position] == 'E') {
	 fprintf(file, "%s", "DEL");
	 genome[snp.chr][snp.position] = 'E';
	 } else if (genome[snp.chr][snp.position] == 'P' || genome[snp.chr][snp.position] == 'O') {
	 fprintf(file, "%s", "DUP");
	 genome[snp.chr][snp.position] = 'O';
	 } else if (genome[snp.chr][snp.position] == 'C' || genome[snp.chr][snp.position] == 'B') {
	 fprintf(file, "%s", "CNV");
	 genome[snp.chr][snp.position] = 'B';
	 } else if (genome[snp.chr][snp.position] == 'V' || genome[snp.chr][snp.position] == 'X') {
	 fprintf(file, "%s", "INV");
	 genome[snp.chr][snp.position] = 'X';
	 } else if (genome[snp.chr][snp.position] == 'N' || genome[snp.chr][snp.position] == 'M') {
	 fprintf(file, "%s", "NA");
	 genome[snp.chr][snp.position] = 'M';
	 }
	 fprintf(file, "%c", '\t');
	 if (((nuc == 'T' || nuc == 'I') || (nuc == 'U' || nuc == 'D')) || ((nuc == 'P' || nuc == 'C') || (nuc == 'V' || nuc == 'N'))) {
	 fprintf(file, "%s", "UNI");
	 } else {
	 fprintf(file, "%s", "REP");
	 }
	 fprintf(file, "%c", '\n');
	 }
	 }
	 myfile.getline(buffer, buffer_size);
	 }
	 myfile.close();*/
	fclose(file);
}

void overlap_snpsGWASDB(std::string svs_file, std::string snp_file, int max_dist, int min_svs, int allele, std::string output) {

	Parameter::Instance()->min_freq = (double) allele / 100;   //TODO parameter?
	cout << "Freq: " << Parameter::Instance()->min_freq << endl;
	std::map<std::string, std::string> genome = parse_vcf(svs_file);

	std::vector<strvcfentry> entries = parse_vcf(svs_file, min_svs);
	/*srand(time(NULL));
	 for (size_t i = 0; i < entries.size(); i++) {
	 int chr = rand() % genome.size();
	 int num = 0;
	 for (std::map<std::string, std::string>::iterator j = genome.begin(); j != genome.end(); j++) {
	 if (num == chr) {
	 cout << "Prev\t" << entries[i].start.chr << "\t" << entries[i].start.pos << "\t" << entries[i].stop.pos << "\t" << entries[i].type << endl;
	 entries[i].start.chr = (*j).first;
	 entries[i].stop.chr = (*j).first;
	 int len = entries[i].stop.pos - entries[i].start.pos;
	 entries[i].start.pos = rand() % ((*j).second.size() - len);
	 entries[i].stop.pos = entries[i].start.pos + len;
	 cout << "New\t" << entries[i].start.chr << "\t" << entries[i].start.pos << "\t" << entries[i].stop.pos << "\t" << entries[i].type << endl;
	 break;
	 }
	 num++;
	 }
	 }
	 exit(0);*/

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
	//fprintf(file, "%s", buffer);
//	fprintf(file, "%c", '\t');
	fprintf(file, "%s", "CHR\tPOS\tPval\tSVs");
	fprintf(file, "%c", '\n');
	myfile.getline(buffer, buffer_size);

	while (!myfile.eof()) {
		//cout << "parse" << endl;
		int pop = 0;
		int count = 0;
		std::vector<breakpoint_str> snps;

		snps = parse_chrs(buffer);

		double pvalue = 10000;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {

			if (count == 1 && buffer[i - 1] == '\t') {
				parse_pos(&buffer[i], snps);
			}
			if (count == 4 && buffer[i - 1] == '\t') {
				pvalue = atof(&buffer[i]);
			}
			if (count == 6 && buffer[i] != 'N' && buffer[i - 1] == '\t') {
				pop = atoi(&buffer[i]);
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}

		if (pop > 1000) {
			//cout<<snps[0].chr<<" "<< snps[0].position<<endl;
			//cout << "print: " << snps.size() << endl;
			for (size_t i = 0; i < snps.size(); i++) {
				if (genome.find(snps[i].chr) != genome.end() && genome[snps[i].chr].size() > snps[i].position) {

					breakpoint_str snp = snps[i];
					fprintf(file, "%s", snps[i].chr.c_str());
					fprintf(file, "%c", '\t');
					fprintf(file, "%i", snps[i].position);
					fprintf(file, "%c", '\t');
					fprintf(file, "%e", pvalue);
					fprintf(file, "%c", '\t');
					std::string SV = "NA";
					if (genome[snp.chr][snp.position] == 'T') {
						SV = "TRA";
					} else if (genome[snp.chr][snp.position] == 'I') {
						SV = "INS";
					} else if (genome[snp.chr][snp.position] == 'U') {
						SV = "UNK";
					} else if (genome[snp.chr][snp.position] == 'D') {
						SV = "DEL";
					} else if (genome[snp.chr][snp.position] == 'P') {
						SV = "DUP";
					} else if (genome[snp.chr][snp.position] == 'C') {
						SV = "CNV";
					} else if (genome[snp.chr][snp.position] == 'V') {
						SV = "INV";
					}
					fprintf(file, "%s", SV.c_str());
					fprintf(file, "%c", '\n');
				}
			}
		}
		//	cout << "fin" << endl;
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
}

std::map<std::string, std::string> parse_fasta(std::string genome_file) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(genome_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SNP Parser: could not open file: " << genome_file.c_str() << std::endl;
		exit(0);
	}

	std::map<std::string, std::string> genome;
	myfile.getline(buffer, buffer_size);

	std::string name="";
	std::string sequence="";
	while (!myfile.eof()) {
		if (buffer[0] == '>') {
			if (!sequence.empty()) {
				genome[name] = sequence;
			}
			name.clear();
			sequence.clear();
			//skip >
			for (size_t i = 1; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n' && buffer[i] != ' '; i++) {
				name += buffer[i];
			}
		} else {
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n' && buffer[i] != ' '; i++) {
				sequence += buffer[i];
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	if (!sequence.empty()) {
		genome[name] = sequence;
	}
	myfile.close();
	return genome;
}

void generate_random_regions(std::string genome_file, std::string svs_vcf, int min_svs, std::string output) {
	std::map<std::string, std::string> genome = parse_fasta(genome_file);
	cout << "CHR detected: " << genome.size() << endl;
	Parameter::Instance()->min_freq = 0.05;
	std::vector<strvcfentry> entries = parse_vcf(svs_vcf, min_svs);
	cout << "Parsed SVS: " << entries.size() << endl;
	srand(time(NULL));
	FILE *file;
	file = fopen(output.c_str(), "w");

	for (size_t i = 0; i < entries.size(); i++) {

		bool valid_region = false;
		while (!valid_region) {
			int chr = rand() % genome.size();
			//cout<<"CHR: "<<chr;
			int num = 0;
			for (std::map<std::string, std::string>::iterator j = genome.begin(); j != genome.end(); j++) {
				if (num == chr) {
					entries[i].start.chr = (*j).first;
					entries[i].stop.chr = (*j).first;
				//	cout<<" "<<chr<<" "<<entries[i].start.chr <<endl;
					int len = entries[i].stop.pos - entries[i].start.pos;
					entries[i].stop.pos = (*j).second.size();
					entries[i].start.pos = rand() % ((int) (*j).second.size() - len*2);
					entries[i].stop.pos = entries[i].start.pos + len;

					double n_count = 0;
					for (size_t pos = entries[i].start.pos; pos < entries[i].stop.pos; pos++) {
						if (toupper(genome[entries[i].start.chr][pos]) == 'N') {
							n_count++;
						}
					}
				//	cout<<"Valid: "<<n_count<<endl;
					valid_region = (bool) (n_count / (double) len < 0.05);
					break;
				}
				num++;

			}

		}
		for (size_t pos = entries[i].start.pos; pos < entries[i].stop.pos; pos++) {
			genome[entries[i].start.chr][pos]='N'; //masking region such that its not chosen again.
		}

		//print new coords:
		fprintf(file, "%s", entries[i].start.chr.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", entries[i].start.pos);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", entries[i].stop.pos);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", entries[i].type);
		fprintf(file, "%c", '\n');
	//	cout<<"Done: "<<i<<endl;

	}
	fclose(file);

}

