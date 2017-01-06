/*
 * Pac_Simulator.cpp
 *
 *  Created on: Feb 1, 2016
 *      Author: fsedlaze
 */

#include "Pac_Simulator.h"
char ins() {
	switch (rand() % 4) {
	case 0:
		return 'A';
		break;
	case 1:
		return 'C';
		break;
	case 2:
		return 'G';
		break;
	case 3:
		return 'T';
		break;
	}
	return 'N'; //should not happen
}
void print_sam(FILE*& file, std::string name, int pos, std::string new_seq, std::string cigar, bool strand) {
	fprintf(file, "%s", name.c_str());
	fprintf(file, "%c", '_');
	fprintf(file, "%i", pos);
	if (strand) {
		fprintf(file, "%s", "\t0\t");
	} else {
		fprintf(file, "%s", "\t16\t");
	}
	fprintf(file, "%s", name.c_str());
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", pos + 1);
	fprintf(file, "%s", "\t255\t");
	fprintf(file, "%s", cigar.c_str());
	fprintf(file, "%s", "\t*\t0\t0\t");
	fprintf(file, "%s", new_seq.c_str());
	fprintf(file, "%c", '\t');
	for (size_t pos = 0; pos < new_seq.size(); pos++) {
		fprintf(file, "%c", 'H');
	}
	fprintf(file, "%s", "\tNM:i:4");
	fprintf(file, "%c", '\n');
}

char comp(char base) {
	switch (base) {
	case 'A':
		return 'T';
		break;
	case 'C':
		return 'G';
		break;
	case 'G':
		return 'C';
		break;
	case 'T':
		return 'A';
		break;
	}
	return base;
}
void rev_comp(std::string & read) {
	std::string new_read;
	for (std::string::reverse_iterator i = read.rbegin(); i != read.rend(); i++) {
		new_read += comp((*i));
	}
	read.clear();
	read = new_read;
}
void simulate_reads(std::string name, std::string seq, FILE*& file, FILE*& sam, FILE*& file2) {
	size_t i = 0;
	int len = 20000 + rand() % 1000;
	while (i < seq.size()) {
		if (i + len < seq.size()) {
			fprintf(file, "%c", '@');
			fprintf(file, "%s", name.c_str());
			fprintf(file, "%c", '_');
			fprintf(file, "%i", (int) i);
			fprintf(file, "%c", '\n');

			fprintf(file2, "%c", '>');
			fprintf(file2, "%s", name.c_str());
			fprintf(file2, "%c", '_');
			fprintf(file2, "%i", (int) i);
			fprintf(file2, "%c", '\n');

			std::string read = seq.substr(i, len * 2);
			bool strand = true;
			if (rand() % 100 < 50) {
				strand = false;
				rev_comp(read);
			}
			std::string new_seq="";
			int tmp = 1;
			std::stringstream ss;

			//14M2D3M2D3M     *       0       0
			//AGCTTTTCATTCTA--CGC--A
			//	14M1D2M1D3M     *       0       0
			//AGCTTTTCATTCTA CG CA
			char mod = ' ';
			size_t pos = 0;
			while (new_seq.size() < len && pos < read.size()) {
				if (rand() % 100 < 15 && (pos > 0 && pos < read.size() - 1)) { //why 4??
					if (rand() % 100 < 40) {
						if (mod != 'D' && mod != ' ') {
							ss << tmp;
							ss << mod;
							tmp = 0;
						}
						tmp++;
						mod = 'D';
						pos++;
						//deletion
					} else {
						if (mod != 'I' && mod != ' ') {
							ss << tmp;
							ss << mod;
							tmp = 0;
						}
						//insertion
						tmp++;
						new_seq += ins();
						mod = 'I';
					}
				} else {
					if (mod != 'M' && mod != ' ') {
						ss << tmp;
						ss << mod;
						tmp = 0;
					}
					mod = 'M';
					tmp++;
					new_seq += read[pos];
					pos++;
				}
			}
			if (tmp - 1 > 0) {
				ss << tmp - 1;
				ss << mod;
			}
			for (size_t pos = 0; pos < new_seq.size(); pos++) {
				fprintf(file, "%c", new_seq[pos]);
				fprintf(file2, "%c", new_seq[pos]);
			}
			fprintf(file, "%s", "\n+\n");
			fprintf(file2, "%s", "\n");
			for (size_t pos = 0; pos < new_seq.size(); pos++) {
				fprintf(file, "%c", 'H');
			}
			fprintf(file, "%c", '\n');

			print_sam(sam, name, i, new_seq, ss.str(), strand);
		}
		i += 286;
	}
}

void simulate_pac(std::string genome, std::string output) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(genome.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SAM Parser: could not open file: " << genome.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	std::string name;
	std::string seq;
	std::string out = output;
	out += ".fq";
	FILE * file;
	file = fopen(out.c_str(), "w");
	FILE * file2;
	out = output;
	out += ".fa";
	file2 = fopen(out.c_str(), "w");
	FILE * sam;
	out = output;
	out += ".sam";
	sam = fopen(out.c_str(), "w");

	srand(time(NULL));
	while (!myfile.eof()) {
		if (buffer[0] == '>') {
			if (!seq.empty()) {
				simulate_reads(name, seq, file, sam, file2);
				name.clear();
				seq.clear();
			}
			for (size_t i = 1; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n' && buffer[i] != ' '; i++) {
				name += buffer[i];
			}
		} else {
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n' && buffer[i] != ' '; i++) {
				seq += buffer[i];
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	if (!seq.empty()) {
		simulate_reads(name, seq, file, sam, file2);
	}
	myfile.close();
	fclose(file);
	fclose(sam);

}
