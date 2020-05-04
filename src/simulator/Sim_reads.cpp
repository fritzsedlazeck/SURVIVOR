/*
 * Nanopore_sim.cpp
 *
 *  Created on: May 30, 2017
 *      Author: sedlazec
 */

#include "Sim_reads.h"
std::map<std::string, std::string> parse_genome(std::string genome_file, int min_length) {
	size_t buffer_size;
	char*buffer;
	ifstream myfile;

	myfile.open(genome_file.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Fasta Parser: could not open file: " << genome_file.c_str() << endl;
		exit(0);
	}

	buffer_size = 20000;
	buffer = new char[buffer_size];

	myfile.getline(buffer, buffer_size);
	string seq = "";
	string name;

//	std::default_random_engine generator;
//	std::poisson_distribution<int> distribution(2);

	std::map<std::string, std::string> genome;
	while (!myfile.eof()) {
		if (buffer[0] == '>') {
			if ((int) seq.size() > min_length) {

				//	int border=distribution(generator);
				//cout<<"max copies: "<<border<<std::endl;
				//for (size_t i = 0; i < border; i++) {
				stringstream ss;
				ss << name;
				//ss << "_";
				//ss << 0;
				//	ss << i;
				//	std::cout<<ss.str()<<std::endl;
				genome[ss.str()] = seq;
				//}

			}
			name.clear();
			seq.clear();
			for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != ' ' && buffer[i] != '\0'; i++) {
				name += toupper(buffer[i]);
			}
		} else {
			for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
				seq += toupper(buffer[i]);
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	if ((int) seq.size() > min_length) {
		stringstream ss;
		ss << name;
		//	ss << i;
		//	std::cout<<ss.str()<<std::endl;
		genome[ss.str()] = seq;
	}
	name.clear();
	seq.clear();
	myfile.close();
	return genome;
}
char complementbp(char old) {

	switch (old) {
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
	default:
		return 'N';
		break;
	}
}

char new_nuc(char old) {
	std::string nucs = "ACGT";
	int index = rand() % 4;
	switch (old) {
	case 'A':
		while (index != 0) {
			index = rand() % 4;
		}
		return nucs[index];
		break;
	case 'C':
		while (index != 1) {
			index = rand() % 4;
		}
		return nucs[index];
		break;
	case 'G':
		while (index != 2) {
			index = rand() % 4;
		}
		return nucs[index];
		break;
	case 'T':
		while (index != 3) {
			index = rand() % 4;
		}
		return nucs[index];
		break;

	default:
		return nucs[index];
		break;
	}
}

double rand_normal(double mean, double stddev) {					//Box muller method
	static double n2 = 0.0;
	static int n2_cached = 0;
	if (!n2_cached) {
		double x, y, r;
		do {
			x = 2.0 * rand() / RAND_MAX - 1;
			y = 2.0 * rand() / RAND_MAX - 1;

			r = x * x + y * y;
		} while (r == 0.0 || r > 1.0);
		{
			double d = sqrt(-2.0 * log(r) / r);
			double n1 = x * d;
			n2 = y * d;
			double result = n1 * stddev + mean;
			n2_cached = 1;
			return result;
		}
	} else {
		n2_cached = 0;
		return n2 * stddev + mean;
	}
}

std::vector<read_position> parse_error_profile(std::string error_profile_file) {
	char*buffer;
	ifstream myfile;

	myfile.open(error_profile_file.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Fasta Parser: could not open file: " << error_profile_file.c_str() << endl;
		exit(0);
	}

	size_t buffer_size = 2000;
	buffer = new char[buffer_size];

	myfile.getline(buffer, buffer_size);					//avoid header
	myfile.getline(buffer, buffer_size);

	std::vector<read_position> error_profile;
	while (!myfile.eof()) {
		int count = 0;
		read_position tmp;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 1 && buffer[i - 1] == '\t') {
				tmp.total = atof(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				tmp.match = atof(&buffer[i]);
			}
			if (count == 3 && buffer[i - 1] == '\t') {
				tmp.mismatch = atof(&buffer[i]);
			}
			if (count == 4 && buffer[i - 1] == '\t') {
				tmp.ins = atof(&buffer[i]);
			}
			if (count == 5 && buffer[i - 1] == '\t') {
				tmp.del = atof(&buffer[i]);
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		error_profile.push_back(tmp);
		myfile.getline(buffer, buffer_size);
	}

	return error_profile;
}

void simulate_reads(std::string genome_file, std::string error_profile_file, int coverage, std::string output) {
	srand(time(NULL));

	std::vector<read_position> error_profile = parse_error_profile(error_profile_file);
	int avg_readlen = 0;
	///int max_size = error_profile.size() - 1;
	for (size_t i = 0; i < error_profile.size() && error_profile[i].total < 0.6; i++) {
		avg_readlen++;
	}

	std::map<std::string, std::string> genome = parse_genome(genome_file, avg_readlen * 3);
	std::cout << "\tParsing done: " << genome.size() << " chrs " << std::endl;
	long genome_size = 0;
	int hits = 0;
	for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
		hits++;
		genome_size += (*i).second.size();
	}

	std::cout << "\tAVG read length: " << avg_readlen << std::endl;
	int num_reads = ((double) genome_size / (double) avg_readlen) * (double) coverage;
	std::cout << "\tNum of reads: " << num_reads << std::endl;
	FILE *file;
	file = fopen(output.c_str(), "w");
	double prev = 0;
//	cout << "\t\tReads simulated: " << prev << "%" << std::endl;
	for (int i = 0; i < num_reads; i++) { //start to simulate reads:
		double bp = (rand() % 1000000);
	//	if (bp / 1000000 < 0.0000001) {
	//		cout << "BP " << bp << " " << bp / 1000000 << endl;
	//	}
		bp = bp / 1000000;
		size_t size = 0;
		while (size < error_profile.size()) { //1: Pick a read size based on the profile
			if (size > 10 && bp < error_profile[size].total) {
				break;
			}
			size++;
		}

		//	cout << size << ": " << bp << endl;
		//std::cout<<"Read init: "<<size<<std::endl;
		size_t len = 0;
		std::string chr;
		while (size > len) {  //2: Pick a chromosome:
			int pos = rand() % (int) (genome.size());
			chr = (*genome.begin()).first; //check that again...
			for (std::map<std::string, std::string>::iterator j = genome.begin(); j != genome.end() && pos >= 0; j++) {
				if (pos == 0) {
					chr = (*j).first;
				}
				pos--;
			}
			len = genome[chr].size();
		}
		//cout << "choose subregion" << chr << endl;
		std::string read = "";
		double num_N = 0;
		size_t start_pos = 0;

		//Obtaining the sequence (avoid large regions of N's)
		int attempts = 10;
		std::string final_read;
		do { //2: Pick a position on chromosome and obtain sequence avoiding too many N's:
			final_read = "";
			start_pos = (size_t) rand() % (genome[chr].size() - size + 2);

			//new attempt to make it faster:
			//	cout << " Start2: " << chr << " " << start_pos  << endl;
			read = genome[chr].substr(start_pos, size);
			if (read.empty()) {
				cerr << "ERROR! " << endl;
			}

			for (size_t j = 0; j < read.size(); j++) {
				if (read[j] == 'N' || read[j] == 'n') {
					num_N++;
					if (num_N / (double) read.size() > 0.1) {
						break;
					}
					read[j] = new_nuc('N');
				}
				double bp = (rand() % 1000000); //bp probability
				bp = bp / 1000000;
				if (bp < error_profile[j].match) {
					final_read += read[j];
				} else if (bp < error_profile[j].mismatch + error_profile[j].match) {
					final_read += new_nuc(read[j]);
				} else if (bp < error_profile[j].ins + error_profile[j].match + error_profile[j].mismatch) {
					final_read += read[j];
					final_read += new_nuc('N');
				}
			}
			attempts--;
			//	cout << "Ns: " << num_N <<" size: "<<(double) read.size() << endl;
		} while (num_N / (double) read.size() > 0.1 && attempts != 0);
		if (attempts != 0) {
			bool flag = true;
			if (rand() % 100 < 51) {
				flag = false;
				//reverse read.
				std::string new_read;
				for (std::string::reverse_iterator ri = final_read.rbegin(); ri != final_read.rend(); ri++) {
					new_read += complementbp((*ri));
				}
				final_read = new_read;
			}

			std::stringstream name;
			name << chr;
			name << "_";
			name << start_pos;
			if (flag) {
				name << "_+";
			} else {
				name << "_-";
			}

			fprintf(file, "%s", name.str().c_str());
			fprintf(file, "%c", '\n');
			fprintf(file, "%s", final_read.c_str());
			fprintf(file, "%c", '\n');

			if (i % 10000 == 0 && prev < (i * 100) / num_reads) {
				prev = (i * 100) / num_reads;
				std::cout << "\t\tReads simulated: " << prev << "%" << '\r' << std::flush;
				//cout << "\t\tReads simulated: " << prev << "%" << std::endl;
			}
		}
	}

	fclose(file);
}
