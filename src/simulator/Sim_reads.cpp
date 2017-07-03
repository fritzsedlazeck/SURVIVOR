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
			if (seq.size() > min_length) {

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
			for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != ' ' &&buffer[i] != '\0'; i++) {
				name += toupper(buffer[i]);
			}
		} else {
			for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
				seq += toupper(buffer[i]);
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	if (seq.size() > min_length) {
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

	int min_length = 500;					//think about!
	std::map<std::string, std::string> genome = parse_genome(genome_file, min_length);
	std::cout << "\tParsing done: " << genome.size() << " chrs " << std::endl;
	long genome_size = 0;
	int hits = 0;
	for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
	//	size_t found = (*i).first.find("_0");
	//	if (found != std::string::npos) {
			hits++;
			genome_size += (*i).second.size();
	//	}
	}
	std::vector<read_position> error_profile = parse_error_profile(error_profile_file);
	int avg_readlen = 0;
	///int max_size = error_profile.size() - 1;
	for (size_t i = 0; i < error_profile.size() && error_profile[i].total < 0.51; i++) {
		avg_readlen++;
	}

	std::cout << "\tAVG read length: " << avg_readlen << std::endl;
	int num_reads = ((double) genome_size / (double) avg_readlen) * (double) coverage;
	std::cout << "\tNum of reads: " << num_reads << std::endl;
	FILE *file;
	file = fopen(output.c_str(), "w");

	for (size_t i = 0; i < num_reads; i++) {
		double bp = (rand() % 1000000);
		bp=bp/1000000;
		int size = 0;
		while(size<error_profile.size()){
			if(bp < error_profile[size].total){
				break;
			}
			size++;
		}

		//cout<<size<<": "<<bp<<endl;
		//std::cout<<"Read init: "<<size<<std::endl;
		int pos = rand() % (int) (genome.size());
	//	cout<<"chr: "<<pos;
		//std::cout<<" pos: "<<pos<<" "<<genome.size()<<std::endl;
		std::string chr = (*genome.begin()).first; //check that again...
		for (std::map<std::string, std::string>::iterator j = genome.begin(); j != genome.end() && pos >= 0; j++) {
			chr = (*j).first;
			pos--;
		}
		//std::cout<<"chr: "<<chr<<std::endl;
		std::string read = "N";
		double num_N = 1;
		int start_pos = 0;
		while (num_N / (double) read.size() > 0.1) { //optaining the sequence (avoid large regions of N's)
			read = "";
			num_N = 0;
			start_pos = rand() % (int) genome[chr].size();
			while ((start_pos + size) >= genome[chr].size()) { // check such that we dont always just get the ends.
				start_pos = rand() % (int) genome[chr].size();
			}

			for (size_t j = 0; j < size && j + (size_t) start_pos < genome[chr].size(); j++) {
				read += genome[chr][start_pos + j];
				if (genome[chr][start_pos + j] == 'N' || genome[chr][start_pos + j] == 'n' ) {
					num_N++;
				}
			}
		}

		size_t t = 0;
		std::string final_read = "";
		while (t < read.size()) {

			double bp = (rand() % 1000000); //bp probability
			bp=bp/1000000;

			//if (bp < error_profile[t].total) { //break the read!
			//	cout<<t<<": "<<bp<<endl;
			//	break;
			//}
			if (bp < error_profile[t].match) {
				final_read += read[t];
			} else if (bp - error_profile[t].match < error_profile[t].mismatch) {
				final_read += new_nuc(read[t]);
			} else if (bp - (error_profile[t].match + error_profile[t].mismatch) < error_profile[t].ins) {
				final_read += new_nuc('N');
			} //else is a deletion!
			t++;
		}

		std::stringstream name;
		name << chr;
		name << "_";
		name << start_pos;
		//	std::cout<<"Read: "<<read.size()<<std::endl;
		if (read.size() > min_length) {
			fprintf(file, "%s", name.str().c_str());
			fprintf(file, "%c", '\n');
			fprintf(file, "%s", final_read.c_str());
			fprintf(file, "%c", '\n');
		}
		if (i % 1000 == 0) {
			cout << "\t\tReads simulated: " << (i*100)/num_reads<<"%" << std::endl;
		}
	}

	fclose(file);
}
