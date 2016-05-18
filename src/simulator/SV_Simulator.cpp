/*
 * SV_Simulator.cpp
 *
 *  Created on: Jan 30, 2016
 *      Author: fsedlaze
 */

#include "SV_Simulator.h"

int parse_value(char* buffer, size_t buffer_size) { //required for parameter!
	int count = 0;
	for (size_t i = 1; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
		if (count == 1) {
			return atoi(&buffer[i]);
		}
		if (buffer[i] == ' ') {
			count++;
		}
	}
	return -1;
}

void simulate(std::map<std::string, std::string> genome, std::vector<struct_var> svs, std::string output_prefix) {

	std::cout << "apply SV" << std::endl;

}

parameter parse_param(std::string filename) {
	parameter tmp;
	size_t buffer_size = 200000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.dup_max = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.dup_min = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.dup_num = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);

	tmp.indel_min = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.indel_max = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.indel_num = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);

	tmp.translocations_min = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.translocations_max = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.translocations_num = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);

	tmp.inv_min = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.inv_max = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.inv_num = parse_value(buffer, buffer_size);
	myfile.close();
	return tmp;
}

std::map<std::string, std::string> read_fasta(std::string ref_file, int min_length) {
	size_t buffer_size = 200000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(ref_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: " << ref_file.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	std::map<std::string, std::string> genome;
	std::string seq;
	std::string name;
	while (!myfile.eof()) {
		if (buffer[0] == '>') {
			if (seq.size() > min_length) {
				genome[name] = seq;
			}
			name.clear();
			seq.clear();

			for (size_t i = 1; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0' && buffer[i] != ' '; i++) {
				name += (buffer[i]);
			}
		} else {
			for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
				seq += toupper(buffer[i]);
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
		seq += toupper(buffer[i]);
	}
	myfile.close();
	if (seq.size() > min_length) { //TODO do I really need this?
		genome[name] = seq;
	}
	seq.clear();
	name.clear();
	std::cout << "# Chrs passed size threshold:" << genome.size() << std::endl;
	return genome;
}

void sort_svs(std::vector<struct_var> svs) {

	std::map<std::string, std::vector<struct_var> > svs_tmp;
	for (size_t i = 0; i < svs.size(); i++) {
		svs_tmp[svs[i].pos.chr].push_back(svs[i]); //sort by chr:
	}

	for (std::map<std::string, std::vector<struct_var> >::iterator i = svs_tmp.begin(); i != svs_tmp.end(); i++) {
		std::vector<struct_var> tmp;
		if (!(*i).second.empty()) {
			tmp.push_back((*i).second[0]);
			for (size_t j = 0; j < (*i).second.size(); i++) {
				std::cout<<(*i).second[j].pos.chr.c_str()<<" "<<(*i).second[j].pos.start<<std::endl;
				size_t t=0;
				while(tmp[t].pos.start <(*i).second[j].pos.start){
					t++;
				}
				tmp.insert(tmp.begin()+t,(*i).second[j]);
			}
		}
		for(size_t j=0;j<tmp.size();j++){
			std::cout<<tmp[j].pos.chr.c_str()<<" "<<tmp[j].pos.start<<std::endl;
		}
	}
}
position get_pos(std::map<std::string, std::string> genome, int min_pos, int max_pos) {
	position tmp;
	std::map<std::string, std::string>::iterator chr = genome.begin();
	int id = rand() % genome.size(); //chose rand chr.
	int count = 0;
	while (chr != genome.end() && id != count) { //fast forward
		count++;
		chr++;
	}

	tmp.chr = (*chr).first;

	tmp.start = rand() % ((*chr).second.size() - max_pos); //choose random start pos within chr length
	if (max_pos == -1) {
		//insertion, translocation:
		tmp.stop = max_pos;
	} else {
		tmp.stop = tmp.start + min_pos + (rand() % (max_pos - min_pos)); // choose stop location
	}

	int num = 0;

	while ((*chr).second.size() < tmp.stop && num < 50) { //choose chr,start, stop such that the mutation fits. Allow max 50 iterations!
		tmp.start = rand() % (*chr).second.size(); //choose random start pos within chr length
		tmp.stop = tmp.start + min_pos + (rand() % (max_pos - min_pos)); // choose stop location
		num++;
	}
	if (num == 50) {
		std::cerr << "Simulations are hindered by the two small chr size. " << std::endl;
		tmp.stop = -2;
	}
	return tmp;
}
std::vector<struct_var> generate_mutations(std::string parameter_file, std::map<std::string, std::string> genome) {
	parameter par = parse_param(parameter_file);

	std::vector<struct_var> svs;
//duplications
	struct_var mut;
	for (int i = 0; i < par.dup_num; i++) {
		mut.type = 0;
		//get_start location;
		mut.pos = get_pos(genome, par.dup_min, par.dup_max);
		//get_opposit location
		svs.push_back(mut);
	}
	//indels
	for (int i = 0; i < par.indel_num; i++) {
		std::cout << "indel" << std::endl;
		if (rand() % 100 <= 50) {
			mut.type = 1;
			mut.pos = get_pos(genome, par.indel_min, par.indel_max);
			mut.target = mut.pos;
		} else {
			mut.type = 4;
			mut.pos = get_pos(genome, par.indel_min, par.indel_max);
			mut.target = mut.pos;
		}
		svs.push_back(mut);
	}
	//inv
	for (int i = 0; i < par.inv_num; i++) {
		std::cout << "inv" << std::endl;
		mut.type = 2;
		mut.pos = get_pos(genome, par.inv_min, par.inv_max);
		mut.target = mut.pos;
		mut.target.start = mut.pos.stop;
		mut.target.stop = mut.pos.start;
		svs.push_back(mut);
	}
	//tra
	for (int i = 0; i < par.translocations_num; i++) {
		std::cout << "tra" << std::endl;
		mut.type = 3;
		mut.pos = get_pos(genome, par.translocations_min, par.translocations_max);
		mut.target = get_pos(genome, par.translocations_min, par.translocations_max);
		while (strcmp(mut.target.chr.c_str(), mut.pos.chr.c_str()) == 0) {
			mut.target = get_pos(genome, par.translocations_min, par.translocations_max);
		}
		svs.push_back(mut);
	}

//	sort_svs(svs);
	return svs;
}
void store_ins(std::vector<insertions> & ins, insertions tmp) {

	std::vector<insertions>::iterator i = ins.begin();
	while (i != ins.end() && tmp.target.start > (*i).target.start) {
		i++;
	}
	ins.insert(i, tmp);
}
char complement(char nuc) {
	switch (nuc) {
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
		return nuc;
		break;
	}
	return nuc;
}
void invert(std::string &seq) {
	std::string tmp;
	for (std::string::reverse_iterator i = seq.rbegin(); i != seq.rend(); i++) {
		tmp += complement((*i));
	}
	seq.clear();

	seq = tmp;
}
std::string rand_seq(int length) {
	std::string tmp;
	tmp.resize((size_t) length);
	for (int i = 0; i < length; i++) {
		switch (rand() % 4) {
		case 0:
			tmp[i] = 'A';
			break;
		case 1:
			tmp[i] = 'C';
			break;
		case 2:
			tmp[i] = 'G';
			break;
		case 3:
			tmp[i] = 'T';
			break;
		}
	}
	return tmp;
}
void apply_mutations(std::map<std::string, std::string> &genome, std::vector<struct_var> svs, bool directions) {
	srand(time(NULL));
	std::vector<insertions> ins;
	insertions in;
	std::string seq1;
	std::string seq2;
	int pos;
	//all mutations that do not change the coordinates are later applied.
	//all others are directly applied (e.g. INV, TRA)
	for (size_t i = 0; i < svs.size(); i++) {
		std::string tmp;
		switch (svs[i].type) {
		case 0:
			//duplication
			svs[i].seq = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			in.seq = svs[i].seq;
			in.target = svs[i].pos;	//check
			store_ins(ins, in);
			break;
		case 1:
			//insertion:
			svs[i].seq = rand_seq(svs[i].target.stop - svs[i].target.start);
			in.seq = svs[i].seq;
			in.target = svs[i].target;
			store_ins(ins, in);
			break;
		case 2:
			//inversion
			tmp = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			//			cout << tmp << endl;
			invert(tmp);
			//			cout << tmp << endl;
			genome[svs[i].pos.chr].erase(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			genome[svs[i].pos.chr].insert(svs[i].pos.start, tmp);
			break;
		case 3:
			//translocations
			seq1 = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			seq2 = genome[svs[i].target.chr].substr(svs[i].target.start, (svs[i].target.stop - svs[i].target.start));

			pos = 0;
			for (int j = svs[i].target.start; j < svs[i].target.stop; j++) {
				genome[svs[i].target.chr][j] = seq1[pos];
				pos++;
			}
			pos = 0;
			for (int j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
				genome[svs[i].pos.chr][j] = seq2[pos];
				pos++;
			}
			break;
		case 4:
			svs[i].type = 4;
			//deletion: //just mark those regions
			std::cout << "size: " << genome[svs[i].pos.chr].size() << " " << svs[i].pos.start << " " << (svs[i].pos.stop - svs[i].pos.start) << std::endl;
			genome[svs[i].pos.chr].erase(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			genome[svs[i].pos.chr].insert(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start), 'X');
			break;
		default:
			break;
		}
	}
	if (directions) {
		for (std::vector<insertions>::reverse_iterator i = ins.rbegin(); i != ins.rend(); i++) {
			genome[(*i).target.chr].insert((*i).target.start, (*i).seq);
		}
	} else {
		for (std::vector<insertions>::iterator i = ins.begin(); i != ins.end(); i++) {
			genome[(*i).target.chr].insert((*i).target.start, (*i).seq);
		}
	}
}
void write_fasta(std::string output_prefix, std::map<std::string, std::string> genome) {
	std::string out = output_prefix;
	out += ".fasta";
	FILE *file2;
	file2 = fopen(out.c_str(), "w");
	if (file2 == NULL) {
		std::cout << "Error in printing: The file or path that you set " << output_prefix.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
		exit(0);
	}

	for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
		fprintf(file2, "%c", '>');
		fprintf(file2, "%s", (*i).first.c_str());
		fprintf(file2, "%c", '\n');
		int len = 0;
		for (size_t j = 1; j < (*i).second.size() + 1; j++) {
			if ((*i).second[j - 1] != 'X') {
				fprintf(file2, "%c", (*i).second[j - 1]);
				len++;
			}
			if (len % 100 == 0) {
				fprintf(file2, "%c", '\n');
			}

		}
		if (len % 100 != 0) {
			fprintf(file2, "%c", '\n');
		}
	}
	std::cout << std::endl;
	fclose(file2);
}
void write_sv(std::string output_prefix, std::vector<struct_var> svs) {
	std::string out = output_prefix;
	out += ".bed";
	FILE *file2;
	file2 = fopen(out.c_str(), "w");
	if (file2 == NULL) {
		std::cout << "Error in printing: The file or path that you set " << out.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
		exit(0);
	}

	out = output_prefix;
	out += ".insertions.fa";
	FILE *file;
	file = fopen(out.c_str(), "w");
	if (file == NULL) {
		std::cout << "Error in printing: The file or path that you set " << out.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
		exit(0);
	}

	for (size_t i = 0; i < svs.size(); i++) {
		if (svs[i].type == 1) { //write inserted sequeces to fasta file!
			fprintf(file, "%c", '>');
			fprintf(file, "%s", svs[i].pos.chr.c_str());
			fprintf(file, "%c", '_');
			fprintf(file, "%i", svs[i].pos.start);
			fprintf(file, "%c", '\n');
			fprintf(file, "%s", svs[i].seq.c_str());
			fprintf(file, "%c", '\n');
		}
		//write pseudo bed:
		if (svs[i].type == 3) {
			fprintf(file2, "%s", svs[i].pos.chr.c_str());
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%i", svs[i].pos.start);
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%s", svs[i].target.chr.c_str());
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%i", svs[i].target.start);
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%s", "TRA\n");

			fprintf(file2, "%s", svs[i].pos.chr.c_str());
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%i", svs[i].pos.stop);
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%s", svs[i].target.chr.c_str());
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%i", svs[i].target.stop);
			fprintf(file2, "%c", '\t');

		} else {
			fprintf(file2, "%s", svs[i].pos.chr.c_str());
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%i", svs[i].pos.start);
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%s", svs[i].pos.chr.c_str());
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%i", svs[i].pos.stop);
			fprintf(file2, "%c", '\t');
		}

		switch (svs[i].type) {
		case 0:
			fprintf(file2, "%s", "DUP");
			break;
		case 1:
			fprintf(file2, "%s", "INS");
			break;
		case 2:
			fprintf(file2, "%s", "INV");
			break;
		case 3:
			fprintf(file2, "%s", "TRA");
			break;
		case 4:
			fprintf(file2, "%s", "DEL");
			break;
		default:
			break;
		}
		fprintf(file2, "%c", '\n');
	}
	fclose(file);
	fclose(file2);
}

void simulate_SV(std::string ref_file, std::string parameter_file, bool coordinates, std::string output_prefix) {

	//read in list of SVs over vcf?
	//apply vcf to genome?
	srand(time(NULL));
	parameter par = parse_param(parameter_file);
	int min_chr_len = std::max(std::max(par.dup_max, par.indel_max), std::max(par.inv_max, par.translocations_max));
	std::map<std::string, std::string> genome = read_fasta(ref_file, min_chr_len);
	std::cout << "generate SV" << std::endl;
	std::vector<struct_var> svs = generate_mutations(parameter_file, genome);
	for (size_t i = 0; i < svs.size(); i++) {
		std::cout << svs[i].type << " " << svs[i].pos.chr.c_str() << " " << svs[i].pos.start << " " << svs[i].pos.stop << std::endl;
	}
	apply_mutations(genome, svs, coordinates);	//problem: We need two different coordinates. Simulate once for one and then for the other???
	std::cout << "write genome" << std::endl;
	write_fasta(output_prefix, genome);
	std::cout << "write SV" << std::endl;
	write_sv(output_prefix, svs);
}

void generate_parameter_file(std::string parameter_file) {
	FILE *file2;
	file2 = fopen(parameter_file.c_str(), "w");
	if (file2 == NULL) {
		std::cerr << "Error in printing: The file or path that you set " << parameter_file.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
		exit(0);
	}
	fprintf(file2, "%s", "PARAMETER FILE: DO JUST MODIFY THE VALUES AND KEEP THE SPACES!\n");
	fprintf(file2, "%s", "DUPLICATION_minimum_length: 100\n");
	fprintf(file2, "%s", "DUPLICATION_maximum_length: 10000\n");
	fprintf(file2, "%s", "DUPLICATION_number: 3\n");

	fprintf(file2, "%s", "INDEL_minimum_length: 20\n");
	fprintf(file2, "%s", "INDEL_maximum_length: 500\n");
	fprintf(file2, "%s", "INDEL_number: 1\n");

	fprintf(file2, "%s", "TRANSLOCATION_minimum_length: 1000\n");
	fprintf(file2, "%s", "TRANSLOCATION_maximum_length: 3000\n");
	fprintf(file2, "%s", "TRANSLOCATION_number: 2\n");

	fprintf(file2, "%s", "INVERSION_minimum_length: 600\n");
	fprintf(file2, "%s", "INVERSION_maximum_length: 800\n");
	fprintf(file2, "%s", "INVERSION_number: 4\n");

	fclose(file2);
}
