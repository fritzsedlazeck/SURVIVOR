/*
 * SV_Simulator.cpp
 *
 *  Created on: Jan 30, 2016
 *      Author: fsedlaze
 */

#include "SV_Simulator.h"

bool is_valid(char base) {
	return (((base == 'A' || base == 'C') || (base == 'R' || base == 'X')) || ((base == 'T' || base == 'G') || (base == 'N' || base == 'M')));
}

void check_genome(std::map<std::string, std::string> &genome, std::string msg) {
	std::cout << msg << " Genome checking:" << std::endl;

	for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
		for (size_t j = 1; j < (*i).second.size() + 1; j++) {
			//if (!is_valid((*i).second[j - 1])) {
			//std::cout << "err! " << (*i).second[j - 1] << std::endl;
			//}
		}
	}
}
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

float parse_value_float(char* buffer, size_t buffer_size) { //required for parameter!
	int count = 0;
	for (size_t i = 1; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
		if (count == 1) {
			return atof(&buffer[i]);
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
	tmp.dup_min = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.dup_max = parse_value(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);
	tmp.dup_max_amp= parse_value(buffer, buffer_size);
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
	myfile.getline(buffer, buffer_size);
	tmp.inv_del_num = 0;

	if (!myfile.eof()) {
		tmp.inv_del_min = parse_value(buffer, buffer_size);
		myfile.getline(buffer, buffer_size);
		tmp.inv_del_max = parse_value(buffer, buffer_size);
		myfile.getline(buffer, buffer_size);
		tmp.inv_del_num = parse_value(buffer, buffer_size);
		myfile.getline(buffer, buffer_size);
	}
	tmp.inv_dup_num = 0;
	if (!myfile.eof()) {
		tmp.inv_dup_min = parse_value(buffer, buffer_size);
		myfile.getline(buffer, buffer_size);
		tmp.inv_dup_max = parse_value(buffer, buffer_size);
		myfile.getline(buffer, buffer_size);
		tmp.inv_dup_num = parse_value(buffer, buffer_size);
		myfile.getline(buffer, buffer_size);
	}
	tmp.diploid=false;
	if (!myfile.eof()) {
		tmp.diploid=(bool)( parse_value(buffer,buffer_size)==2);
		myfile.getline(buffer, buffer_size);
		tmp.hom_rate= parse_value_float(buffer,buffer_size);
		myfile.getline(buffer, buffer_size);
	}

	tmp.intrachr_num = 0;
	/*if (!myfile.eof()) {
	 tmp.intrachr_min = parse_value(buffer, buffer_size);
	 myfile.getline(buffer, buffer_size);
	 tmp.intrachr_max = parse_value(buffer, buffer_size);
	 myfile.getline(buffer, buffer_size);
	 tmp.intrachr_num = parse_value(buffer, buffer_size);
	 //std::cout<<"NUM: "<<tmp.intrachr_num<<std::endl;
	 }*/
	myfile.close();
	return tmp;
}

std::map<std::string, std::string> read_fasta(std::string ref_file, int min_length) {
	//size_t buffer_size = 200000;
	//char*buffer = new char[buffer_size];
	std::string buffer;
	std::ifstream myfile;

	myfile.open(ref_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: " << ref_file.c_str() << std::endl;
		exit(0);
	}

	getline(myfile, buffer);
	std::map<std::string, std::string> genome;
	std::string seq;
	std::string name;
	while (!myfile.eof()) {
		if (buffer[0] == '>') {
			if ((int) seq.size() > min_length) {
				genome[name] = seq;
			}
			name.clear();
			seq.clear();

			for (size_t i = 1; i < buffer.size() && buffer[i] != '\n' && buffer[i] != '\0' && buffer[i] != ' '; i++) {
				name += (buffer[i]);
			}
		} else {
			for (size_t i = 0; i < buffer.size() && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
				seq += toupper(buffer[i]);
			}
		}
		getline(myfile, buffer);
	}
	for (size_t i = 0; i < buffer.size() && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
		seq += toupper(buffer[i]);
	}
	myfile.close();
	if ((int) seq.size() > min_length) {
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
				std::cout << (*i).second[j].pos.chr.c_str() << " " << (*i).second[j].pos.start << std::endl;
				size_t t = 0;
				while (tmp[t].pos.start < (*i).second[j].pos.start) {
					t++;
				}
				tmp.insert(tmp.begin() + t, (*i).second[j]);
			}
		}
		for (size_t j = 0; j < tmp.size(); j++) {
			std::cout << tmp[j].pos.chr.c_str() << " " << tmp[j].pos.start << std::endl;
		}
	}
}
float percent_N(std::string seq) {
	double n = 0;
	double size = (double) seq.size();
	for (size_t i = 0; i < seq.size(); i++) {
		if (seq[i] == 'N') {
			n++;
		}
	}
//	std::cout<<"Percent: "<<n/size << " "<<n <<" " <<size<<std::endl;
	return n / size;
}
position get_pos(std::map<std::string, std::string> genome, int min_pos, int max_pos) {
	//std::cout<<max_pos<<std::endl;
	position tmp;
	std::string seq = "N";
//	std::cout<<"Start:"<<std::endl;
	for (size_t i = 0; i < 100 && percent_N(seq) > 0.05; i++) {
		std::map<std::string, std::string>::iterator chr = genome.begin();
		int id = rand() % genome.size(); //chose rand chr.
		int count = 0;
		while (chr != genome.end() && id != count) { //fast forward
			count++;
			chr++;
		}

		tmp.chr = (*chr).first;

		tmp.start = rand() % (((*chr).second.size() - max_pos)); //choose random start pos within chr length
		if (max_pos == -1) {
			//insertion, translocation:
			tmp.stop = max_pos;
		} else {
			tmp.stop = tmp.start + min_pos + (rand() % (max_pos - min_pos)); // choose stop location
		}

		int num = 0;

		while ((int) (*chr).second.size() < tmp.stop && num < 100) { //choose chr,start, stop such that the mutation fits. Allow max 50 iterations!
			tmp.start = rand() % (((*chr).second.size() - max_pos)); //choose random start pos within chr length
			tmp.stop = tmp.start + min_pos + (rand() % (max_pos - min_pos)); // choose stop location
			num++;
		}
		if (num == 100) {
			std::cerr << "Simulations are hindered by the two small chr size. " << std::endl;
			tmp.stop = -2;
		}
		if (max_pos != -1) {
			seq = (*chr).second.substr(tmp.start, tmp.stop - tmp.start);
		}
		//std::cout<<(*chr).first<<" "<<tmp.start<<" "<<tmp.stop<<std::endl;
	}
	//std::cout<<"end:"<<std::endl;
	return tmp;
}

bool is_overlapping(position curr, std::vector<struct_var> svs) {

	for (size_t i = 0; i < svs.size(); i++) {
		if (strcmp(svs[i].pos.chr.c_str(), curr.chr.c_str()) == 0) {
			if (svs[i].pos.stop >= curr.start && svs[i].pos.start <= curr.stop) {
				return true;
			}
		}
	}
	return false;
}
position choose_pos(std::map<std::string, std::string> genome, int min, int max, std::vector<struct_var>& svs) {
	position pos = get_pos(genome, min, max);
	int num = 0;
	while (is_overlapping(pos, svs) && num < 30) {
		pos = get_pos(genome, min, max);
		num++;
	}
	if (num == 30) {
		std::cerr << "Terminate program as it could not find a non overlapping region" << std::endl;
		exit(0);
	}
	return pos;
}
std::vector<struct_var> generate_mutations(std::string parameter_file, std::map<std::string, std::string> genome) {
	parameter par = parse_param(parameter_file);
	std::vector<struct_var> svs;
//duplications
	struct_var mut;
	for (int i = 0; i < par.dup_num; i++) {
		mut.type = 0;
		//get_start location;
		mut.pos = choose_pos(genome, par.dup_min, par.dup_max, svs);
		//get_opposit location
		mut.copy_num=1;
		mut.print = true;
		svs.push_back(mut);
	}
	//indels
	for (int i = 0; i < par.indel_num; i++) {
		//std::cout << "indel" << std::endl;
		if (rand() % 100 <= 50) {
			mut.type = 1; //insertion
		} else {
			mut.type = 4; //deletion
		}
		mut.pos = choose_pos(genome, par.indel_min, par.indel_max, svs);
		mut.target = mut.pos;
		mut.copy_num=1;
		mut.print = true;
		svs.push_back(mut);
	}
	//inv
	for (int i = 0; i < par.inv_num; i++) {
		//	std::cout << "inv" << std::endl;
		mut.type = 2;
		mut.pos = choose_pos(genome, par.inv_min, par.inv_max, svs);
		mut.target = mut.pos;
		mut.target.start = mut.pos.stop;
		mut.target.stop = mut.pos.start;
		mut.copy_num=1;
		mut.print = true;
		svs.push_back(mut);
	}

	//tra inter
	/*std::cout << par.intrachr_num << std::endl;
	 for (int i = 0; i < par.intrachr_num; i++) {
	 mut.type = 6;
	 mut.pos = choose_pos(genome, par.translocations_min, par.translocations_max, svs);
	 //std::cout<<i<<": "<<mut.pos.chr<<" "<<mut.pos.start<<" size: "<<mut.pos.stop-mut.pos.start<<std::endl;
	 mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs); //TRA has to be of the same size!
	 while (strcmp(mut.target.chr.c_str(), mut.pos.chr.c_str()) != 0) {
	 mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs);
	 }

	 //I need to be sure about the same lenght of the tra!:
	 int size1 = mut.pos.stop - mut.pos.start;
	 int size2 = mut.target.stop - mut.target.start;

	 mut.pos.stop = mut.pos.start + std::min(size1, size2);
	 mut.target.stop = mut.target.start + std::min(size1, size2);
	 svs.push_back(mut);
	 }*/
	//tra
	for (int i = 0; i < par.translocations_num; i++) {
		//	std::cout << "tra" << std::endl;
		mut.type = 3;
		mut.pos = choose_pos(genome, par.translocations_min, par.translocations_max, svs);
		//std::cout<<"size: "<<mut.pos.stop-mut.pos.start<<std::endl;
		mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs); //TRA has to be of the same size!
		while (strcmp(mut.target.chr.c_str(), mut.pos.chr.c_str()) == 0) {
			mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs);
		}

		//I need to be sure about the same lenght of the tra!:
		int size1 = mut.pos.stop - mut.pos.start;
		int size2 = mut.target.stop - mut.target.start;

		mut.pos.stop = mut.pos.start + std::min(size1, size2);
		mut.target.stop = mut.target.start + std::min(size1, size2);
		mut.copy_num=1;
		mut.print = true;
		svs.push_back(mut);
	}
	//complex inv_del
	for (int i = 0; i < par.inv_del_num; i++) {
		//1. sim:
		mut.type = 2;
		mut.pos = choose_pos(genome, par.inv_del_min, par.inv_del_max, svs);
		//2. determin size of del:
		int len = (mut.pos.stop - mut.pos.start) / 10; //dels are ~20% ofthe size!
		mut.pos.start += len;
		mut.pos.stop -= len;
		mut.target = mut.pos;
		mut.target.start = mut.pos.stop;
		mut.target.stop = mut.pos.start;
		mut.print = true;
		mut.copy_num=1;
		svs.push_back(mut);

		struct_var del;
		//the del infront:
		del.type = 4;
		del.pos.chr = mut.pos.chr;
		del.pos.stop = mut.pos.start;
		del.pos.start = del.pos.stop - len;
		del.target = del.pos;
		mut.print = true;
		mut.copy_num=1;
		svs.push_back(del);

		//the del behind:
		del.pos.start = mut.pos.stop;
		del.pos.stop = del.pos.start + len;
		del.target = del.pos;
		mut.print = true;
		mut.copy_num=1;
		svs.push_back(del);
	}
	//inv dup
	for (int i = 0; i < par.inv_dup_num; i++) {
		mut.type = 5;
		//get_start location;
		mut.pos = choose_pos(genome, par.inv_dup_min, par.inv_dup_max, svs);
		//get_opposit location
		mut.print = true;
		mut.copy_num=1;
		svs.push_back(mut);
	}
//	sort_svs(svs);
	return svs;
}

std::vector<struct_var> generate_mutations_ref(std::string parameter_file, std::map<std::string, std::string> genome) {
	parameter par = parse_param(parameter_file);
	std::vector<struct_var> svs;
//duplications
	struct_var mut;
	//indels
	for (int i = 0; i < par.indel_num; i++) {
		//std::cout << "indel" << std::endl;
		if (rand() % 100 <= 50) {
			mut.type = 1; //insertion
		} else {
			mut.type = 4; //deletion
		}
		mut.pos = choose_pos(genome, par.indel_min, par.indel_max, svs);
		mut.target = mut.pos;
		mut.copy_num=1;
		svs.push_back(mut);
	}
	//inv
	for (int i = 0; i < par.inv_num; i++) {
		//	std::cout << "inv" << std::endl;
		mut.type = 2;
		mut.pos = choose_pos(genome, par.inv_min, par.inv_max, svs);
		mut.target = mut.pos;
		mut.target.start = mut.pos.stop;
		mut.target.stop = mut.pos.start;
		mut.copy_num=1;
		svs.push_back(mut);
	}
	//tra
	for (int i = 0; i < par.translocations_num; i++) {
		//	std::cout << "tra" << std::endl;
		mut.type = 3;
		mut.pos = choose_pos(genome, par.translocations_min, par.translocations_max, svs);
		//std::cout<<"size: "<<mut.pos.stop-mut.pos.start<<std::endl;
		mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs); //TRA has to be of the same size!
		while (strcmp(mut.target.chr.c_str(), mut.pos.chr.c_str()) == 0) {
			mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs);
		}

		//I need to be sure about the same lenght of the tra!:
		int size1 = mut.pos.stop - mut.pos.start;
		int size2 = mut.target.stop - mut.target.start;

		mut.pos.stop = mut.pos.start + std::min(size1, size2);
		mut.target.stop = mut.target.start + std::min(size1, size2);
		mut.copy_num=1;
		svs.push_back(mut);
	}

	return svs;
}

void store_sorted(std::vector<struct_var> &svs, struct_var tmp) {
	std::vector<struct_var>::iterator i = svs.begin();
	while (i != svs.end() && (tmp.target.start > (*i).target.start)) {
		i++;
	}
	svs.insert(i, tmp);
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
	//tmp.resize((size_t) length);
	for (int i = 0; i < length; i++) {
		switch (rand() % 4) {
		case 0:
			tmp += 'A';
			break;
		case 1:
			tmp += 'C';
			break;
		case 2:
			tmp += 'G';
			break;
		case 3:
			tmp += 'T';
			break;
		}
	}
	return tmp;
}
void apply_mutations(std::map<std::string, std::string> &genome, std::vector<struct_var>& svs,int max_dup_amp) {
	srand(time(NULL));
	std::vector<insertions> ins;
	insertions in;
	std::string seq1;
	std::string seq2;
	int pos;


	std::vector<struct_var> new_svs; //Thanks to the invdup we need this.
	//all mutations that do not change the coordinates are later applied.
	//all others are directly applied (e.g. INV, TRA)
	for (size_t i = 0; i < svs.size(); i++) {
		std::string tmp;
		int num=1;
		switch (svs[i].type) {
		case 0:
			//duplication
			svs[i].seq = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			in.seq = svs[i].seq;
			in.target = svs[i].pos;	//check
			num+=rand() % max_dup_amp + 1;
			svs[i].copy_num=num;
			while(num>0){
				store_ins(ins, in);
				num--;
			}

			break;
		case 1:
			//insertion:
			svs[i].seq = rand_seq(svs[i].target.stop - svs[i].target.start);
			svs[i].ref = genome[svs[i].pos.chr][svs[i].pos.start];
			in.seq = svs[i].seq;
			in.target = svs[i].target;
			store_ins(ins, in);
			break;
		case 2:
			//inversion
			tmp = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			//std::cout<<"INV: "<<tmp.size()<<std::endl;
			invert(tmp);
			genome[svs[i].pos.chr].erase(svs[i].pos.start, tmp.size());
			genome[svs[i].pos.chr].insert(svs[i].pos.start, tmp);
			break;
		case 3:
			//translocations
			seq1 = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			seq2 = genome[svs[i].target.chr].substr(svs[i].target.start, (svs[i].target.stop - svs[i].target.start));
			//	std::cout<<"TRA: "<<seq1.size()<<" "<<seq2.size()<<std::endl;

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
		case 4: //deletion: //just mark those regions
			//std::cout<<"DEL: "<<svs[i].pos.chr<<" "<<svs[i].pos.start <<" "<<svs[i].pos.stop<<" g: "<< genome[svs[i].pos.chr].size()<<std::endl;
			//	std::cout << "size: " << genome[svs[i].pos.chr].size() << " " << svs[i].pos.start << " " << (svs[i].pos.stop - svs[i].pos.start) << std::endl;
			svs[i].seq = genome[svs[i].pos.chr][svs[i].pos.start];
			svs[i].ref = "";
			for (int j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
				svs[i].ref += genome[svs[i].pos.chr][j];
				genome[svs[i].pos.chr][j] = 'X';
			}
			break;

		case 5:
			//invduplication first a duplication and then an inversion
			svs[i].seq = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			in.seq = svs[i].seq;
			invert(in.seq);
			in.target = svs[i].pos;	//check
			store_ins(ins, in);
			svs[i].type = 0;	//dup
			new_svs.push_back(svs[i]);
			svs[i].type = 2;	//inv
			break;
		case 6:
			//inter tra:
			svs[i].seq = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
			// = genome[svs[i].target.chr].substr(svs[i].target.start, (svs[i].target.stop - svs[i].target.start));
			//	std::cout<<"TRA: "<<seq1.size()<<" "<<seq2.size()<<std::endl;

			in.seq = svs[i].seq;
			in.target = svs[i].target;
			store_ins(ins, in);
			pos = 0;
			for (int j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
				genome[svs[i].pos.chr][j] = 'X';
				pos++;
			}
			break;
		default:
			break;
		}
	}

	//so far coorindates remain the same!
	for (std::vector<insertions>::reverse_iterator i = ins.rbegin(); i != ins.rend(); i++) { /// sequence is getting added!
		genome[(*i).target.chr].insert((*i).target.start, (*i).seq);
	}
	for (size_t i = 0; i < new_svs.size(); i++) {
		svs.push_back(new_svs[i]);
	}

}

void apply_mutations_ref(std::map<std::string, std::string> &genome, std::vector<struct_var> &svs) {
	std::cout << "apply mut ref!" << std::endl;
	srand(time(NULL));
	std::vector<insertions> ins;
	insertions in;
	std::string seq1;
	std::string seq2;
	int pos;
	//all mutations that do not change the coordinates are later applied.
	//all others are directly applied (e.g. INV, TRA)
	//for (size_t i = 0; i < svs.size(); i++) {
	std::vector<struct_var> tmp_svs;
	for (size_t i = 0; i < svs.size(); i++) {
		store_sorted(tmp_svs, svs[i]);
	}
	svs = tmp_svs;	//stored according to postions 0-> max

	//switch ins vs. del
	for (size_t i = 0; i < svs.size(); i++) {
		if (svs[i].type == 4) {
			svs[i].type = 1;
		} else if (svs[i].type == 1) {
			svs[i].type = 4;
		}
	}

	for (std::vector<struct_var>::iterator i = svs.begin(); i != svs.end(); i++) {
		std::cout << "apply: " << (*i).pos.chr << " " << (*i).pos.start << " " << (*i).type << std::endl;
		std::string tmp;
		switch ((*i).type) {
		case 1:
			genome[(*i).pos.chr].erase((*i).pos.start, (*i).pos.stop - (*i).pos.start);
			break;
		case 2:
			//inversion
			tmp = genome[(*i).pos.chr].substr((*i).pos.start, ((*i).pos.stop - (*i).pos.start));
			//std::cout<<"INV: "<<tmp.size()<<std::endl;
			invert(tmp);
			genome[(*i).pos.chr].erase((*i).pos.start, tmp.size());
			genome[(*i).pos.chr].insert((*i).pos.start, tmp);
			break;

		case 4:
			//deletions: (simulated insertions)
			(*i).seq = rand_seq((*i).target.stop - (*i).target.start);
			in.seq = (*i).seq;
			in.target = (*i).target;
			genome[in.target.chr].insert(in.target.start, in.seq);
			break;

		default:
			break;
		}
	}
	for (std::vector<struct_var>::iterator i = svs.begin(); i != svs.end(); i++) {
		if ((*i).type == 3) {
			//translocations
			seq1 = genome[(*i).pos.chr].substr((*i).pos.start, ((*i).pos.stop - (*i).pos.start));
			seq2 = genome[(*i).target.chr].substr((*i).target.start, ((*i).target.stop - (*i).target.start));
			//	std::cout<<"TRA: "<<seq1.size()<<" "<<seq2.size()<<std::endl;
			pos = 0;
			for (int j = (*i).target.start; j < (*i).target.stop; j++) {
				genome[(*i).target.chr][j] = seq1[pos];
				pos++;
			}
			pos = 0;
			for (int j = (*i).pos.start; j < (*i).pos.stop; j++) {
				genome[(*i).pos.chr][j] = seq2[pos];
				pos++;
			}
		}
	}

}
void write_fasta(std::string output_prefix, std::map<std::string, std::string> genome, bool append) {
	std::string out = output_prefix;
	out += ".fasta";
	FILE *file2;
	if (!append) {
		file2 = fopen(out.c_str(), "w");
	} else {
		file2 = fopen(out.c_str(), "a");
	}
	if (file2 == NULL) {
		std::cout << "Error in printing: The file or path that you set " << output_prefix.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
		exit(0);
	}

	bool flag = false;
	for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
		fprintf(file2, "%c", '>');
		fprintf(file2, "%s", (*i).first.c_str());
		fprintf(file2, "%c", '\n');
		int len = 0;
		for (size_t j = 1; j < (*i).second.size() + 1; j++) {
			//if (!is_valid((*i).second[j - 1])) {
			//		std::cout << "err! " << (*i).second[j - 1] << std::endl;
			//	}
			if ((*i).second[j - 1] != 'X') {
				fprintf(file2, "%c", (*i).second[j - 1]);
				len++;
				flag = true;
			}
			if (len % 100 == 0 && flag) {
				flag = false;
				fprintf(file2, "%c", '\n');
			}

		}
		if (len % 100 != 0) {
			fprintf(file2, "%c", '\n');
		}
	}
	//std::cout << std::endl;
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
			if (svs[i].type == 3) {
				fprintf(file2, "%s", "TRA\n");
			}
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
		case 5:
			fprintf(file2, "%s", "INVDUP");
			break;
		default:
			break;
		}
		fprintf(file2, "%c", '\n');
	}
	fclose(file);
	fclose(file2);
}

char mut_char(char old) {
	std::string nucs = "ACGT";
	int index = rand() % 4;
	switch (old) {
	case 'A':
		while (index == 0) {
			index = rand() % 4;
		}
		return nucs[index];
		break;
	case 'C':
		while (index == 1) {
			index = rand() % 4;
		}
		return nucs[index];
		break;
	case 'G':
		while (index == 2) {
			index = rand() % 4;
		}
		return nucs[index];
		break;
	case 'T':
		while (index == 3) {
			index = rand() % 4;
		}
		return nucs[index];
		break;
	default:
		return old;
		break;
	}
	return old;
}
const std::string currentDateTime2() {
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
	return buf;
}
void print_vcf_header2(FILE *&file, std::map<std::string, std::string> genome) {
	fprintf(file, "%s", "##fileformat=VCFv4.2\n");
	fprintf(file, "%s", "##source=Sniffles\n");
	std::string time = currentDateTime2();
	fprintf(file, "%s", "##fileDate=");
	fprintf(file, "%s", time.c_str());

	//REport over all chrs:
	for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
		fprintf(file, "%s", "\n");
		fprintf(file, "%s", "##contig=<ID=");
		fprintf(file, "%s", (*i).first.c_str());
		fprintf(file, "%s", ",length=");
		fprintf(file, "%i", (int) (*i).second.size());
		fprintf(file, "%c", '>');
	}

	fprintf(file, "%s", "\n");
	fprintf(file, "%s", "##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(file, "%s", "##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(file, "%s", "##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(file, "%s", "##ALT=<ID=INVDUP,Description=\"InvertedDUP with unknown boundaries\">\n");
	fprintf(file, "%s", "##ALT=<ID=TRA,Description=\"Translocation\">\n");
	fprintf(file, "%s", "##ALT=<ID=INS,Description=\"Insertion\">\n");
	fprintf(file, "%s", "##FILTER=<ID=PASS,Description=\"PASS variant.\">\n");
	fprintf(file, "%s", "##FILTER=<ID=UNRESOLVED,Description=\"An insertion that is longer than the read and thus we cannot predict the full size.\">\n");

	fprintf(file, "%s", "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n");
	fprintf(file, "%s", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=AF,Number=.,Type=Integer,Description=\"Allele Frequency.\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");

	fprintf(file, "%s", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n");
}
void print_snp_vcf(std::string chr, int pos, char old_allele, char new_allele, FILE *&file, int id, std::string gt) {
	std::ostringstream convert;   // stream used for the conversion
	convert << chr;
	convert << "\t";
	convert << pos;      // insert the textual representation of 'Number' in the characters in the stream
	convert << "\t";
	convert << "SNP";
	convert << id;
	convert << "SURVIVOR\t";
	convert << old_allele;
	convert << "\t";
	convert << new_allele;

	convert << "\t.\tPASS\tPRECISE;SVMETHOD=SURVIVOR_sim;SVLEN=1\tGT\t";
	convert << gt;

	fprintf(file, "%s", convert.str().c_str());
	fprintf(file, "%c", '\n');
}

std::string print_vcf_sv(std::string chr, int pos, int type_id, std::string end_chr, int end_pos, int id, std::string seq, std::string ref, std::string gt,int dup_amp) {
	std::string type = "";
	switch (type_id) {
	case 0:
		type = "DUP";
		break;
	case 1:
		type = "INS";
		break;
	case 2:
		type = "INV";
		break;
	case 4:
		type = "DEL";
		break;
	case 5:
		type = "INVDUP";
		break;
	default:
		break;
	}

	std::ostringstream convert;   // stream used for the conversion
	convert << chr;
	convert << "\t";
	convert << pos;      // insert the textual representation of 'Number' in the characters in the stream
	convert << "\t";
	convert << type;
	convert << id;
	convert << "SURVIVOR\t";
	if (strncmp(type.c_str(), "INS", 3) == 0 || strncmp(type.c_str(), "DEL", 3) == 0) {
		convert << ref;
		convert << "\t";
		convert << seq;
//convert << "\t";
	} else {
		convert << "N\t<";
		convert << type;
		convert << ">";
	}

	convert << "\t.\tPASS\tPRECISE;SVTYPE=";
	convert << type;
	convert << ";SVMETHOD=SURVIVOR_sim;CHR2=";
	convert << end_chr;
	convert << ";END=";
	convert << end_pos;
	convert << ";SVLEN=";
	convert << end_pos - pos;
	if(type_id==0){
		convert << ";dup_num=";
		convert << dup_amp;
	}
	convert << "\tGT\t";
	convert << gt;
	convert << "\n";
	return convert.str();
}
void print_vcf_svs(FILE *& file, std::vector<struct_var> svs, int id) {
	for (size_t i = 0; i < svs.size(); i++) {
		//write pseudo bed:
		if (svs[i].type == 3) {
			fprintf(file, "%s", print_vcf_sv(svs[i].pos.chr, svs[i].pos.start, svs[i].type, svs[i].target.chr, svs[i].target.start, id, svs[i].seq, svs[i].ref, "1/1",svs[i].copy_num).c_str());
			fprintf(file, "%s", print_vcf_sv(svs[i].pos.chr, svs[i].pos.stop, svs[i].type, svs[i].target.chr, svs[i].target.stop, id, svs[i].seq, svs[i].ref, "1/1",svs[i].copy_num).c_str());
		} else { // all other types:

			fprintf(file, "%s", print_vcf_sv(svs[i].pos.chr, svs[i].pos.start, svs[i].type, svs[i].pos.chr, svs[i].pos.stop, id, svs[i].seq, svs[i].ref, "1/1",svs[i].copy_num).c_str());
		}
		id++;
	}
}
std::vector<struct_var> mut_snv(std::map<std::string, std::string> & genome, float snp_freq, FILE *&file2, int &id, int hap, float hom_rat) {
	std::cout << "generate SNP" << std::endl;
	std::vector<struct_var> snv;
	struct_var tmp;
	for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
		for (size_t pos = 0; pos < (*i).second.size(); pos++) {
			if ((*i).second[pos] != 'X' && (*i).second[pos] != 'N') {  // Do not mutate the nuc when it is X or N.
				float x = ((float) rand() / (float) (RAND_MAX));
				if (x < snp_freq) {
					char new_nuc = mut_char(toupper((*i).second[pos]));

					if (hap != -1) { //otherwise we dont need this!
						x = ((float) rand() / (float) (RAND_MAX));
						if (x < hom_rat) {
							tmp.pos.chr = (*i).first;
							tmp.pos.start = pos;
							tmp.seq = new_nuc;
							tmp.ref = (*i).second[pos];
							snv.push_back(tmp);
							print_snp_vcf(tmp.pos.chr, tmp.pos.start + 1, tmp.ref[0], tmp.seq[0], file2, id, "1/1");
						} else {
							print_snp_vcf((*i).first, pos + 1, (*i).second[pos], new_nuc, file2, id, "0|1"); // The SNP positions in VCF file are 1-based.
						}
					} else {
						print_snp_vcf((*i).first, pos + 1, (*i).second[pos], new_nuc, file2, id, "1/1");
					}
					id++;
					(*i).second[pos] = new_nuc;
				}
				if ((*i).second[pos] == ' ') {
					std::cout << "ERR" << std::endl;
				}
			}
		}
	}
	return snv;
}
void simulate_SV(std::string ref_file, std::string parameter_file, float snp_freq, bool coordinates, std::string output_prefix) {


	//read in list of SVs over vcf?
	//apply vcf to genome?
	srand(time(NULL));
	parameter par = parse_param(parameter_file);
	int min_chr_len = std::max(std::max(par.dup_max, par.indel_max), std::max(par.inv_max, par.translocations_max));
	std::map<std::string, std::string> genome = read_fasta(ref_file, min_chr_len);

	if (par.translocations_num > 0 && genome.size() < 2) {
		std::cerr << "We cannot simulate translocations without a second chromosome" << std::endl;
		exit(1);
	}

	if (genome.empty()) {
		std::cerr << "No Chromosomes found that are larger than min length of " << min_chr_len << ". Check parameter file to reduce the size of SV" << std::endl;
		exit(1);
	}

	std::string out = output_prefix;
	out += ".vcf";
	FILE *file2;
	file2 = fopen(out.c_str(), "w");
	if (file2 == NULL) {
		std::cout << "Error in printing: The file or path that you set " << out.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
		exit(0);
	}
	print_vcf_header2(file2, genome);

	/// Start sim of SV:
	std::cout << "generate SV over chrs:" << genome.size() << std::endl;
	float snv_ration = (snp_freq * par.hom_rate);
	std::vector<struct_var> svs;
	std::vector<struct_var> snv;
	int id = 0;
	if (coordinates) {	//simulate reads
		svs = generate_mutations(parameter_file, genome); // select coordinates and mark genome.
		if (!par.diploid) {

			snv = mut_snv(genome, (snv_ration * 2), file2, id, -1, -1); // apply SNV (no indels) , div 2 since this is only 1 hap.
		} else {
			snv = mut_snv(genome, snp_freq, file2, id, 1, par.hom_rate);
		}
		apply_mutations(genome, svs,par.dup_max_amp);	// apply changes to the genome (e.g. insertions).
	} else { //real read fake ref!
		svs = generate_mutations_ref(parameter_file, genome);
		apply_mutations_ref(genome, svs);
		mut_snv(genome, snp_freq, file2, id, -1, -1); // apply SNV (no indels) -1: cannot be diplod!
	}
	//check_genome(genome, "Post SV simulation");

	/// Start writing out:
	std::cout << "write genome" << std::endl;
	write_fasta(output_prefix, genome, false);
	std::cout << "write SV" << std::endl;

	if (par.diploid) {
		std::cout << "Simulate diploid!" << std::endl;

		genome = read_fasta(ref_file, min_chr_len);
		if (!coordinates) {
			std::cerr << "This is not possible" << std::endl;
			exit(0);
		}
		std::vector<struct_var> svs2;
		svs2 = generate_mutations(parameter_file, genome); // select coordinates and mark genome.

		//TODO swap a few svs vs. svs2!

		for (size_t i = 0; i < svs2.size(); i++) {
			float x = ((float) rand() / (float) (RAND_MAX));
			if (x < (par.hom_rate)) { // homozygous SV
				svs2[i] = svs[i];
				//svs[i].gt = "1/1";
				if (svs[i].type == 3) {
					fprintf(file2, "%s", print_vcf_sv(svs[i].pos.chr, svs[i].pos.start, svs[i].type, svs[i].target.chr, svs[i].target.start, id, svs[i].seq, svs[i].ref, "1/1",svs[i].copy_num).c_str());
					fprintf(file2, "%s", print_vcf_sv(svs[i].pos.chr, svs[i].pos.stop, svs[i].type, svs[i].target.chr, svs[i].target.stop, id, svs[i].seq, svs[i].ref, "1/1",svs[i].copy_num).c_str());
				} else { // all other types:
					fprintf(file2, "%s", print_vcf_sv(svs[i].pos.chr, svs[i].pos.start, svs[i].type, svs[i].pos.chr, svs[i].pos.stop, id, svs[i].seq, svs[i].ref, "1/1",svs[i].copy_num).c_str());
				}
				svs2[i].print=false;
			} else { // het SVs!
				if (svs[i].type == 3) {
					fprintf(file2, "%s", print_vcf_sv(svs[i].pos.chr, svs[i].pos.start, svs[i].type, svs[i].target.chr, svs[i].target.start, id, svs[i].seq, svs[i].ref, "0|1",svs[i].copy_num).c_str());
					fprintf(file2, "%s", print_vcf_sv(svs[i].pos.chr, svs[i].pos.stop, svs[i].type, svs[i].target.chr, svs[i].target.stop, id, svs[i].seq, svs[i].ref, "0|1",svs[i].copy_num).c_str());
				} else { // all other types:
					fprintf(file2, "%s", print_vcf_sv(svs[i].pos.chr, svs[i].pos.start, svs[i].type, svs[i].pos.chr, svs[i].pos.stop, id, svs[i].seq, svs[i].ref, "0|1",svs[i].copy_num).c_str());
				}
				id++;

			}
		}

		//SNV:
		for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
			for (size_t pos = 0; pos < (*i).second.size(); pos++) {
				if ((*i).second[pos] != 'X' && (*i).second[pos] != 'N') {  // Do not mutate the nuc when it is X or N.
					float x = ((float) rand() / (float) (RAND_MAX));
					if (x < (snv_ration)) { // Div 2 since this is only on the other hap!
						char new_nuc = mut_char(toupper((*i).second[pos]));
						print_snp_vcf((*i).first, pos + 1, (*i).second[pos], new_nuc, file2, id, "1|0"); // The SNP positions in VCF file are 1-based.
						id++;
						(*i).second[pos] = new_nuc;
					}
					if ((*i).second[pos] == ' ') {
						std::cout << "ERR" << std::endl;
					}
				}
			}
		}
		for (size_t i = 0; i < snv.size(); i++) {
			genome[snv[i].pos.chr][snv[i].pos.start] = snv[i].seq[0];
			//print_snp_vcf(snv[i].pos.chr, snv[i].pos.start + 1, snv[i].ref[0], snv[i].seq[0], file2, id, "1/1");
		}

		//	mut_snv(genome, snp_freq, file2, id, 2,hom_rate); // apply SNV (no indels)
		apply_mutations(genome, svs2,par.dup_max_amp);

		for (size_t i = 0; i < svs2.size(); i++) {
			if (svs2[i].print) {
				if (svs2[i].type == 3) {
					fprintf(file2, "%s", print_vcf_sv(svs2[i].pos.chr, svs2[i].pos.start, svs2[i].type, svs2[i].target.chr, svs2[i].target.start, id, svs2[i].seq, svs2[i].ref, "1|0",svs2[i].copy_num).c_str());
					fprintf(file2, "%s", print_vcf_sv(svs2[i].pos.chr, svs2[i].pos.stop, svs2[i].type, svs2[i].target.chr, svs2[i].target.stop, id, svs2[i].seq, svs2[i].ref, "1|0",svs2[i].copy_num).c_str());
				} else { // all other types:
					fprintf(file2, "%s", print_vcf_sv(svs2[i].pos.chr, svs2[i].pos.start, svs2[i].type, svs2[i].pos.chr, svs2[i].pos.stop, id, svs2[i].seq, svs2[i].ref, "1|0",svs2[i].copy_num).c_str());
				}
				id++;
			}

		}

		std::cout << "write genome" << std::endl;
		write_fasta(output_prefix, genome, true);
		std::cout << "write SV" << std::endl;
		//	write_sv(output_prefix, svs2);
		write_sv(output_prefix, svs);
		write_sv(output_prefix, svs2);
	} else {
		write_sv(output_prefix, svs);
		print_vcf_svs(file2, svs, id);

	}
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
	fprintf(file2, "%s", "DUPLICATION_maximum_num: 5\n");
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

	fprintf(file2, "%s", "INV_del_minimum_length: 600\n");
	fprintf(file2, "%s", "INV_del_maximum_length: 800\n");
	fprintf(file2, "%s", "INV_del_number: 2\n");

	fprintf(file2, "%s", "INV_dup_minimum_length: 600\n");
	fprintf(file2, "%s", "INV_dup_maximum_length: 800\n");
	fprintf(file2, "%s", "INV_dup_number: 2\n");

	fprintf(file2, "%s", "Number_haploid: 1\n");
	fprintf(file2, "%s", "homozygous_ratio: 0.6\n");

	fclose(file2);
}
