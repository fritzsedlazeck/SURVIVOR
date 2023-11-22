/*
 * combine_svs.cpp
 *
 *  Created on: Jul 6, 2016
 *      Author: fsedlaze
 */

#include "combine_svs.h"

bool is_valid_id(std::string id) {
	if (id.empty()) {
		return false;
	}
//	int count = 0;
	for (size_t i = 0; i < id.size(); i++) {
		if (id[i] == '-') {
			return false;
		}
	}
	return true;
}

bool is_valid_alleles(std::pair<std::string, std::string> alleles) {

	//search for a <> tag in second:
	for (size_t i = 0; i < alleles.second.size(); i++) {
		if (alleles.second[i] == '<' || alleles.second[i] == ':') {
			return false;
		}
	}
	return true;
}

breakpoint_str convert_position(strcoordinate pos) {
	breakpoint_str tmp;
	tmp.chr = pos.chr;
	tmp.position = pos.pos;
	return tmp;
}

std::string get_support_vec(std::vector<Support_Node *> caller_info) {

	std::string ss;
	ss.resize(Parameter::Instance()->max_caller, '0');
	for (size_t i = 0; i < caller_info.size(); i++) {
		//if (strncmp(caller_info[i]->genotype.c_str(), "0/0", 3) != 0) {
		ss[caller_info[i]->id] = '1';
		//}
	}
	return ss;
}
int get_support(std::vector<Support_Node *> caller_info) {
	return caller_info.size();

	int support = 0;
	for (size_t i = 0; i < caller_info.size(); i++) {
		//	std::cout<<"GO:"<<caller_info[i]->genotype<<std::endl;
		if (strncmp(caller_info[i]->genotype.c_str(), "0/0", 3) != 0) {
			//if (!caller_info[i]->starts.empty()) {
			support++;
		}
	}
	return support;
}
const std::string currentDateTime() {
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
	return buf;
}
void print_header(FILE *& file, std::vector<std::string> names, std::map<std::string, int> chrs) {
	std::map<std::string, int> names_checked;
	for (size_t i = 0; i < names.size(); i++) { // To aviod reporting the same identifier twice!
		if (names_checked.find(names[i]) == names_checked.end()) {
			names_checked[names[i]] = 1;

		} else {
			std::stringstream new_name;
			new_name << names[i];
			new_name << "_";
			new_name << names_checked[names[i]];
			names_checked[names[i]]++;
			names[i] = new_name.str();
		}
	}
	fprintf(file, "%s", "##fileformat=VCFv4.1\n");
	fprintf(file, "%s", "##source=SURVIVOR\n");
	std::string time = currentDateTime();
	fprintf(file, "%s", "##fileDate=");
	fprintf(file, "%s", time.c_str());
	fprintf(file, "%s", "\n"); //TODO change!
	for (std::map<std::string, int>::iterator i = chrs.begin(); i != chrs.end(); i++) {
		//	std::cout<<"Header: "<<(*i).first<<" "<<(*i).second<<std::endl;
		fprintf(file, "%s", "##contig=<ID=");
		fprintf(file, "%s", (*i).first.c_str());
		if ((*i).second > 0) {
			fprintf(file, "%s", ",length=");
			fprintf(file, "%i", (*i).second);
		}
		fprintf(file, "%s", ">\n");
	}
	fprintf(file, "%s", "##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(file, "%s", "##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(file, "%s", "##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(file, "%s", "##ALT=<ID=BND,Description=\"Translocation\">\n");
	fprintf(file, "%s", "##ALT=<ID=INS,Description=\"Insertion\">\n");
	fprintf(file, "%s", "##INFO=<ID=CIEND,Number=2,Type=String,Description=\"PE confidence interval around END\">\n");
	fprintf(file, "%s", "##INFO=<ID=CIPOS,Number=2,Type=String,Description=\"PE confidence interval around POS\">\n");
	fprintf(file, "%s", "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n");
	fprintf(file, "%s", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">\n");
	fprintf(file, "%s", "##INFO=<ID=RE,Number=1,Type=Integer,Description=\"read support\">\n");
	fprintf(file, "%s", "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Method for generating this merged VCF file.\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV.\">\n");
	fprintf(file, "%s", "##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description=\"Vector of supporting samples.\">\n");
	fprintf(file, "%s", "##INFO=<ID=SUPP,Number=1,Type=String,Description=\"Number of samples supporting the variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint.\">\n");

	fprintf(file, "%s", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=PSV,Number=1,Type=String,Description=\"Previous support vector\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=LN,Number=1,Type=Integer,Description=\"predicted length\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DR,Number=2,Type=Integer,Description=\"# supporting reference,variant reads in that order\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=ST,Number=1,Type=String,Description=\"Strand of SVs\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=QV,Number=1,Type=String,Description=\"Quality values: if not defined a . otherwise the reported value.\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=TY,Number=1,Type=String,Description=\"Types\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=ID,Number=1,Type=String,Description=\"Variant ID from input.\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=RAL,Number=1,Type=String,Description=\"Reference allele sequence reported from input.\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=AAL,Number=1,Type=String,Description=\"Alternative allele sequence reported from input.\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=CO,Number=1,Type=String,Description=\"Coordinates\">\n");
	fprintf(file, "%s", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (size_t i = 0; i < names.size(); i++) {
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", names[i].c_str());
	}
	fprintf(file, "%c", '\n');
}
double get_avglen(std::vector<Support_Node *> caller_info) {
	double size = 0;
	double num = 0;
	for (size_t i = 0; i < caller_info.size(); i++) {
		size += caller_info[i]->len;
		if (caller_info[i]->len != 0) {
			num++;
		}

	}
	return size / num;
}
int get_index_medpos(std::vector<Support_Node *> caller_info, std::pair<int, int> &cpos, std::pair<int, int> &cend) {
	std::vector<int> positions;

	for (size_t i = 0; i < caller_info.size(); i++) {
		for (size_t j = 0; j < caller_info[i]->starts.size(); j++) {
			size_t pos = 0;
			while (pos < positions.size() && positions[pos] < caller_info[i]->starts[j]) {
				pos++;
			}
			if (pos == positions.size()) {
				positions.push_back(i);
			} else {
				positions.insert(positions.begin() + j, i);
			}
		}

	}
	int index = positions[positions.size() / 2];
	for (size_t i = 0; i < caller_info.size(); i++) {

		for (size_t j = 0; j < caller_info[i]->starts.size(); j++) {
			if (i == 0 && j == 0) {

				//TODO think about it! some of these come out as negative.
				cpos.first = caller_info[i]->starts[j] - caller_info[index]->starts[0];
				cpos.second = caller_info[i]->starts[j] - caller_info[index]->starts[0];
				cend.first = caller_info[i]->stops[j] - caller_info[index]->stops[0];
				cend.second = caller_info[i]->stops[j] - caller_info[index]->stops[0];
			}
			cpos.first = std::min(caller_info[i]->starts[j] - caller_info[index]->starts[0], cpos.first);
			cpos.second = std::max(caller_info[i]->starts[j] - caller_info[index]->starts[0], cpos.second);
			cend.first = std::min(caller_info[i]->stops[j] - caller_info[index]->stops[0], cend.first);
			cend.second = std::max(caller_info[i]->stops[j] - caller_info[index]->stops[0], cend.second);
		}
	}
	return index;
}

std::string print_strands(std::pair<bool, bool> strands) {
	std::string strand = "";
	if (strands.first) {
		strand += "+";
	} else {
		strand += "-";
	}
	if (strands.second) {
		strand += "+";
	} else {
		strand += "-";
	}
	return strand;
}

void print_GTs(std::ostringstream & convert, SVS_Node * entry) {
	convert << "\tGT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO";
	int pos = 0;
	//std::cout<<"Check: "<<Parameter::Instance()->max_caller <<" vs "<<entry->caller_info.size()<<std::endl;
	for (size_t i = 0; i < Parameter::Instance()->max_caller; i++) {
		convert << "\t";
		if (pos < entry->caller_info.size() && i == entry->caller_info[pos]->id) {
			//	std::cout<<"hit: "<<i<<std::endl;
			convert << entry->caller_info[pos]->genotype;
			convert << ":";
			if (!entry->caller_info[pos]->pre_supp_vec.empty()) {
				convert << entry->caller_info[pos]->pre_supp_vec;
			} else {
				convert << "NA";
			}
			convert << ":";
			convert << entry->caller_info[pos]->len;
			convert << ":";
			convert << entry->caller_info[pos]->num_support.first;   //ref
			convert << ",";
			convert << entry->caller_info[pos]->num_support.second;   //alt
			convert << ":";
			convert << print_strands(entry->caller_info[pos]->strand);
			convert << ":";

			for (size_t j = 0; j < entry->caller_info[pos]->quality.size(); j++) { //report quality
				if (j > 0) {
					convert << ",";
				}
				if (entry->caller_info[pos]->quality[j] != -1) {
					convert << entry->caller_info[pos]->quality[j];
				} else {
					convert << ".";
				}
			}
			convert << ":";

			if (entry->caller_info[pos]->types.empty()) {
				convert << "NaN";
			} else {
				for (size_t j = 0; j < entry->caller_info[pos]->types.size(); j++) {
					if (j > 0) {
						convert << ",";
					}
					convert << trans_type(entry->caller_info[pos]->types[j]);
				}
			}
			convert << ":";   //ID VCF
			for (size_t j = 0; j < entry->caller_info[pos]->vcf_ID.size(); j++) {
				if (entry->caller_info[pos]->vcf_ID[j] != ':') {
					convert << entry->caller_info[pos]->vcf_ID[j];
				} else {
					convert << "_";
				}
			}

			if (is_valid_alleles(entry->caller_info[pos]->alleles)) {
				convert << ":";   //RAL
				convert << entry->caller_info[pos]->alleles.first;

				convert << ":";   //AAL
				convert << entry->caller_info[pos]->alleles.second;
			} else {
				convert << ":";   //RAL
				convert << "NA";

				convert << ":";   //AAL
				convert << "NA";
			}
			convert << ":";   //CALLER INFO;
			if (entry->caller_info[pos]->starts.empty()) {
				convert << "NaN";
			}
			for (size_t j = 0; j < entry->caller_info[pos]->starts.size(); j++) {
				if (j > 0) {
					convert << ",";
				}
				convert << entry->first.chr;
				convert << "_";
				if (entry->caller_info[pos]->starts[j] != 0) {
					convert << entry->caller_info[pos]->starts[j];
				} else {
					convert << "1";
				}
				convert << "-";
				convert << entry->second.chr;
				convert << "_";
				if (entry->caller_info[pos]->starts[j] != 0) {
					if (entry->caller_info[pos]->types[j] == 4) {   //TODO for mike A.
						convert << std::max(entry->caller_info[pos]->stops[j] - entry->caller_info[pos]->sv_lengths[j], 1);
					} else {
						convert << entry->caller_info[pos]->stops[j];   //entry->second.position;
					}
					//convert << entry->caller_info[pos]->stops[j];
				} else {
					convert << "1";
				}
			}
			pos++;
		} else { //check len!
				 //GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO
			convert << "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN";
		}

	}
}

void print_entry_overlap_BND(FILE *& file, SVS_Node * entry, int id) {
	if (entry->type != 3) {
		print_entry_overlap(file, entry, id);
	} else { // Just for TRA!
		std::ostringstream convert;   // stream used for the conversion
		std::pair<int, int> cipos;
		std::pair<int, int> ciend;
		int index = get_index_medpos(entry->caller_info, cipos, ciend);

		convert << entry->first.chr;   //caller_info[index]->starts[0].chr;
		convert << "\t";
		convert << entry->caller_info[index]->starts[0];  //entry->first.position;
		convert << "\t";
		convert << trans_type(entry->caller_info[index]->types[0]);
		convert << "00";
		convert << id;
		convert << "SUR\tN\t";
		if (!entry->strand.first) { //&&
			convert << "]";
			convert << entry->second.chr;
			convert << ":";
			convert << entry->second.position;
			convert << "]N";
		} else {
			convert << "N[";
			convert << entry->second.chr;
			convert << ":";
			convert << entry->second.position;
			convert << "[";
		}

		convert << "\t";
		int max_qual = -1;
		for (size_t i = 0; i < entry->caller_info.size(); i++) {
			for (size_t j = 0; j < entry->caller_info[i]->quality.size(); j++) {
				if (max_qual < entry->caller_info[i]->quality[j]) {
					max_qual = entry->caller_info[i]->quality[j];
				}
			}
		}
		if (max_qual == -1) {
			convert << '.';
		} else {
			convert << max_qual;
		}
		convert << "\t";
		convert << "PASS\t";
		convert << "SUPP=";
		convert << get_support(entry->caller_info);
		convert << ";SUPP_VEC=";
		convert << get_support_vec(entry->caller_info); //todo make aware of prev_supp/ supp vec
		convert << ";SVLEN=";
		if (entry->type != 3) {
			convert << get_avglen(entry->caller_info);
		} else {
			convert << "100000";   // TODO think about it.
		}
		/*convert << ";med_start=";
		 convert << get_start_medpos(entry->caller_info);
		 convert << ";med_stop=";
		 convert << get_stop_medpos(entry->caller_info);*/
		convert << ";SVTYPE=BND";
		convert << ";SVMETHOD=SURVIVORv2;CHR2=";
		convert << entry->second.chr;   //caller_info[index]->stops[0].chr;
		convert << ";END=";
		convert << entry->caller_info[index]->stops[0];   //entry->second.position;

		convert << ";CIPOS=";
		convert << cipos.first;
		convert << ",";
		convert << cipos.second;
		convert << ";CIEND=";
		convert << ciend.first;
		convert << ",";
		convert << ciend.second;

		//if (Parameter::Instance()->use_strand) {
		convert << ";STRANDS=";
		convert << print_strands(entry->caller_info[index]->strand);
		//}
		print_GTs(convert, entry);

		fprintf(file, "%s", convert.str().c_str());
		fprintf(file, "%c", '\n');

	}
}

void print_entry_overlap(FILE *& file, SVS_Node * entry, int id) {
	std::ostringstream convert;   // stream used for the conversion
	std::pair<int, int> cipos;
	std::pair<int, int> ciend;
	int index = get_index_medpos(entry->caller_info, cipos, ciend);

	convert << entry->first.chr;   //caller_info[index]->starts[0].chr;
	convert << "\t";
	if (entry->caller_info[index]->starts[0] != 0) {
		convert << entry->caller_info[index]->starts[0];  //entry->first.position;
	} else {
		convert << "1";
	}
	convert << "\t";
	if (is_valid_id(entry->caller_info[index]->vcf_ID)) {
		for (size_t i = 0; entry->caller_info[index]->vcf_ID[i] != ';' && i < entry->caller_info[index]->vcf_ID.size(); i++) {
			convert << entry->caller_info[index]->vcf_ID[i];
		}
	} else {
		convert << trans_type(entry->caller_info[index]->types[0]);
		convert << "00";
		convert << id;
		convert << "SUR";
	}

	convert << "\t";

	//REF / ALT field:
	if (entry->type == 3) { //FOR TRA for the BND annotation
		convert << "N\t";
		if (!entry->strand.first) { //&&
			convert << "]";
			convert << entry->second.chr;
			convert << ":";
			if (entry->second.position != 0) {
				convert << entry->second.position;
			} else {
				convert << "1";
			}
			convert << "]N";
		} else {
			convert << "N[";
			convert << entry->second.chr;
			convert << ":";
			if (entry->second.position != 0) {
				convert << entry->second.position;
			} else {
				convert << "1";
			}
			convert << "[";
		}
	} else {
		if (entry->caller_info[index]->alleles.first.empty()) {

			convert << "N";
		} else {
			convert << entry->caller_info[index]->alleles.first;
		}
		convert << "\t";

		if (entry->caller_info[index]->alleles.second.empty()) {

			convert << "<";
			convert << trans_type(entry->caller_info[index]->types[0]);
			convert << ">";

		} else {
			convert << entry->caller_info[index]->alleles.second;
		}
	}
	convert << "\t";

	int max_qual = -1;
	for (size_t i = 0; i < entry->caller_info.size(); i++) {
		for (size_t j = 0; j < entry->caller_info[i]->quality.size(); j++) {
			if (max_qual < entry->caller_info[i]->quality[j]) {
				max_qual = entry->caller_info[i]->quality[j];
			}
		}
	}
	if (max_qual == -1) {
		convert << '.';
	} else {
		convert << max_qual;
	}
	convert << "\tPASS\t";
	//INFO FIELD
	convert << "SUPP=";
	convert << get_support(entry->caller_info);
	convert << ";SUPP_VEC=";
	convert << get_support_vec(entry->caller_info); //todo make aware of prev_supp/ supp vec
	convert << ";SVLEN=";
	if (entry->type == 0) {
		convert<< (int)entry->caller_info[index]->len *-1;
	//	convert << (int) round(get_avglen(entry->caller_info)) * -1;
	} else if (entry->type != 3) {
		convert<< (int)entry->caller_info[index]->len;
		//convert << (int) round(get_avglen(entry->caller_info));
	} else {
		convert << "0";   // TODO think about it.
	}
	/*convert << ";med_start=";
	 convert << get_start_medpos(entry->caller_info);
	 convert << ";med_stop=";
	 convert << get_stop_medpos(entry->caller_info);*/
	convert << ";SVTYPE=";
	convert << trans_type(entry->caller_info[index]->types[0]);
	convert << ";SVMETHOD=SURVIVOR";
	convert << Parameter::Instance()->version;
	convert << ";CHR2=";
	convert << entry->second.chr;   //caller_info[index]->stops[0].chr;
	convert << ";END=";
	if (entry->caller_info[index]->stops[0] != 0) {
		if (entry->caller_info[index]->types[0] == 4) {   //TODO for mike A.
			convert << std::max((entry->caller_info[index]->stops[0] - (int)entry->caller_info[index]->sv_lengths[0]), 1);
		} else {
			convert << entry->caller_info[index]->stops[0];   //entry->second.position;
		}
	} else {
		convert << "1";
	}
	convert << ";CIPOS=";

	convert << cipos.first;
	convert << ",";
	convert << cipos.second;
	convert << ";CIEND=";
	convert << ciend.first;
	convert << ",";
	convert << ciend.second;

	//if (Parameter::Instance()->use_strand) {
	convert << ";STRANDS=";
	convert << print_strands(entry->caller_info[index]->strand);
	//}
	print_GTs(convert, entry);
	fprintf(file, "%s", convert.str().c_str());
	fprintf(file, "%c", '\n');
}
bool mysort(SVS_Node* i, SVS_Node* j) {

	return (i->first.position < j->first.position);

}

bool compareFunction(std::string a, std::string b) {
	return a < b;
}
bool is_not_digit(char c) {
	return !std::isdigit(c);
}

bool numeric_string_compare(const std::string& s1, const std::string& s2) {
	// handle empty strings...

	std::string::const_iterator it1 = s1.begin(), it2 = s2.begin();

	if (std::isdigit(s1[0]) && std::isdigit(s2[0])) {
		int n1, n2;
		std::stringstream ss(s1);
		ss >> n1;
		ss.clear();
		ss.str(s2);
		ss >> n2;

		if (n1 != n2)
			return n1 < n2;

		it1 = find_if(s1.begin(), s1.end(), is_not_digit);
		it2 = find_if(s2.begin(), s2.end(), is_not_digit);
	}

	return std::lexicographical_compare(it1, s1.end(), it2, s2.end());
}

std::string parse_name(char * buffer) {
	size_t i = 0;
	std::string name = "";
	while (buffer[i] != ',' && (buffer[i] != '>' && buffer[i] != '\n')) {
		name += buffer[i];
		i++;
	}
	return name;
}
void parse_vcf_header(std::map<std::string, int> &chrs, std::string filename) {

	std::string buffer;
	std::ifstream myfile;
	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cerr << "VCF Parser: could not open file: " << filename.c_str() << std::endl;
		exit(1);
	}
	getline(myfile, buffer);
	while (!myfile.eof() && buffer[0] == '#') {
		if (strncmp(&buffer[0], "##contig=", 9) == 0) {
			int count = 0;
			std::string name = "";
			for (size_t i = 12; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				//##contig=<ID=AB325691,length=20000>
				if (buffer[i] == '=' && count == 0) {
					//parse name
					name = parse_name(&buffer[i + 1]);
					if (chrs.find(name) != chrs.end()) {
						break;
					}
					count++;
				} else if (buffer[i] == '=' && count == 1) {
					//parse pos:
					chrs[name] = atoi(&buffer[i + 1]);
					break;
				}
			}
			if (chrs.find(name) == chrs.end()) {
				chrs[name] = -1;
			}
		}
		getline(myfile, buffer);
	}
}

void combine_calls_svs(std::string files, double max_dist, int min_support, int type_save, int strand_save, int dynamic_size, int min_svs, std::string output) {
	std::vector<std::string> names = parse_filename(files);

	Parameter::Instance()->max_caller = names.size();
	Parameter::Instance()->max_dist = max_dist;
	Parameter::Instance()->use_type = (type_save == 1);

	Parameter::Instance()->use_strand = (strand_save == 1);
	Parameter::Instance()->min_length = min_svs;
	Parameter::Instance()->dynamic_size = false;   //(dynamic_size==1);
	Parameter::Instance()->min_support = min_support;

//	std::cout<<"Min: "<<	Parameter::Instance()->min_support<<std::endl;
	IntervallTree bst;
	TNode *root = NULL;
	std::map<std::string, int> chrs;
	for (size_t id = 0; id < names.size(); id++) {
		parse_vcf_header(chrs, names[id]);
		std::vector<strvcfentry> entries = parse_vcf(names[id], min_svs);
		std::cout << "merging entries: " << entries.size() << std::endl;
		for (size_t j = 0; j < entries.size(); j++) {
			breakpoint_str start = convert_position(entries[j].start);
			breakpoint_str stop = convert_position(entries[j].stop);
			if (entries[j].type == 4) { //TODO for mike A.
				stop.position += entries[j].sv_len;
			}
			meta_data_str tmp;
			tmp.caller_id = (int) id;
			tmp.genotype = entries[j].genotype;
			tmp.QV = entries[j].quality;
			tmp.num_reads = entries[j].num_reads;
			tmp.sv_len = entries[j].sv_len;
			tmp.pre_supp_vec = entries[j].prev_support_vec;
			tmp.vcf_ID = entries[j].sv_id;
			tmp.allleles = entries[j].alleles;

			if(start.position==151164){
				std::cout<<"FOUND "<<std::endl;
			}
			bst.insert(start, stop, entries[j].type, entries[j].strands, tmp, root);
		}
		entries.clear();
	}

	std::map<std::string, std::vector<SVS_Node *> > union_set;
	bst.get_breakpoints(root, union_set);

	FILE * file;
	file = fopen(output.c_str(), "w");
	print_header(file, names, chrs);

	int id = 0;

	std::vector<std::string> keys;
	for (std::map<std::string, std::vector<SVS_Node *> >::iterator i = union_set.begin(); i != union_set.end(); i++) {
		keys.push_back((*i).first);
	}

//	std::cout<<"Min: "<<	Parameter::Instance()->min_support<<std::endl;
	std::sort(keys.begin(), keys.end(), numeric_string_compare);

	for (size_t i = 0; i < keys.size(); i++) {
		std::vector<SVS_Node *> points = union_set[keys[i]];
		for (std::vector<SVS_Node *>::reverse_iterator i = points.rbegin(); i != points.rend(); i++) {
			int support = get_support((*i)->caller_info);
			int len = 100000;
			if ((*i)->type != 3) {
				len = get_avglen((*i)->caller_info);
			}
			short type = (*i)->type;
			if ((*i)->type == -1) {
				type = 5;
			}

			if (support >= min_support && len > min_svs) {
				print_entry_overlap(file, (*i), id);
			}

			id++;
		}
	}
	fclose(file);
}

size_t get_len_id(int sv_len) {
	if (sv_len < 50) {
		return 0;
	} else if (sv_len < 100) {
		return 1;
	} else if (sv_len < 1000) {
		return 2;
	} else if (sv_len < 10000) {
		return 3;
	}
	return 4;
}

void summarize_VCF_files(std::string filename, int min_size, std::string output) {
	std::vector<std::string> names = parse_filename(filename);

	Parameter::Instance()->max_caller = names.size();
	Parameter::Instance()->min_length = min_size;
	Parameter::Instance()->dynamic_size = false;   //(dynamic_size==1);

	std::vector<std::vector<std::vector<int> > > svs;   //first: caller, second: type, third: size
	std::vector<std::vector<int> > caller;
	std::vector<int> tmp;
	tmp.assign(5, 0);   //sizes
	caller.assign(6, tmp);   //types
	svs.assign(names.size(), caller);

	for (size_t id = 0; id < names.size(); id++) {
		std::vector<strvcfentry> entries = parse_vcf(names[id], min_size);
		std::cout << "processing entries: " << entries.size() << std::endl;
		for (size_t j = 0; j < entries.size(); j++) {
			int len = get_len_id(entries[j].sv_len);
			if (entries[j].type == -1) {
				entries[j].type = 5;
			}
			if ((entries[j].type >= 0 && entries[j].type < 6) && (len >= 0 && len < 5)) {
				if (entries[j].type == 3) { //TRA
					len = 4;
				}
				svs[id][entries[j].type][len]++;
			}
		}
		entries.clear();
	}
	FILE * file;
	file = fopen(output.c_str(), "w");
	fprintf(file, "%s", "Caller\tDEL20\tDEL50bp\tDEL100\tDEL1k\tDEL10k\tDUP20\tDUP50bp\tDUP100\tDUP1k\tDUP10k\t");
	fprintf(file, "%s", "INV20\tINV50bp\tINV100\tINV1k\tINV10k\tTRA10k\t");
	fprintf(file, "%s", "INS20\tINS50bp\tINS100\tINS1k\tINS10k\tUNK20\tUNK50bp\tUNK100\tUNK1k\tUNK10k\n");
	for (size_t i = 0; i < svs.size(); i++) {

		fprintf(file, "%s", names[i].c_str());
		fprintf(file, "%s", "\t");
		for (size_t types = 0; types < svs[i].size(); types++) {
			if (types != 3) {
				for (size_t len = 0; len < svs[i][types].size(); len++) {
					fprintf(file, "%i", svs[i][types][len]);
					fprintf(file, "%s", "\t");
				}
			} else {
				fprintf(file, "%i", svs[i][types][4]);
				fprintf(file, "%s", "\t");
			}

		}
		fprintf(file, "%s", "\n");
	}
	fclose(file);
}
