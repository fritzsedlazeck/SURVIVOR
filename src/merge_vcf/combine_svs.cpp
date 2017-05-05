/*
 * combine_svs.cpp
 *
 *  Created on: Jul 6, 2016
 *      Author: fsedlaze
 */

#include "combine_svs.h"

breakpoint_str convert_position(strcoordinate pos) {
	breakpoint_str tmp;
	tmp.chr = pos.chr;
	tmp.position = pos.pos;
	return tmp;
}

int get_support(std::vector<Support_Node *> caller_info) {
	int support = 0;
	for (size_t i = 0; i < caller_info.size(); i++) {
		if (!caller_info[i]->starts.empty()) {
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
void print_header(FILE *& file, std::vector<std::string> names) {
	fprintf(file, "%s", "##fileformat=VCFv4.1\n");
	std::string time = currentDateTime();
	fprintf(file, "%s", "##fileDate=");
	fprintf(file, "%s", time.c_str());
	fprintf(file, "%s", "\n"); //TODO change!
	fprintf(file, "%s", "##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(file, "%s", "##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(file, "%s", "##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(file, "%s", "##ALT=<ID=TRA,Description=\"Translocation\">\n");
	fprintf(file, "%s", "##ALT=<ID=INS,Description=\"Insertion\">\n");
	fprintf(file, "%s", "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n");
	fprintf(file, "%s", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">\n");
	fprintf(file, "%s", "##INFO=<ID=RE,Number=1,Type=Integer,Description=\"read support\">\n");
	fprintf(file, "%s", "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=LN,Number=1,Type=Integer,Description=\"predicted length\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# supporting reference,variant reads in that order\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=ST,Number=1,Type=String,Description=\"Strand of SVs\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=TY,Number=1,Type=String,Description=\"Types\">\n");
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
long get_start_medpos(std::vector<Support_Node *> caller_info) {
	std::vector<int> positions;

	for (size_t i = 0; i < caller_info.size(); i++) {
		for (size_t j = 0; j < caller_info[i]->starts.size(); j++) {
			size_t pos = 0;
			while (pos < positions.size() && positions[pos] < caller_info[i]->starts[j].position) {
				pos++;
			}
			if (pos == positions.size()) {
				positions.push_back(caller_info[i]->starts[j].position);
			} else {
				positions.insert(positions.begin() + j, caller_info[i]->starts[j].position);
			}
		}

	}
	return positions[positions.size() / 2];
}

long get_stop_medpos(std::vector<Support_Node *> caller_info) {
	std::vector<int> positions;
	for (size_t i = 0; i < caller_info.size(); i++) {
		for (size_t j = 0; j < caller_info[i]->starts.size(); j++) {
			size_t pos = 0;
			while (pos < positions.size() && positions[pos] < caller_info[i]->starts[j].position) {
				pos++;
			}
			if (pos == positions.size()) {
				positions.push_back(caller_info[i]->starts[j].position);
			} else {
				positions.insert(positions.begin() + j, caller_info[i]->starts[j].position);
			}
		}

	}
	return positions[positions.size() / 2];
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

void print_entry_overlap(FILE *& file, SVS_Node * entry, int id) {
	std::ostringstream convert;   // stream used for the conversion
	convert << entry->first.chr;
	convert << "\t";
	convert << get_start_medpos(entry->caller_info);//entry->first.position;
	convert << "\t";
	convert << trans_type(entry->type);
	convert << "00";
	convert << id;
	convert << "SUR\tN\t<";
	convert << trans_type(entry->type);
	convert << ">\t.\t";
	convert << "PASS\t";
	convert << "SUPP=";
	convert << get_support(entry->caller_info);
	convert << ";AVGLEN=";
	if (entry->type != 3) {
		convert << get_avglen(entry->caller_info);
	} else {
		convert << "100000";   // TODO think about it.
	}
	/*convert << ";med_start=";
	convert << get_start_medpos(entry->caller_info);
	convert << ";med_stop=";
	convert << get_stop_medpos(entry->caller_info);*/
	convert << ";SVTYPE=";
	convert << trans_type(entry->type);
	convert << ";SVMETHOD=SURVIVORv2;CHR2=";
	convert << entry->second.chr;
	convert << ";END=";
	convert <<  get_stop_medpos(entry->caller_info);//entry->second.position;
	//if (Parameter::Instance()->use_strand) {
	convert << ";STRANDS=";
	convert << print_strands(entry->strand);
	//}
	convert << "\tGT:LN:DR:ST:TY:CO\t";

	for (size_t i = 0; i < entry->caller_info.size(); i++) {
		convert << entry->caller_info[i]->genotype;
		convert << ":";
		convert << entry->caller_info[i]->len;
		convert << ":";
		convert << entry->caller_info[i]->num_support.first;   //ref
		convert << ",";
		convert << entry->caller_info[i]->num_support.second;   //alt
		convert << ":";
		convert << print_strands(entry->caller_info[i]->strand);
		convert << ":";
		if (entry->caller_info[i]->types.empty()) {
			convert << "NaN";
		}
		for (size_t j = 0; j < entry->caller_info[i]->types.size(); j++) {
			if (j > 0) {
				convert << ",";
			}
			convert << trans_type(entry->caller_info[i]->types[j]);
		}
		convert << ":";
		if (entry->caller_info[i]->starts.empty()) {
			convert << "NaN";
		}
		for (size_t j = 0; j < entry->caller_info[i]->starts.size(); j++) {
			if (j > 0) {
				convert << ",";
			}
			convert << entry->caller_info[i]->starts[j].chr;
			convert << "_";
			convert << entry->caller_info[i]->starts[j].position;
			convert << "-";
			convert << entry->caller_info[i]->stops[j].chr;
			convert << "_";
			convert << entry->caller_info[i]->stops[j].position;
		}
		convert << "\t";
	}
	fprintf(file, "%s", convert.str().c_str());
	fprintf(file, "%c", '\n');
}

void combine_calls_svs(std::string files, int max_dist, int min_support, int type_save, int strand_save, int dynamic_size, int min_svs, std::string output) {
	std::vector<std::string> names = parse_filename(files);

	Parameter::Instance()->max_caller = names.size();
	Parameter::Instance()->max_dist = max_dist;
	Parameter::Instance()->use_type = (type_save == 1);
	Parameter::Instance()->use_strand = (strand_save == 1);
	Parameter::Instance()->min_length = min_svs;
	Parameter::Instance()->dynamic_size = false;   //(dynamic_size==1);

	IntervallTree bst;
	TNode *root = NULL;

	for (size_t id = 0; id < names.size(); id++) {
		std::vector<strvcfentry> entries = parse_vcf(names[id], min_svs);
		std::cout << "merging entries: " << entries.size() << std::endl;
		for (size_t j = 0; j < entries.size(); j++) {
			breakpoint_str start = convert_position(entries[j].start);
			breakpoint_str stop = convert_position(entries[j].stop);
			bst.insert(start, stop, entries[j].type, entries[j].num_reads, (int) id, entries[j].genotype, entries[j].strands, entries[j].sv_len, root);

		}
		entries.clear();
	}

	std::vector<SVS_Node *> points;
	bst.get_breakpoints(root, points);

	FILE * file;
	file = fopen(output.c_str(), "w");
	print_header(file, names);

	std::vector<int> hist;
	std::vector<std::vector<std::vector<int> > > svs_summary; //first: support, second: type, third: size
	std::vector<std::vector<int> > support_vec;
	std::vector<int> tmp;
	tmp.assign(5, 0); //sizes
	support_vec.assign(6, tmp);
	svs_summary.assign(25, support_vec);

	for (size_t i = 0; i < points.size(); i++) {
		int support = get_support(points[i]->caller_info);
		int len = 100000;
		if (points[i]->type != 3) {
			len = get_avglen(points[i]->caller_info);
		}
		short type = points[i]->type;
		if (points[i]->type == -1) {
			type = 5;
		}
		while(support>=svs_summary.size()){
			svs_summary.push_back(support_vec);
		}
		if (len < 50) {
			if (type >= 0 && type < 6) {
				svs_summary[support][type][0]++;
			}
		} else if (len < 100) {
			if (type >= 0 && type < 6) {
				svs_summary[support][type][1]++;
			}
		} else if (len< 1000) {
			if (type >= 0 && type < 6) {
				svs_summary[support][type][2]++;
			}
		} else if (len < 10000) {
			if (type >= 0 && type < 6) {
				svs_summary[support][type][3]++;
			}
		} else {
			if (type >= 0 && type < 6) {
				svs_summary[support][type][4]++;
			}
		}

		if (support >= min_support) {
			print_entry_overlap(file, points[i], i);
		}
		while (support >= hist.size()) {
			hist.push_back(0);
		}
		hist[support]++;
	}
	std::cout << "Histogram over the # of callers overlapping per SVs: " << std::endl;
	for (size_t i = 1; i < hist.size(); i++) {
		std::cout << i << " " << hist[i] << std::endl;
	}
	fclose(file);
	std::string out = output;
	out += "_venn";
	file = fopen(out.c_str(), "w");

	fprintf(file, "%s", "Identifier");
	for (size_t i = 0; i < names.size(); i++) {
		fprintf(file, "%s", names[i].c_str());
		fprintf(file, "%s", "\t");
	}
	fprintf(file, "%s", "Length");
	fprintf(file, "%s", "\n");
	for (size_t i = 0; i < points.size(); i++) {
		std::ostringstream convert;
		if (points[i]->type > 5) {
			std::vector<int> determine_types;
			determine_types.assign(5, 0);
			for (size_t j = 0; j < points[i]->caller_info.size(); j++) {
				for (size_t t = 0; t < points[i]->caller_info[j]->types.size(); t++) {
					if (points[i]->caller_info[j]->types[t] < 5) {
						determine_types[points[i]->caller_info[j]->types[t]]++;
					}
				}
			}
			short id = -1;
			int max = -1;
			for (size_t j = 0; j < determine_types.size(); j++) {
				if (determine_types[j] > max) {
					max = determine_types[j];
					id = (short) j;
				}
			}
			points[i]->type = id;
		}
		convert << trans_type(points[i]->type);
		convert << "00";
		convert << i;
		convert << "SUR";
		fprintf(file, "%s", convert.str().c_str());
		fprintf(file, "%s", "\t");
		for (size_t j = 0; j < points[i]->caller_info.size(); j++) {
			if (!points[i]->caller_info[j]->starts.empty()) {
				fprintf(file, "%i", 1);
			} else {
				fprintf(file, "%i", 0);
			}
			fprintf(file, "%s", "\t");
		}
		if (points[i]->type != 3) {
			fprintf(file, "%f", get_avglen(points[i]->caller_info));
		} else {
			fprintf(file, "%i", 100000);
		}
		fprintf(file, "%s", "\n");
	}
	fclose(file);

	out = output;
	out += "_venn_summary";
	file = fopen(out.c_str(), "w");
	fprintf(file, "%s", "Support\tDEL20\tDEL50bp\tDEL100\tDEL1k\tDEL10k\tDUP20\tDUP50bp\tDUP100\tDUP1k\tDUP10k\t");
	fprintf(file, "%s", "INV20\tINV50bp\tINV100\tINV1k\tINV10k\tTRA20\tTRA50bp\tTRA100\tTRA1k\tTRA10\t");
	fprintf(file, "%s", "INS20\tINS50bp\tINS100\tINS1k\tINS10k\tUNK20\tUNK50bp\tUNK100\tUNK1k\tUNK10k\n");
	for (size_t i = 1; i < svs_summary.size(); i++) {
		fprintf(file, "%i", (int) i);
		fprintf(file, "%s", "\t");
		for (size_t types = 0; types < svs_summary[i].size(); types++) {
			for (size_t len = 0; len < svs_summary[i][types].size(); len++) {
				fprintf(file, "%i", svs_summary[i][types][len]);
				fprintf(file, "%s", "\t");
			}
		}
		fprintf(file, "%s", "\n");
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

	IntervallTree bst;
	TNode *root = NULL;

	std::vector<std::vector<std::vector<int> > > svs;   //first: caller, second: type, third: size
	std::vector<std::vector<int> > caller;
	std::vector<int> tmp;
	tmp.assign(5, 0);   //sizes
	caller.assign(6, tmp);   //types
	svs.assign(names.size(), caller);

	for (size_t id = 0; id < names.size(); id++) {
		std::vector<strvcfentry> entries = parse_vcf(names[id], min_size);
		std::cout << "merging entries: " << entries.size() << std::endl;
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
	fprintf(file, "%s", "INV20\tINV50bp\tINV100\tINV1k\tINV10k\tTRA20\tTRA50bp\tTRA100\tTRA1k\tTRA10\t");
	fprintf(file, "%s", "INS20\tINS50bp\tINS100\tINS1k\tINS10k\tUNK20\tUNK50bp\tUNK100\tUNK1k\tUNK10k\n");
	for (size_t i = 0; i < svs.size(); i++) {

		fprintf(file, "%i", (int) i);
		fprintf(file, "%s", "\t");
		for (size_t types = 0; types < svs[i].size(); types++) {
			for (size_t len = 0; len < svs[i][types].size(); len++) {
				fprintf(file, "%i", svs[i][types][len]);
				fprintf(file, "%s", "\t");
			}
		}
		fprintf(file, "%s", "\n");
	}
	fclose(file);
}
