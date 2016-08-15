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
	fprintf(file, "%s", "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# supporting variant reads\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=TY,Number=1,Type=String,Description=\"Types\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=CO,Number=1,Type=String,Description=\"Coordinates\">\n");

	fprintf(file, "%s", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (size_t i = 0; i < names.size(); i++) {
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", names[i].c_str());
	}
	fprintf(file, "%c", '\n');
}
double get_avglen(std::vector<Support_Node *>caller_info){
	double size=0;
	double num=0;
	for(size_t i=0;i<caller_info.size();i++){
		size+=caller_info[i]->len;
		if(caller_info[i]->len!=0){
			num++;
		}
	}
	return size/num;
}
void print_entry_overlap(FILE *& file, SVS_Node * entry, int id) {
	std::ostringstream convert;   // stream used for the conversion
	convert << entry->first.chr;
	convert << "\t";
	convert << entry->first.position;
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
	convert << get_avglen(entry->caller_info);
	convert << ";SVTYPE=";
	convert << trans_type(entry->type);
	convert << ";SVMETHOD=SURVIVORv2;CHR2=";
	convert << entry->second.chr;
	convert << ";END=";
	convert << entry->second.position;
	convert << "\tGT:LN:DV:TY:CO\t";

	for (size_t i = 0; i < entry->caller_info.size(); i++) {
		convert << "./.:";
		convert << entry->caller_info[i]->len;
		convert << ":";
		convert << entry->caller_info[i]->num_support;
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

void combine_calls_svs(std::string files, int max_dist, int min_support, int type_save, std::string output) {
	std::vector<std::string> names = parse_filename(files);

	Parameter::Instance()->max_caller = names.size();
	Parameter::Instance()->max_dist = max_dist;
	Parameter::Instance()->use_type = (type_save == 1);
	IntervallTree bst;
	TNode *root = NULL;

	for (size_t id = 0; id < names.size(); id++) {
		std::vector<strvcfentry> entries = parse_vcf(names[id]);
		std::cout << "merging entries: " << entries.size() << std::endl;
		for (size_t j = 0; j < entries.size(); j++) {
			bst.insert(convert_position(entries[j].start), convert_position(entries[j].stop), entries[j].type,entries[j].num_reads, (int) id, root);
		}
		entries.clear();
	}

	std::vector<SVS_Node *> points;
	bst.get_breakpoints(root, points);

	FILE * file;
	file = fopen(output.c_str(), "w");
	print_header(file, names);

	std::vector<int> hist;

	for (size_t i = 0; i < points.size(); i++) {
		int support=get_support(points[i]->caller_info);
		if (support >= min_support) {
			print_entry_overlap(file, points[i], i);
		}
		while(support>=hist.size()){
			hist.push_back(0);
		}
		hist[support]++;
	}
	std::cout<<"Histogram over the # of callers overlapping per SVs: "<<std::endl;
	for(size_t i=1;i<hist.size();i++){
		std::cout<<i<<" "<<hist[i]<<std::endl;
	}
	fclose(file);
	std::string out=output;
	out+="_venn";
	file = fopen(out.c_str(), "w");


	for(size_t i=0;i<names.size();i++){
		fprintf(file,"%s",names[i].c_str());
		fprintf(file,"%s","\t");
	}
	fprintf(file,"%s","\n");
	for (size_t i = 0; i < points.size(); i++) {

		for (size_t j = 0; j < points[i]->caller_info.size(); j++) {
			if(!points[i]->caller_info[j]->starts.empty()){
				fprintf(file,"%i",1);
			}else{
				fprintf(file,"%i",0);
			}
			fprintf(file,"%s","\t");
		}
		fprintf(file,"%s","\n");
	}
	fclose(file);
}
