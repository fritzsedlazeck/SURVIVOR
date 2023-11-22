/*
 * combine_svs.h
 *
 *  Created on: Jul 6, 2016
 *      Author: fsedlaze
 */

#ifndef MERGE_VCF_COMBINE_SVS_H_
#define MERGE_VCF_COMBINE_SVS_H_
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iosfwd>
#include <algorithm>
#include "../simulator/Eval_vcf.h"
#include "Paramer.h"



struct breakpoint_str {
	std::string chr;
	int position;
};

struct meta_data_str{
	int caller_id;
	short type;
	std::string genotype;
	int sv_len;
	std::string pre_supp_vec;
	int QV;
	std::pair<int, int> num_reads;
	std::string vcf_ID;
	std::pair<std::string, std::string> allleles ; //first=REF; second=ALT
};


class Support_Node{
public:
	Support_Node(){
		id=0;
		len=0;
		num_support.first=0;
		num_support.second=0;
		strand.first=false;
		strand.second=false;
		genotype="./.";
		pre_supp_vec="";
	}
	~Support_Node(){

	}
	int id;
	int len;
	std::vector<int> quality;
	std::vector<short> types;
	std::vector<short> sv_lengths;
	std::vector<int> starts;
	std::vector<int> stops;
	std::pair<int,int> num_support;
	std::pair<bool,bool> strand;
	std::string genotype;
	std::string pre_supp_vec;
	std::pair<std::string,std::string> alleles;
	std::string vcf_ID;
};

class SVS_Node {
public:
	//just for testing!

	SVS_Node() {

		type=-1;
		num_support.first=-1;
		num_support.second=-1;
		strand.first=false;
		strand.second=false;
		caller_info.clear();
		genotype="./.";

		types[0]=false; //DEL
		types[1]=false; //DUP
		types[2]=false; //INV
		types[3]=false; //TRA
		types[4]=false; //UNK

		strands[0]=false; //+
		strands[1]=false; //-
		strands[2]=false; //+
		strands[3]=false; //-

	}
	~SVS_Node() {
		caller_info.clear();
	}
	//TODO change that to getter and setter!
	short type;
	breakpoint_str first;
	breakpoint_str second;
	std::vector<Support_Node *> caller_info;
	std::string entry;
	std::pair<int,int> num_support;
	std::pair<bool,bool> strand;
	std::string genotype;
	bool types[5];
	bool strands[4];
};

#include "../structs.h"
#include "IntervallTree.h"
#include "../vcfs/Merge_VCF.h"
void parse_vcf_header(std::map<std::string, int> &chrs, std::string filename);
void combine_calls_svs(std::string file, double max_dist, int min_support, int type_save, int strand_save,int dynamic_size,int min_svs, std::string output);
breakpoint_str convert_position(strcoordinate pos);
void summarize_VCF_files(std::string filename, int min_size, std::string output);
void print_entry_overlap(FILE *& file, SVS_Node * entry, int id);
#endif /* MERGE_VCF_COMBINE_SVS_H_ */
