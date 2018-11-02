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
		quality=-1;
	}
	~Support_Node(){

	}
	int id;
	int len;
	int quality;
	std::vector<short> types;
	std::vector<short> sv_lengths;
	std::vector<int> starts;
	std::vector<int> stops;
	std::pair<int,int> num_support;
	std::pair<bool,bool> strand;
	std::string genotype;
	std::string pre_supp_vec;
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


};

#include "../structs.h"
#include "IntervallTree.h"
#include "../vcfs/Merge_VCF.h"

void combine_calls_svs(std::string file, int max_dist, int min_support, int type_save, int strand_save,int dynamic_size,int min_svs, std::string output);
breakpoint_str convert_position(strcoordinate pos);
void summarize_VCF_files(std::string filename, int min_size, std::string output);
void print_entry_overlap(FILE *& file, SVS_Node * entry, int id);
#endif /* MERGE_VCF_COMBINE_SVS_H_ */
