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
#include "../simulator/Eval_vcf.h"
#include "Paramer.h"

struct breakpoint_str {
	std::string chr;
	int position;
};


class Support_Node{
public:
	Support_Node(){
		len=0;
		num_support=0;
	}
	~Support_Node(){

	}
	int len;
	std::vector<short> types;
	std::vector<breakpoint_str> starts;
	std::vector<breakpoint_str> stops;
	int num_support;
};

class SVS_Node {
public:
	//just for testing!

	SVS_Node() {
		type=-1;
		num_support=-1;
		caller_info.clear();
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
	int num_support;

};

#include "../structs.h"
#include "IntervallTree.h"
#include "../vcfs/Merge_VCF.h"

void combine_calls_svs(std::string file, int max_dist, int min_support, int type_save, std::string output);
#endif /* MERGE_VCF_COMBINE_SVS_H_ */
