/*
 * SV_Simulator.h
 *
 *  Created on: Jan 30, 2016
 *      Author: fsedlaze
 */

#ifndef SIMULATOR_SV_SIMULATOR_H_
#define SIMULATOR_SV_SIMULATOR_H_

#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>

struct parameter {
	int dup_min;
	int dup_max;
	int dup_max_amp;
	int dup_num;

	int indel_min;
	int indel_max;
	int indel_num;

	int translocations_min;
	int translocations_max;
	int translocations_num;

	int inv_min;
	int inv_max;
	int inv_num;

	int inv_del_min;
	int inv_del_max;
	int inv_del_num;

	int inv_dup_min;
	int inv_dup_max;
	int inv_dup_num;

	int intrachr_num;
	int intrachr_min;
	int intrachr_max;

	bool diploid;
	float hom_rate;
};

struct position {
	std::string chr;
	int start;
	int stop;
};

struct struct_var {
	int type; //0:dup;1:del;2:ins;3:inv;4:tra
	position pos;
	position target;
	std::string seq; //not mandadory!
	std::string ref;
	bool print;
	int copy_num;
};

struct insertions {
	position target;
	std::string seq;
};

void simulate_SV(std::string ref_file, std::string parameter_file,float snp_freq, bool coordinates, std::string output_prefix);
void generate_parameter_file(std::string parameter_file);
#endif /* SIMULATOR_SV_SIMULATOR_H_ */
