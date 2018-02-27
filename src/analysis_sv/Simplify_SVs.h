/*
 * Simplify_SVs.h
 *
 *  Created on: Nov 28, 2017
 *      Author: sedlazec
 */

#ifndef ANALYSIS_SV_SIMPLIFY_SVS_H_
#define ANALYSIS_SV_SIMPLIFY_SVS_H_
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
#include"../vcfs/Merge_VCF.h"
using namespace std;

struct sv_simple_str {
	strcoordinate start;
	strcoordinate stop;
	std::string svtype;
	std::vector< std::string > accessions;
	pair<bool,bool> strands;
};

void simplify_svs(std::string file, std::string pop_file, int min_size, std::string output);


#endif /* ANALYSIS_SV_SIMPLIFY_SVS_H_ */
