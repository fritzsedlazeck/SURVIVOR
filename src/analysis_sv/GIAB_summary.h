/*
 * GIAB_summary.h
 *
 *  Created on: Apr 13, 2017
 *      Author: fsedlaze
 */

#ifndef GIAB_SUMMARY_H_
#define GIAB_SUMMARY_H_
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

struct svstruct{
	std::vector<int> support;
	std::string type;
	double size;
};

void summary_giab(std::string venn_file, std::string output);


#endif /* GIAB_SUMMARY_H_ */
