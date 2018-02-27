/*
 * Select_samples.h
 *
 *  Created on: Feb 27, 2018
 *      Author: sedlazec
 */

#ifndef ANALYSIS_SV_SELECT_SAMPLES_H_
#define ANALYSIS_SV_SELECT_SAMPLES_H_
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
//#include"../vcfs/Merge_VCF.h"
using namespace std;


void select_greedy(std::string vcf_file,  std::string output);

#endif /* ANALYSIS_SV_SELECT_SAMPLES_H_ */
