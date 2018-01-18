/*
 * MUMmer_overlap.h
 *
 *  Created on: Dec 27, 2017
 *      Author: sedlazec
 */

#ifndef ANALYSIS_SV_MUMMER_OVERLAP_H_
#define ANALYSIS_SV_MUMMER_OVERLAP_H_

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

void overlapp_mummer(std::string vcf_SVs_file, std::string mummer_files, int max_dist, std::string output);



#endif /* ANALYSIS_SV_MUMMER_OVERLAP_H_ */
