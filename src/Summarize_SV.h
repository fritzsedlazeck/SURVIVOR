/*
 * Summarize_SV.h
 *
 *  Created on: Nov 18, 2015
 *      Author: fsedlaze
 */

#ifndef SUMMARIZE_SV_H_
#define SUMMARIZE_SV_H_
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include "vcfs/Merge_VCF.h"

using namespace std;
void summary_SV(std::string filename, int min_size, int max_size, int min_reads,std::string output);
void summary_venn(std::string filename, std::string output);


#endif /* SUMMARIZE_SV_H_ */
