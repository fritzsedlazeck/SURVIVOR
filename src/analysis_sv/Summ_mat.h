/*
 * Summ_mat.h
 *
 *  Created on: Jul 5, 2017
 *      Author: sedlazec
 */

#ifndef ANALYSIS_SV_SUMM_MAT_H_
#define ANALYSIS_SV_SUMM_MAT_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
using namespace std;

void summarize_svs_table_window (std::string venn_file,int window,std::string output);
void summarize_svs_table_window_stream(int window, std::string output);

#endif /* ANALYSIS_SV_SUMM_MAT_H_ */
