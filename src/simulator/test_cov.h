/*
 * test_cov.h
 *
 *  Created on: Jun 17, 2016
 *      Author: fsedlaze
 */

#ifndef TEST_COV_H_
#define TEST_COV_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
//#include <random>

using namespace std;

///void est_cov(int coverage,int read_length,int num_SV,int genome,int min_overlap,int min_support);
void est_cov(int read_length, int num_SV, int min_overlap, int min_support,int cov);
void count_valid_reads(double allowed_n_ratio);


#endif /* TEST_COV_H_ */
