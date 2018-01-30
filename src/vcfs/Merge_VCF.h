/*
 * Merge_VCF.h
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */

#ifndef MERGE_VCF_H_
#define MERGE_VCF_H_
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "../merge_vcf/Paramer.h"
#include "../merge_vcf/combine_svs.h"

#include "../structs.h"
std::vector<strvcfentry> parse_vcf(std::string filename,int min_svs);
strvcfentry parse_vcf_entry(std::string buffer);
strcoordinate parse_stop(const char * buffer);
void merge_vcf(std::string filenames, int max_dist, int min_observed, std::string outputfile);
int overlap(strvcfentry tmp, std::vector<strvcfentry> & final_vcf,int max_dist);
strcoordinate parse_stop(const char * buffer);
std::pair <bool,bool>parse_strands(const char * buffer);
std::vector<std::string> parse_filename(std::string filename);
short get_type(std::string type);

#endif /* MERGE_VCF_H_ */
