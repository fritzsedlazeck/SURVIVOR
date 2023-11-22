/*
 * Convert_VCF_to_BED.h
 *
 *  Created on: Mar 3, 2015
 *      Author: fsedlaze
 */

#ifndef CONVERT_VCF_TO_BED_H_
#define CONVERT_VCF_TO_BED_H_

#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "../merge_vcf/combine_svs.h"
#include "../vcfs/Merge_VCF.h"
#include "../simulator/Eval_vcf.h"
#include "../vcfs/Annotate_vcf.h"

void convert_vcf(std::string vcf_file, std::string output);
void convert_vcf_bede(std::string vcffile,int min_length, std::string output);
void process_bed_file(std::string bedfile,std::string type,std::string output);
void parse_VCF_to_bed(std::string vcffile,int min_length,int max_length, std::string output);
void change_insert_pos(std::string vcffile, std::string output);
void prepare_svviz(std::string vcffile, std::string bam, std::string ref, std::string output);
#endif /* CONVERT_VCF_TO_BED_H_ */
