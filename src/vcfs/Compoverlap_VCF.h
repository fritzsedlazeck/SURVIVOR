/*
 * Compoverlap_VCF.h
 *
 *  Created on: Feb 27, 2015
 *      Author: fsedlaze
 */

#ifndef COMPOVERLAP_VCF_H_
#define COMPOVERLAP_VCF_H_
#include "Merge_VCF.h"
#include "../simulator/Eval_vcf.h"

void comp_overlap_vcf(std::string vcf1, std::string vcf2,int max_dis,std::string output);
void print_entry(strvcfentry entry, FILE *& out);
void print_header(std::string vcf_file, FILE *& out);

#endif /* COMPOVERLAP_VCF_H_ */
