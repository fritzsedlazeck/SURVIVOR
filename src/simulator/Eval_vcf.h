/*
 * Eval_vcf.h
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */

#ifndef EVAL_VCF_H_
#define EVAL_VCF_H_
#include <sstream>

#include "../vcfs/Merge_VCF.h"
#include "../structs.h"
struct strreport{
	short del;
	short dup;
	short inv;
	short tra;
	short ins;
	short other;
};
void eval_vcf(std::string vcf_file,std::string bed_file,int max_allowed_dist,std::string output);
bool match_coords(strsimul c1, strvcfentry c2, int max_allowed_dist);
std::string trans_type(short type);

void eval_paper(std::string vcf_file,std::string bed_file,int max_allowed_dist);

#endif /* EVAL_VCF_H_ */
