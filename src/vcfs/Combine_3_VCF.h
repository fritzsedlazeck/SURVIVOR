/*
 * Combine_3_VCF.h
 *
 *  Created on: Mar 10, 2015
 *      Author: fsedlaze
 */

#ifndef COMBINE_3_VCF_H_
#define COMBINE_3_VCF_H_
#include "../simulator/Eval_vcf.h"
#include "Compoverlap_VCF.h"

void combine_calls(std::string vcf_delly, std::string vcf_lumpy,std::string vcf_pindel,int max_dist,std::string output);
void combine_calls_new(std::string files, int max_dist,int min_caller, std::string output);

#endif /* COMBINE_3_VCF_H_ */
