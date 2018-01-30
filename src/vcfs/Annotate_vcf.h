/*
 * Annotate_vcf.h
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */

#ifndef ANNOTATE_VCF_H_
#define ANNOTATE_VCF_H_
#include "../structs.h"
#include "Merge_VCF.h"
#include "../simulator/Eval_vcf.h"
struct LTR_reg{
	std::string chr;
	int start;
	int stop;
};
struct SV_reg{
	std::string header;
	std::string chr;
	strcoordinate start;
	strcoordinate stop;
	int type;
};
void generate_gene_list(std::string vcf_file, std::string annotation,int max_distance, std::string output);

int get_num_strains(strvcfentry entry);
void overlap_gtf(std::string vcf_file, std::string gtf_file,int max_distance,int min_num_occurance,int max_num_occurance,int type, std::string output);

void gene_overlap(std::string SV_file,std::string LTR_file,std::string gtf_file, int min_dist_LTR, int min_dist_gene,std::string output);
#endif /* ANNOTATE_VCF_H_ */
