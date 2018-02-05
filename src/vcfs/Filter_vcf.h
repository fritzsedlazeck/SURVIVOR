/*
 * Filter_vcf.h
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */

#ifndef FILTER_VCF_H_
#define FILTER_VCF_H_
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "../structs.h"
#include "Merge_VCF.h"

void filter_vcf(std::string vcf_file,std::string genomic_regions,int min_size, int max_size,int min_reads,std::string outputvcf);

void filter_vcf_sniffles(std::string vcf_file,int min_lenght, std::string outputvcf);
void summarize_paper_gaib(std::string venn_file);
#endif /* FILTER_VCF_H_ */
