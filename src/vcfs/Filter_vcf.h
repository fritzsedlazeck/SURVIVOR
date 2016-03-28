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

void filter_vcf(std::string vcf_file,std::string genomic_regions,int min_alternative_pairs, float min_alt_ref_ratio, int max_genotype,std::string outputvcf);

#endif /* FILTER_VCF_H_ */
