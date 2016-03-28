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
#include "../vcfs/Merge_VCF.h"
#include "../simulator/Eval_vcf.h"
#include "../vcfs/Annotate_vcf.h"

void convert_vcf(std::string vcf_file, std::string output);
void convert_vcf_bede(std::string vcffile,int min_length, std::string output);

#endif /* CONVERT_VCF_TO_BED_H_ */
