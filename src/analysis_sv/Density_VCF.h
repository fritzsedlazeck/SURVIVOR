/*
 * Density_VCF.h
 *
 *  Created on: May 15, 2020
 *      Author: sedlazec
 */

#ifndef ANALYSIS_SV_DENSITY_VCF_H_
#define ANALYSIS_SV_DENSITY_VCF_H_

#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iosfwd>
#include <algorithm>
#include "../simulator/Eval_vcf.h"
#include "../merge_vcf/Paramer.h"
#include "../merge_vcf/combine_svs.h"
using namespace std;

void density_VCF(std::string vcf_file, int window, std::string output);


#endif /* ANALYSIS_SV_DENSITY_VCF_H_ */
