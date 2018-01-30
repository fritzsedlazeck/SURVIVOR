/*
 * Overlap_snps.h
 *
 *  Created on: Jul 18, 2017
 *      Author: sedlazec
 */

#ifndef SNP_OVERLAP_OVERLAP_SNPS_H_
#define SNP_OVERLAP_OVERLAP_SNPS_H_

#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iosfwd>
#include <algorithm>
#include <sstream>
#include <algorithm>
#include "../simulator/Eval_vcf.h"
#include "../merge_vcf/Paramer.h"
using namespace std;
void overlap_snpsGWASDB(std::string svs_file, std::string snp_file, int max_dist, int min_svs, int allele, std::string output);
void overlap_snps(std::string svs_file, std::string snp_file, int max_dist, int min_svs, int allele, std::string output);
void overlap_snps_gwas(std::string svs_file, std::string random_SV,int max_dist, int min_svs, std::string output);
void generate_random_regions(std::string genome_file, std::string svs_vcf, int min_svs, std::string output);
#endif /* SNP_OVERLAP_OVERLAP_SNPS_H_ */
