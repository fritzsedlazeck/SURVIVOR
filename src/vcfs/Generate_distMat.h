/*
 * Generate_distMat.h
 *
 *  Created on: Jul 17, 2017
 *      Author: sedlazec
 */

#ifndef VCFS_GENERATE_DISTMAT_H_
#define VCFS_GENERATE_DISTMAT_H_
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include "Merge_VCF.h"

void generate_dist_mat(std::string svs_vcf, std::string snp_vcf, std::string weighted_file, std::string output);

#endif /* VCFS_GENERATE_DISTMAT_H_ */
