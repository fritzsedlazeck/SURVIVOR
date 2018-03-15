/*
 * Convert_hapcut2.h
 *
 *  Created on: Mar 14, 2018
 *      Author: sedlazec
 */

#ifndef CONVERT_CONVERT_HAPCUT2_H_
#define CONVERT_CONVERT_HAPCUT2_H_


#include "../vcfs/Merge_VCF.h"
#include "../structs.h"
#include "../simulator/Eval_vcf.h"
#include <math.h>
#include <iosfwd>
using namespace std;


void process_hapcut(std::string orig_snp, std::string hapcut2, std::string output);
#endif /* CONVERT_CONVERT_HAPCUT2_H_ */
