/*
 * Process_Coverage.h
 *
 *  Created on: Apr 13, 2017
 *      Author: fsedlaze
 */

#ifndef CONVERT_PROCESS_COVERAGE_H_
#define CONVERT_PROCESS_COVERAGE_H_
#include "../vcfs/Merge_VCF.h"
#include "../structs.h"
#include "../simulator/Eval_vcf.h"
#include <math.h>
#include <iosfwd>
void summarize_badcoverage(std::string filename,int win_size,int min_cov, std::string output);


#endif /* CONVERT_PROCESS_COVERAGE_H_ */
