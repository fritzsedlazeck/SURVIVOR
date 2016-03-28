/*
 * Extract_Seq.h
 *
 *  Created on: Mar 18, 2015
 *      Author: fsedlaze
 */

#ifndef EXTRACT_SEQ_H_
#define EXTRACT_SEQ_H_

#include "vcfs/Merge_VCF.h"
#include "simulator/Eval_vcf.h"
void extract_breakpoint_seq(std::string vcf_file, std::string reference_file, int len,std::string output);


#endif /* EXTRACT_SEQ_H_ */
