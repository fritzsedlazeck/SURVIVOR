/*
 * Detect_nested.h
 *
 *  Created on: Apr 27, 2017
 *      Author: fsedlaze
 */

#ifndef VCFS_DETECT_NESTED_H_
#define VCFS_DETECT_NESTED_H_

#include "Merge_VCF.h"
#include "../simulator/Eval_vcf.h"

struct nested_sv{
	std::string chr;
	int id;
	int del;
	int inv;
	int dup;
	int others;
};

void detect_nested(std::string vcf_file, std::string output );

#endif /* VCFS_DETECT_NESTED_H_ */
