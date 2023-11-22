/*
 * process_Lumpy.h
 *
 *  Created on: Feb 24, 2015
 *      Author: fsedlaze
 */

#ifndef PROCESS_LUMPY_H_
#define PROCESS_LUMPY_H_

#include "../vcfs/Merge_VCF.h"
#include "../structs.h"
#include "../simulator/Eval_vcf.h"
#include <math.h>
#include <iosfwd>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

void print_header(std::string name, std::string output);
void print_entries(std::string output, std::vector<strvcfentry>& entries);
void process_Lumpy( std::string lumpy_bede, std::string output);
void trans_vcf(std::string in_vcf, std::string out_vcf);

#endif /* PROCESS_LUMPY_H_ */
