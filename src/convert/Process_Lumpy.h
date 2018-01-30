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

void print_header(std::string name, std::string output);
void print_entries(std::string output, std::vector<strvcfentry>& entries);
void process_Lumpy( std::string lumpy_bede, int min_number_supporting,float max_eval, std::string output);


#endif /* PROCESS_LUMPY_H_ */
