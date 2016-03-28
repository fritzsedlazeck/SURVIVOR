/*
 * CorrectAllele.h
 *
 *  Created on: Jun 18, 2015
 *      Author: fsedlaze
 */

#ifndef CORRECTALLELE_H_
#define CORRECTALLELE_H_

#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include "structs.h"
#include "vcfs/Merge_VCF.h"

void correct_alleles(std::string vcf_file,std::string table, std::string output);

#endif /* CORRECTALLELE_H_ */
