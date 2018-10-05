/*
 * Phasing_vcf.h
 *
 *  Created on: Sep 26, 2018
 *      Author: sedlazec
 */

#ifndef SRC_PHASING_PHASING_VCF_H_
#define SRC_PHASING_PHASING_VCF_H_
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iosfwd>
#include <algorithm>

struct snp_str{
	std::string chr;
	int position;
	int phase_block;
	bool haplotype; //true =1
	char alt_allele;
	short parental; //0=na ; 1=father; 2=mother;
};

void parental_phasing(std::string parents_vcf, std::string hapcut_output, std::string output);


#endif /* SRC_PHASING_PHASING_VCF_H_ */
