/*
 * Error_scanner.h
 *
 *  Created on: Jun 30, 2017
 *      Author: sedlazec
 */

#ifndef SIMULATOR_ERROR_SCANNER_H_
#define SIMULATOR_ERROR_SCANNER_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
//#include <random>

using namespace std;

struct differences_str {
	int position;
	short type;
};

struct CigarOp {
	char Type;   //!< CIGAR operation type (MIDNSHPX=)
	int Length; //!< CIGAR operation length (number of bases)
};

struct read_position {
	double match;
	double mismatch;
	double ins;
	double del;
	double total;
};

void generate_error_profile(int min_length, std::string output);

#endif /* SIMULATOR_ERROR_SCANNER_H_ */
