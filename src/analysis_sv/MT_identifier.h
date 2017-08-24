/*
 * MT_identifier.h
 *
 *  Created on: Aug 15, 2017
 *      Author: sedlazec
 */

#ifndef ANALYSIS_SV_MT_IDENTIFIER_H_
#define ANALYSIS_SV_MT_IDENTIFIER_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
using namespace std;

struct segment_str{
	int pos;
	std::string chr;
	bool strand;
	int MQ;
	int read_start;
};

#endif /* ANALYSIS_SV_MT_IDENTIFIER_H_ */
