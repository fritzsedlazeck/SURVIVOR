/*
 * DetectDif.h
 *
 *  Created on: Oct 30, 2015
 *      Author: fsedlaze
 */

#ifndef DETECTDIF_H_
#define DETECTDIF_H_


#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <vector>

using namespace std;
struct svs_str{
	std::string chr;
	int start;
	int stop;
	int type;
	bool joined;
};

void detect_divergence(std::string file, float precent_overlap, std::string output);


#endif /* DETECTDIF_H_ */
