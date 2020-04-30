/*
 * Nanopore_sim.h
 *
 *  Created on: May 30, 2017
 *      Author: sedlazec
 */

#ifndef NANOPORE_SIM_H_
#define NANOPORE_SIM_H_
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
#include <stdio.h>
#include <math.h>
#include <random>

#include "Error_scanner.h"
using namespace std;
void simulate_reads(std::string genome,std::string error_profile,int coverage, std::string output);

#endif /* NANOPORE_SIM_H_ */
