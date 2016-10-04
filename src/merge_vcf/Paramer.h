/*
 * Paramer.h
 *
 *  Created on: Aug 20, 2015
 *      Author: fsedlaze
 */

#ifndef PARAMER_H_
#define PARAMER_H_

#include <string.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <ctime>


class Parameter {
private:
	Parameter() {
	}
	~Parameter() {

	}
	static Parameter* m_pInstance;

public:
	int max_dist;
	int max_caller;
	bool use_type;
	bool use_strand;

	static Parameter* Instance() {
		if (!m_pInstance) {
			m_pInstance = new Parameter;
		}
		return m_pInstance;
	}

	double meassure_time(clock_t begin ,std::string msg){
		return 0;
		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		std::cout << msg<<" " << elapsed_secs<<std::endl;
		return elapsed_secs;
	}
};

#endif /* PARAMER_H_ */
