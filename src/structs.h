/*
 * structs.h
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */

#ifndef STRUCTS_H_
#define STRUCTS_H_


#include <string>
#include <vector>
#include <map>


struct strcoordinate{
	int pos;
	std::string chr;
};

struct strvcfentry{
	std::string header;
	strcoordinate start;
	strcoordinate stop;
	short type; //0=DEL,1=DUP,2=INV,3=TRA
	std::map<std::string,std::string> calls;
	int sup_lumpy;
	int caller_id;
	std::vector<int> caller_supports;
	std::pair<bool,bool> strands;
};


struct strentry{
	int valid;
	int not_covered;
	int not_valid;
};

struct strregion{
	strcoordinate start;
	strcoordinate stop;
};

struct strsimul{
	strcoordinate start;
	strcoordinate stop;
	short type;
	bool identified;
};

struct strgene{
	int count;
	strregion region;
	std::string gene_name;
};

#endif /* STRUCTS_H_ */
