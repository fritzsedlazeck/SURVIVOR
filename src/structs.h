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
	std::pair<int,int> num_reads; //ref alt
	std::string genotype;
	int sv_len;
	std::string sv_id;
	double af;
	std::string prev_support_vec;
	int quality;
	std::pair<std::string,std::string> alleles;
	std::pair<int,int> cpos;
	std::pair<int,int> cend;
	int supp;
	//int num_reads;
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
	bool wrong;
};

struct strgene{
	int count;
	strregion region;
	std::string gene_name;
};

#endif /* STRUCTS_H_ */
