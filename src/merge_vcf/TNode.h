/*
 * TNode.h
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#ifndef TREE_TNODE_H_
#define TREE_TNODE_H_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "combine_svs.h"
#include "Paramer.h"
//#include "TNode.h"
//struct svs_str;
//struct breakpoint_str;
// struct support_str;
class TNode {
private:
	SVS_Node * data;
	//int value;
	int height;
	void init() {
		this->parent = NULL;
		this->left = NULL;
		this->right = NULL;
	}
public:
	TNode * parent;
	TNode * left;
	TNode * right;
	TNode() {
		height = 0;
		init();
		this->data = NULL;
	}
	TNode(SVS_Node * point) {
		init();
		this->data = point;
		//this->data->caller_info[caller_id].num_reads= //TODO!
		height = 0;
	}

	TNode(breakpoint_str start, breakpoint_str stop, short type, std::pair<int,int> num_reads, int caller_id,std::string genotype, std::pair<bool,bool> strands,int sv_len) {
		this->data = new SVS_Node();
		this->data->first = start;
		this->data->second = stop;
		this->data->type = type;
		this->data->strand=strands;
		this->data->genotype=genotype;

		init();
		for (int i = 0; i < Parameter::Instance()->max_caller; i++) {
			Support_Node * tmp = new Support_Node();
			tmp->len = 0;
			tmp->num_support.first = 0;
			tmp->num_support.second = 0;
			tmp->types.clear();

			data->caller_info.push_back(tmp);
		}
		//this->data->caller_info[caller_id].num_reads= //TODO!
		this->data->caller_info[caller_id]->starts.push_back(start);
		this->data->caller_info[caller_id]->stops.push_back(stop);
		this->data->caller_info[caller_id]->types.push_back(type);
		this->data->caller_info[caller_id]->len = sv_len;
		if(this->data->caller_info[caller_id]->len==-1){
			this->data->caller_info[caller_id]->len = stop.position - start.position; // take the length of the svs as identifier.
		}
		this->data->caller_info[caller_id]->num_support=num_reads;
		this->data->caller_info[caller_id]->genotype=genotype;
		this->data->caller_info[caller_id]->strand=strands;
		height = 0;
	}

	~TNode() {

	}

	SVS_Node * get_data() {
		return data;
	}
	int get_height() {
		return height;
	}
	void set_height(int val) {
		this->height = val;
	}

	void add(breakpoint_str start, breakpoint_str stop, short type,std::pair<int,int> num_reads, int caller_id,std::string genotype,int svlen,std::pair<bool,bool> strands) {
		this->data->caller_info[caller_id]->starts.push_back(start);
		this->data->caller_info[caller_id]->stops.push_back(stop);
		this->data->caller_info[caller_id]->types.push_back(type);
		this->data->caller_info[caller_id]->num_support.first=std::max(num_reads.first,this->data->caller_info[caller_id]->num_support.first);
		this->data->caller_info[caller_id]->num_support.second=std::max(num_reads.second,this->data->caller_info[caller_id]->num_support.second);
	    this->data->caller_info[caller_id]->genotype=genotype;
		this->data->caller_info[caller_id]->strand=strands;
		if (this->data->caller_info[caller_id]->len == 0) { //first time
			this->data->caller_info[caller_id]->len = svlen;//stop.position-start.position; // take the length of the svs as identifier.
		} else {
			this->data->caller_info[caller_id]->len = std::max( svlen,this->data->caller_info[caller_id]->len);//stop.position-start.position;
		}
	}
};

#endif /* TREE_TNODE_H_ */
