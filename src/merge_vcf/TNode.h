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

	TNode(breakpoint_str start, breakpoint_str stop, short type, int num_reads, int caller_id) {
		this->data = new SVS_Node();
		this->data->first = start;
		this->data->second = stop;
		this->data->type = type;
		init();
		for (int i = 0; i < Parameter::Instance()->max_caller; i++) {
			Support_Node * tmp = new Support_Node();
			tmp->len = 0;
			tmp->num_support = 0;
			tmp->types.clear();

			data->caller_info.push_back(tmp);
		}
		//this->data->caller_info[caller_id].num_reads= //TODO!
		this->data->caller_info[caller_id]->starts.push_back(start);
		this->data->caller_info[caller_id]->stops.push_back(stop);
		this->data->caller_info[caller_id]->types.push_back(type);
		this->data->caller_info[caller_id]->len = stop.position - start.position; // take the length of the svs as identifier.
		this->data->caller_info[caller_id]->num_support=num_reads;
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

	void add(breakpoint_str start, breakpoint_str stop, short type,int num_reads, int caller_id) {
		this->data->caller_info[caller_id]->starts.push_back(start);
		this->data->caller_info[caller_id]->stops.push_back(stop);
		this->data->caller_info[caller_id]->types.push_back(type);
		this->data->caller_info[caller_id]->num_support=num_reads;
		if (this->data->caller_info[caller_id]->len == 0) { //first time
			this->data->caller_info[caller_id]->len = stop.position-start.position; // take the length of the svs as identifier.
		} else {
			this->data->caller_info[caller_id]->len += stop.position-start.position;
		}
	}
};

#endif /* TREE_TNODE_H_ */
