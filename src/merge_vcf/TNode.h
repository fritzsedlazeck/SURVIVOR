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

	TNode(breakpoint_str start, breakpoint_str stop, short type, std::pair<bool, bool> strands, meta_data_str meta_info) {
		this->data = new SVS_Node();
		this->data->first = start;
		this->data->second = stop;
		this->data->type = type;
		this->data->strand = strands;
		this->data->genotype = meta_info.genotype; //do I need this?

		init();
		Support_Node * tmp = new Support_Node();
		if (meta_info.sv_len == -1) {
			tmp->len = stop.position - start.position;
		} else {
			tmp->len = meta_info.sv_len;
		}
		tmp->quality.push_back(meta_info.QV);
		tmp->num_support = meta_info.num_reads;
		tmp->id = meta_info.caller_id;
		tmp->starts.push_back(start.position);
		tmp->stops.push_back(stop.position);
		tmp->types.push_back(type);
		tmp->genotype = meta_info.genotype;
		tmp->strand = strands;
		tmp->pre_supp_vec = meta_info.pre_supp_vec;
		tmp->alleles = meta_info.allleles;
		tmp->vcf_ID = meta_info.vcf_ID;
		data->caller_info.push_back(tmp);
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

	void add(breakpoint_str start, breakpoint_str stop, short type, std::pair<bool, bool> strands, meta_data_str meta_info) {
		int index = -1;
		for (size_t i = 0; i < this->data->caller_info.size(); i++) {
			if (this->data->caller_info[i]->id == meta_info.caller_id) {
				index = i;
			}
		}

		if (index == -1) {
			index = this->data->caller_info.size(); //todo check!
			Support_Node * tmp = new Support_Node();
			tmp->id = meta_info.caller_id;
			this->data->caller_info.push_back(tmp);
		}

		this->data->caller_info[index]->starts.push_back(start.position);
		this->data->caller_info[index]->stops.push_back(stop.position);
		this->data->caller_info[index]->types.push_back(type);
		this->data->caller_info[index]->num_support.first = std::max(meta_info.num_reads.first, this->data->caller_info[index]->num_support.first);
		this->data->caller_info[index]->num_support.second = std::max(meta_info.num_reads.second, this->data->caller_info[index]->num_support.second);
		this->data->caller_info[index]->genotype = meta_info.genotype;
		this->data->caller_info[index]->strand = strands;
		this->data->caller_info[index]->pre_supp_vec = meta_info.pre_supp_vec;
		this->data->caller_info[index]->quality.push_back(meta_info.QV);

		if (meta_info.allleles.first.size() > this->data->caller_info[index]->alleles.first.size() || meta_info.allleles.second.size() > this->data->caller_info[index]->alleles.second.size()) {
			this->data->caller_info[index]->alleles.first = meta_info.allleles.first;
			this->data->caller_info[index]->alleles.second = meta_info.allleles.second;
		}

		if (!this->data->caller_info[index]->vcf_ID.empty()) {
			this->data->caller_info[index]->vcf_ID +=";";// meta_info.vcf_ID;
		}
		this->data->caller_info[index]->vcf_ID = meta_info.vcf_ID;

		if (this->data->caller_info[index]->len == 0) { //first time
			this->data->caller_info[index]->len = meta_info.sv_len; //stop.position-start.position; // take the length of the svs as identifier.
		} else {
			this->data->caller_info[index]->len = std::max(meta_info.sv_len, this->data->caller_info[index]->len); //stop.position-start.position;
		}
	}
};

#endif /* TREE_TNODE_H_ */
