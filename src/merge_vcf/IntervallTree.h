/*
 * IntervallTree.h
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#ifndef TREE_INTERVALLTREE_H_
#define TREE_INTERVALLTREE_H_

#include <vector>
#include <iostream>
#include "TNode.h"
#include "Paramer.h"



class IntervallTree {
private:
	int max(int, int);
	TNode * srl(TNode *&);
	TNode * drl(TNode *&);
	TNode * srr(TNode *&);
	TNode * drr(TNode *&);
	long overlap(breakpoint_str start, breakpoint_str stop,short type, std::pair<bool,bool> strands,SVS_Node * curr_svs);
	bool same_breakpoint(breakpoint_str first, breakpoint_str second,int max_dist);
	void careful_screening(breakpoint_str &start, breakpoint_str& stop ,short type,std::pair<int,int> num_reads, int caller_id, std::string genotype, std::pair<bool,bool> strands,int sv_len, TNode *p);
	long overlap_SNP(breakpoint_str start, SVS_Node * curr_svs);
public:
	void insert(breakpoint_str &start, breakpoint_str &stop ,short type,std::pair<int,int> num_reads,int caller_id, std::string genotype, std::pair<bool,bool> strands, int sv_len, TNode *&p);
	void del(SVS_Node * point, TNode *&);
	int deletemin(TNode *&);
	void find(SVS_Node * point, TNode *&);
	TNode * findmin(TNode*);
	TNode * findmax(TNode*);
	void makeempty(TNode *&);
	void copy(TNode * &, TNode *&);
	TNode * nodecopy(TNode *&);
	void preorder(TNode*);
	void inorder(TNode*,TNode * root);
	void postorder(TNode*);
	int bsheight(TNode*);
	void get_breakpoints(TNode *p,std::vector<SVS_Node *> & points);
	void get_breakpoints(TNode *p, std::map<std::string,std::vector<SVS_Node *> > & points);
	int nonodes(TNode*);
	void collapse_intervalls(TNode *&p);

	std::string findSNP(breakpoint_str &snp, TNode *&p);
};

#endif /* TREE_INTERVALLTREE_H_ */
