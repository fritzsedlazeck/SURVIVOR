/*
 * IntervallTree.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#include "IntervallTree.h"

bool IntervallTree::same_breakpoint(breakpoint_str first, breakpoint_str second){
	return (strcmp(first.chr.c_str(), second.chr.c_str()) == 0 && (abs(first.position - second.position) < Parameter::Instance()->max_dist));
}

long IntervallTree::overlap(breakpoint_str start, breakpoint_str stop,short type, SVS_Node * curr_svs) {

	/*if(Parameter::Instance()->use_type && type != curr_svs->type){ //not nice but will do
		return  (start.position - curr_svs->first.position);
	}*/
	if(!(Parameter::Instance()->use_type && type != curr_svs->type) && (same_breakpoint(start,curr_svs->first) && same_breakpoint(stop, curr_svs->second))){
		return 0; //to be merged
	}
	return (start.position - curr_svs->first.position);
}

// Inserting a node
void IntervallTree::insert(breakpoint_str start, breakpoint_str stop ,short type,int num_reads, int caller_id, TNode *&p) {

	if (p == NULL) {
		p = new TNode(start,stop,type,num_reads,caller_id);
		if (p == NULL) {
			std::cout << "Out of Space\n" << std::endl;
		}
	} else {
		long score = overlap(start,stop,type, p->get_data()); //comparison function
		if (score > 0) {
			insert(start,stop,type,num_reads,caller_id, p->left);
			if ((bsheight(p->left) - bsheight(p->right)) == 2) {
				score = overlap(start,stop,type, p->left->get_data());
				if (score > 0) {
					p = srl(p);
				} else {
					p = drl(p);
				}
			}
		} else if (score < 0) {
			insert(start,stop,type,num_reads,caller_id, p->right);
			if ((bsheight(p->right) - bsheight(p->left)) == 2) {
				score = overlap(start,stop,type, p->right->get_data());
				if (score < 0) {
					p = srr(p);
				} else {
					p = drr(p);
				}
			}
		} else { //overlaps!
			p->add(start,stop,type,num_reads,caller_id);
		}
	}
	int m, n, d;
	m = bsheight(p->left);
	n = bsheight(p->right);
	d = max(m, n);
	p->set_height(d + 1);
}
// Finding the Smallest
TNode * IntervallTree::findmin(TNode * p) {
	if (p == NULL) {
		std::cout << "The tree is empty\n" << std::endl;
		return p;
	} else {
		while (p->left != NULL) {
			p = p->left;
			//return p;
		}
		return p;
	}
}
// Finding the Largest node
TNode * IntervallTree::findmax(TNode * p) {
	if (p == NULL) {
		std::cout << "The tree is empty\n" << std::endl;
		return p;
	} else {
		while (p->right != NULL) {
			p = p->right;
			//return p;
		}
		return p;
	}
}
// Finding an get_value()
void IntervallTree::find(SVS_Node * point, TNode * &p) {
	if (p == NULL) {
		std::cout << "Sorry! get_value() not found\n" << std::endl;
	} else {
		long score = overlap(point->first,point->second,point->type,p->get_data());
		if (score > 0) {
			find(point, p->left);
		} else if (score < 0) {
			find(point, p->right);
		} else {
			std::cout << "get_value() found!\n" << std::endl;
		}

	}
}
// Copy a tree
void IntervallTree::copy(TNode * &p, TNode * &p1) {
	makeempty(p1);
	p1 = nodecopy(p);
}
// Make a tree empty
void IntervallTree::makeempty(TNode * &p) {
	TNode * d;
	if (p != NULL) {
		makeempty(p->left);
		makeempty(p->right);
		d = p;
		free(d);
		p = NULL;
	}
}
// Copy the nodes
TNode * IntervallTree::nodecopy(TNode * &p) {
	TNode * temp;
	if (p == NULL) {
		return p;
	} else {
		temp = new TNode(p->get_data()); //TODO!
		temp->left = nodecopy(p->left);
		temp->right = nodecopy(p->right);
		return temp;
	}
}

// Deleting a node
void IntervallTree::del(SVS_Node * point, TNode * &p) {
	TNode * d;
	if (p == NULL) {
		std::cout << "Sorry! get_value() not found\n" << std::endl;
	} else {
		long score = overlap(point->first, point->second,point->type,p->get_data());
		if (score > 0) {
			del(point, p->left);
		} else if (score < 0) {
			del(point, p->right);
		} else if ((p->left == NULL) && (p->right == NULL)) {
			d = p;
			free(d);
			p = NULL;
			std::cout << "get_value() deleted successfully\n" << std::endl;
		} else if (p->left == NULL) {
			d = p;
			free(d);
			p = p->right;
			std::cout << "get_value() deleted successfully\n" << std::endl;
		} else if (p->right == NULL) {
			d = p;
			p = p->left;
			free(d);
			std::cout << "get_value() deleted successfully\n" << std::endl;
		} else {
			//p->set_value(deletemin(p->right));
		}
	}
}

int IntervallTree::deletemin(TNode * &p) {
	/*	int c;
	std::cout << "inside deltemin\n" << std::endl;
	if (p->left == NULL) {

		p = p->right;
		return c;
	} else {
		c = deletemin(p->left);
		return c;
	}*/
	return 0;
}
void IntervallTree::preorder(TNode * p) {
	if (p != NULL) {
		//std::cout << p->get_data()->to_string() << "\t";
		preorder(p->left);
		preorder(p->right);
	}
}
void IntervallTree::get_breakpoints(TNode *p, std::vector<SVS_Node *> & points) {
	if (p != NULL) {
		get_breakpoints(p->left, points);
		points.push_back(p->get_data());
		get_breakpoints(p->right, points);
	}
}

// Inorder Printing
void IntervallTree::inorder(TNode * p, TNode * root) {
	if (p != NULL) {
		inorder(p->left, root);
		//std::cout << p->get_data()->to_string();
		if (p == root) {
			std::cout << "*\t";
		} else {
			std::cout << "\t";
		}
		inorder(p->right, root);
	}
}

// PostOrder Printing
void IntervallTree::postorder(TNode * p) {
	if (p != NULL) {
		postorder(p->left);
		postorder(p->right);
		//std::cout << p->get_data()->to_string() << "\t";
	}
}

int IntervallTree::max(int value1, int value2) {
	return ((value1 > value2) ? value1 : value2);
}
int IntervallTree::bsheight(TNode * p) {
	int t;
	if (p == NULL) {
		return -1;
	} else {
		t = p->get_height();
		return t;
	}
}

TNode * IntervallTree::srl(TNode * &p1) {
	TNode * p2;
	p2 = p1->left;
	p1->left = p2->right;
	p2->right = p1;
	p1->set_height(max(bsheight(p1->left), bsheight(p1->right)) + 1);
	p2->set_height(max(bsheight(p2->left), p1->get_height()) + 1);
	return p2;
}
TNode * IntervallTree::srr(TNode * &p1) {
	TNode * p2;
	p2 = p1->right;
	p1->right = p2->left;
	p2->left = p1;
	p1->set_height(max(bsheight(p1->left), bsheight(p1->right)) + 1);
	p2->set_height(max(p1->get_height(), bsheight(p2->right)) + 1);
	return p2;
}
TNode * IntervallTree::drl(TNode * &p1) {
	p1->left = srr(p1->left);
	return srl(p1);
}
TNode * IntervallTree::drr(TNode * &p1) {
	p1->right = srl(p1->right);
	return srr(p1);
}

int IntervallTree::nonodes(TNode * p) {
	int count = 0;
	if (p != NULL) {
		nonodes(p->left);
		nonodes(p->right);
		count++;
	}
	return count;
}
