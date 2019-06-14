/*
 * IntervallTree.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#include "IntervallTree.h"

bool IntervallTree::same_breakpoint(breakpoint_str first, breakpoint_str second, int max_dist) {
	//std::cout<<"test: " << first.chr.c_str()<<" "<< second.chr.c_str() <<" : "<<(abs(first.position - second.position)) << " "<< max_dist <<std::endl;
	return (strcmp(first.chr.c_str(), second.chr.c_str()) == 0 && (abs(first.position - second.position) < max_dist));
}

bool is_same_strand(std::pair<bool, bool> first, std::pair<bool, bool> second) {
	return (first.first == second.first && first.second == second.second);
}
bool same_type(short first, short second) {
	if (first == second) {
		return true;
	} else if (((first == 3 || first == 2) && second == 5) || ((second == 3 || second == 2) && first == 5)) { // compare BND to inv or tra -> same type!
		return true;
	}
	return false;
}
int get_dist(long dist, short type) {

	if (type == 3) { //TRA
		return Parameter::Instance()->max_dist; //TODO: change!
	}

	dist = std::max((long) Parameter::Instance()->min_length * 2, dist);

	return std::min((int) (dist * 4), Parameter::Instance()->max_dist);
}

long IntervallTree::overlap_SNP(breakpoint_str start, SVS_Node * curr_svs) {

//	std::cout << "OVL: " << start.chr << " " << start.position << "  SV: " << curr_svs->first.chr << " " << curr_svs->first.position << std::endl;
	//type:
	if (curr_svs->type == 0 || (curr_svs->type == 1 || curr_svs->type == 6)) { // DEL, DUP, CNV?
		if (strcmp(curr_svs->first.chr.c_str(), start.chr.c_str()) == 0) {
			if (start.position >= curr_svs->first.position && start.position <= curr_svs->second.position) {
				std::cout << "FOUND" << std::endl;
				return 0; // is overlapping
			}
		}
	}

	if (strcmp(curr_svs->first.chr.c_str(), start.chr.c_str()) == 0) {
		if (abs(start.position - curr_svs->first.position) <= Parameter::Instance()->max_dist) {
			//	std::cout << "FOUND" << std::endl;
			return 0; // in close proximity to start
		}
	}

	if (strcmp(curr_svs->second.chr.c_str(), start.chr.c_str()) == 0) {
		if (abs(start.position - curr_svs->second.position) <= Parameter::Instance()->max_dist) {
			//		std::cout << "FOUND" << std::endl;
			return 0; // in close proximity to stop
		}
	}
	int dist = (start.position - curr_svs->first.position);
	if (dist == 0) {
		return 1;
	}
	return (dist);
}

long IntervallTree::overlap(breakpoint_str start, breakpoint_str stop, short type, std::pair<bool, bool> strands, SVS_Node * curr_svs) {

	/*if(Parameter::Instance()->use_type && type != curr_svs->type){ //not nice but will do
	 return  (start.position - curr_svs->first.position);
	 }*/

	/*	if (start.position == 112179238 || start.position == 112179329) {
		std::cout << "Comp: " << start.chr << " " << start.position << " " << curr_svs->first.chr << " " << curr_svs->first.position << " " << type << " " << curr_svs->type << std::endl;
	}*/

	int max_dist = Parameter::Instance()->max_dist;
	/*	if (Parameter::Instance()->dynamic_size) {
	 max_dist = std::min(get_dist(stop.position - start.position, type), get_dist(curr_svs->first.position - curr_svs->second.position, curr_svs->type)); //
	 }
	 */
	if (((!Parameter::Instance()->use_strand || is_same_strand(strands, curr_svs->strand)) && (!Parameter::Instance()->use_type || same_type(type, curr_svs->type))) && (same_breakpoint(start, curr_svs->first, max_dist) && same_breakpoint(stop, curr_svs->second, max_dist))) {
	/*	if (start.position == 112179238 || start.position == 112179329) {
			std::cout << "MERGE" << std::endl;
			std::cout << std::endl;
		}*/
		return 0; //to be merged
	}

	/*if (start.position == 112179238 || start.position == 112179329) {
		std::cout << "no MERGE: "<<stop.chr<<" "<<stop.position<<" "<<curr_svs->second.chr<<" "<<curr_svs->second.position  << std::endl;
		if (is_same_strand(strands, curr_svs->strand)) {
			std::cout << "\tsame strand" << std::endl;
		}
		if (same_type(type, curr_svs->type)) {
			std::cout << "\tsame type" << std::endl;
		}
		std::cout << "dist: " << max_dist << std::endl;
		if (same_breakpoint(start, curr_svs->first, max_dist)) {
			std::cout << "\tsame breakpoint start" << std::endl;
		}
		if (same_breakpoint(stop, curr_svs->second, max_dist)) {
			std::cout << "\tsame breakpoint stop" << std::endl;
		}
		std::cout << std::endl;
	}*/
	/*	std::cout<<"no MERGE: "<<start.position <<" "<<curr_svs->first.position <<" vs "<< stop.position << " " << curr_svs->second.position<<std::endl;
	 if( is_same_strand(strands, curr_svs->strand)){
	 std::cout<<"\tsame strand"<<std::endl;
	 }
	 if(same_type(type, curr_svs->type)){
	 std::cout<<"\tsame type"<<std::endl;
	 }
	 std::cout<<"dist: "<<max_dist<<std::endl;
	 if(same_breakpoint(start, curr_svs->first, max_dist)){
	 std::cout<<"\tsame breakpoint start"<<std::endl;
	 }
	 if(same_breakpoint(stop, curr_svs->second, max_dist)){
	 std::cout<<"\tsame breakpoint stop"<<std::endl;
	 }*/
	if (strcmp(start.chr.c_str(), curr_svs->first.chr.c_str()) == 0 && abs(start.position - curr_svs->first.position) < max_dist) {

		return (stop.position - curr_svs->second.position);
	}
	int dist = (start.position - curr_svs->first.position);
	if (dist == 0) {
		return 1;
	}
	return (dist);
}

// Inserting a node SURVIVOR!

void IntervallTree::careful_screening(breakpoint_str &start, breakpoint_str& stop, short type, std::pair<bool, bool> strands, meta_data_str meta_info, TNode *&p) { //maybe I just need the pointer not a ref.
	if (p != NULL && !(start.position == -1 && stop.position == -1)) {
		careful_screening(start, stop, type, strands, meta_info, p->left);
		if (overlap(start, stop, type, strands, p->get_data()) == 0) { //SV type
			p->add(start, stop, type, strands, meta_info);
			start.position = -1;
			stop.position = -1;
			return;
		}
		careful_screening(start, stop, type, strands, meta_info, p->right);
	}
}

void IntervallTree::insert(breakpoint_str &start, breakpoint_str& stop, short type, std::pair<bool, bool> strands, meta_data_str meta_info, TNode *&p) {
	if (start.position == -1 && stop.position == -1) {
		return;
	}
	if (p == NULL) {
		p = new TNode(start, stop, type, strands, meta_info);
		if (p == NULL) {
			std::cout << "Out of Space\n" << std::endl;
		}
	} else {
		long score = overlap(start, stop, type, strands, p->get_data()); //comparison function
		if (score == 0) {
			p->add(start, stop, type, strands, meta_info);
			start.position = -1;
			stop.position = -1;
			return;
		} else if (std::abs(score) < (long) Parameter::Instance()->max_dist) { // if two or more events are too close:
			//std::cout<<"Screen"<<std::endl;
			careful_screening(start, stop, type, strands, meta_info, p);
			if (start.position == -1 && stop.position == -1) {
				return;
			}
		}
		if (score > 0) {
			insert(start, stop, type, strands, meta_info, p->left);
			if ((bsheight(p->left) - bsheight(p->right)) == 2) {
				score = overlap(start, stop, type, strands, p->left->get_data());
				if (score > 0) {
					p = srl(p);
				} else {
					p = drl(p);
				}
			}
		} else if (score < 0) {
			insert(start, stop, type, strands, meta_info, p->right);
			if ((bsheight(p->right) - bsheight(p->left)) == 2) {
				score = overlap(start, stop, type, strands, p->right->get_data());
				if (score < 0) {
					p = srr(p);
				} else {
					p = drr(p);
				}
			}
		}
	}
	int m, n, d;
	m = bsheight(p->left);
	n = bsheight(p->right);
	d = max(m, n);
	p->set_height(d + 1);
}

std::string IntervallTree::findSNP(breakpoint_str &snp, TNode *&p) {
	if (p == NULL) {
		return "NA";
		//	std::cout << "Sorry! get_value() not found\n" << std::endl;
	} else {
		long score = overlap_SNP(snp, p->get_data());
		if (score > 0) {
			return findSNP(snp, p->left);
		} else if (score < 0) {
			return findSNP(snp, p->right);
		} else {
			std::stringstream ss;
			switch (p->get_data()->type) {
			case 0:
				ss << "DEL";
				break;
			case 1:
				ss << "DUP";
				break;
			case 2:
				ss << "INV";
				break;
			case 3:
				ss << "TRA";
				break;
			case 4:
				ss << "INS";
				break;
			case 5:
				ss << "UNK";
				break;
			case 6:
				ss << "CNV";
				break;

			default:
				break;
			}
			return ss.str();
		}
	}

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
		long score = overlap(point->first, point->second, point->type, point->strand, p->get_data());
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
		long score = overlap(point->first, point->second, point->type, point->strand, p->get_data());
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

void IntervallTree::get_breakpoints(TNode *p, std::map<std::string, std::vector<SVS_Node *> > & points) {
	if (p != NULL) {
		get_breakpoints(p->left, points);
		if ((int) p->get_data()->caller_info.size() >= Parameter::Instance()->min_support) {
			points[p->get_data()->first.chr].push_back(p->get_data());
		}
		get_breakpoints(p->right, points);
	}
}

// Inorder Printing
void IntervallTree::inorder(TNode * p, TNode * root) {
	if (p != NULL) {
		inorder(p->left, root);
		//	std::cout<<p->get_data()->first.chr<<" "<<p->get_data()->first.position <<" - "<<p->get_data()->second.position <<std::endl;
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
