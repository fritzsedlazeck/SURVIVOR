/*
 * Detect_nested.cpp
 *
 *  Created on: Apr 27, 2017
 *      Author: fsedlaze
 */

#include "Detect_nested.h"
void detect_nested(std::string vcf_file, std::string output) {
	std::vector<strsimul> simulated;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "BED Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);

	//UTURN
	//INVDEL
	//INVDUP
	nested_sv tmp;
	tmp.del = 0;
	tmp.dup = 0;
	tmp.id = -1;
	tmp.inv = 0;
	tmp.others=0;
	tmp.chr="";

	std::vector<nested_sv> nested_stuff;
	int invdups = 0;
	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			std::string type;
			int id;
			std::string chr;
			bool flag = false;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if(count==0 && buffer[i]!='\t'){
					chr+=buffer[i];
				}
				if (count == 2 && buffer[i - 1] == '\t') {
					id = atoi(&buffer[i]);
				}
				if (count == 2 && buffer[i] != '\t') {
					if (buffer[i] == '_') {
						flag = true;
					}
				}

				if (count == 4 && buffer[i] != '\t') {
					type += buffer[i];
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			if (flag) {
				//std::cout<<"HIT "<<id<<" "<<type<<std::endl;
				size_t i = 0;
				while (i < nested_stuff.size()) {
					if (id == nested_stuff[i].id) {
						break;
					}
					i++;
				}
				if (i == nested_stuff.size()) {
					tmp.id = id;
					nested_stuff.push_back(tmp);
				}
				if (strcmp(type.c_str(), "<DEL>") == 0) {
					nested_stuff[i].del++;
				} else if (strcmp(type.c_str(), "<INV>") == 0) {
					nested_stuff[i].inv++;
				} else if (strcmp(type.c_str(), "<DUP>") == 0) {
					nested_stuff[i].dup++;
				}else{
					nested_stuff[i].others++;
				}
				if(nested_stuff[i].chr.empty()){
					nested_stuff[i].chr=chr;
				}else if(strcmp(nested_stuff[i].chr.c_str(),chr.c_str())!=0){
					nested_stuff[i].others=100;
				}
			}
			//UTURN
			if (strcmp(type.c_str(), "<INVDUP>") == 0) {
				invdups++;
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	int invdel=0;
	int invdup=0;
	for(size_t i=0;i < nested_stuff.size();i++){
		if(nested_stuff[i].others==0 && ( nested_stuff[i].del > 1 &&nested_stuff[i].inv >0 && nested_stuff[i].dup==0) ){
			invdel++;
			std::cout<<"invdel ID: "<<nested_stuff[i].id<<std::endl;
		}else if(nested_stuff[i].others==0 && ( nested_stuff[i].del == 0 &&nested_stuff[i].inv >0 && nested_stuff[i].dup>0) ){
			invdup++;
			std::cout<<"invdup ID: "<<nested_stuff[i].id<<std::endl;
		}
	}
	std::cout<<"Found invdel: "<<invdel<<std::endl;
	std::cout<<"Found invdup: "<<invdup<<std::endl;
	std::cout<<"Found uturn: "<<invdups<<std::endl;

}
