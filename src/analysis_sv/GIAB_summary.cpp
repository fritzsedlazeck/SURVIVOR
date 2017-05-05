/*
 * GIAB_summary.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: fsedlaze
 */

#include "GIAB_summary.h"

void summary_giab(std::string filename, std::string output) {
	size_t buffer_size = 200000;
	char*buffer = new char[buffer_size];
	ifstream myfile;

	myfile.open(filename.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Annotation Parser: could not open file: " << filename.c_str() << endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	myfile.getline(buffer, buffer_size); //skip header

	std::vector<svstruct> svs;
	while (!myfile.eof()) {
		int count = 0;
		svstruct tmp;
		tmp.support.assign(5, 0);
		tmp.type = "";
		tmp.size = 0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
			if (i < 3) {
				tmp.type += buffer[i];
			}
			if (buffer[i - 1] == '\t') {
				if (count < 18) { //illumina
					tmp.support[0] += atoi(&buffer[i]);
				} else if (count < 29) { //Pacbio
					//	cout<<buffer[i]<<",";
					tmp.support[1] += atoi(&buffer[i]);
				} else if (count < 31) {				//GC:
					//	cout<<buffer[i]<<","<<std::endl;
					tmp.support[2] += atoi(&buffer[i]);
				} else if (count < 32) { //10x

					tmp.support[3] += atoi(&buffer[i]);
				} else if (count < 33) {
					tmp.support[4] += atoi(&buffer[i]);
				} else {
					tmp.size = atof(&buffer[i]);
					//std::cout << "ERROR" << std::endl;
				}

			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		svs.push_back(tmp);
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();

	std::vector<std::vector<std::vector<int> > > svs_summary; //first: support, second: type, third: size
	std::vector<std::vector<int> > support_vec;
	std::vector<int> tmp;
	tmp.assign(5, 0); //sizes
	support_vec.assign(6, tmp);
	svs_summary.assign(6, support_vec);

	vector<int> support_sum;
	support_sum.assign(6, 0);
	map<std::string, vector<int> > svs_support_tech;

	vector<vector<int> > svs_support_size;
	svs_support_size.assign(5, support_sum);
	//svs_support_tech.assign(5,0);

	FILE *file2;
	file2 = fopen(output.c_str(), "w");
	fprintf(file2, "%s", "Type\tsize\tillumina\tPacbio\tCG\t10x\tbionano\n");
	for (size_t i = 0; i < svs.size(); i++) {
		fprintf(file2, "%s", svs[i].type.c_str());
		fprintf(file2, "%c", '\t');
		fprintf(file2, "%f", svs[i].size);

		if (svs_support_tech.find(svs[i].type) == svs_support_tech.end()) {
			svs_support_tech[svs[i].type].assign(6, 0);
		}

		int support = 0;
		for (size_t j = 0; j < svs[i].support.size(); j++) {
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%i", svs[i].support[j]);
			if (svs[i].support[j] != 0) {
				support++;
			}
		}
		int type = get_type(svs[i].type);
		if(type==-1){
			type=5;
		}
		support_sum[support]++;
		if (svs[i].size < 50) {
			if (type >= 0 && type < 6) {
				svs_summary[support][type][0]++;
			}
			svs_support_size[0][support]++;
		} else if (svs[i].size < 100) {
			if (type >= 0 && type  < 6) {
				svs_summary[support][type][1]++;
			}
			svs_support_size[1][support]++;
		} else if (svs[i].size < 1000) {

			if (type >= 0 && type < 6) {
				svs_summary[support][type][2]++;
			}
			svs_support_size[2][support]++;
		} else if (svs[i].size < 10000) {
			if (type >= 0 && type < 6) {
				svs_summary[support][type][3]++;
			}
			svs_support_size[3][support]++;
		} else {
			if (type >= 0 && type < 6) {
				svs_summary[support][type][4]++;
			}
			svs_support_size[4][support]++;
		}
		svs_support_tech[svs[i].type][support]++;
		fprintf(file2, "%c", '\n');
	}
	fclose(file2);

	std::cout << "Summary over tech support:" << std::endl;
	for (size_t i = 0; i < support_sum.size(); i++) {
		std::cout << "Tech: " << i << " " << support_sum[i] << std::endl;
	}

	std::cout << "TYPE breakout:" << std::endl;
	for (map<std::string, vector<int> >::iterator i = svs_support_tech.begin(); i != svs_support_tech.end(); i++) {
		for (size_t j = 0; j < (*i).second.size(); j++) {
			std::cout << (*i).first << " " << (*i).second[j] << std::endl;
		}
	}

	std::cout << "SIZE breakout:" << std::endl;
	for (size_t j = 0; j < svs_support_size[0].size(); j++) {

		std::cout << "20-50bp " << svs_support_size[0][j] << std::endl;

		std::cout << "50-100bp " << svs_support_size[1][j] << std::endl;

		std::cout << "100bp-1k " << svs_support_size[2][j] << std::endl;

		std::cout << "1k-10k " << svs_support_size[3][j] << std::endl;

		std::cout << "10k+ " << svs_support_size[4][j] << std::endl;

	}

	FILE * file;
	std::string out = output;
	out += "_summary";
	file = fopen(out.c_str(), "w");
	fprintf(file, "%s", "Support\tDEL20\tDEL50bp\tDEL100\tDEL1k\tDEL10k\tDUP20\tDUP50bp\tDUP100\tDUP1k\tDUP10k\t");
	fprintf(file, "%s", "INV20\tINV50bp\tINV100\tINV1k\tINV10k\tTRA20\tTRA50bp\tTRA100\tTRA1k\tTRA10\t");
	fprintf(file, "%s", "INS20\tINS50bp\tINS100\tINS1k\tINS10k\tUNK20\tUNK50bp\tUNK100\tUNK1k\tUNK10k\n");
	for (size_t i = 1; i < svs_summary.size(); i++) {
		fprintf(file, "%i", (int) i);
		fprintf(file, "%s", "\t");
		for (size_t types = 0; types < svs_summary[i].size(); types++) {
			for (size_t len = 0; len < svs_summary[i][types].size(); len++) {
				fprintf(file, "%i", svs_summary[i][types][len]);
				fprintf(file, "%s", "\t");
			}
		}
		fprintf(file, "%s", "\n");
	}
	fclose(file);
}

