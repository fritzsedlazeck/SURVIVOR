/*
 * GIAB_summary.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: fsedlaze
 */

#include "GIAB_summary.h"

std::vector<int> parse_vec(char * buffer) {
	std::vector<int> ids;
	size_t i = 0;
	while (buffer[i] != ';') {
		if (buffer[i] == '1') {
			ids.push_back(i);
		}
		i++;
	}
	return ids;
}

std::string parse_results(std::string res) {

	std::map<std::string, int> tech;
	tech["10X"] = 0;
	tech["Bio"] = 0;
	tech["CG"] = 0;
	tech["Ill"] = 0;
	tech["PB"] = 0;

	std::map<std::string, int> mendilian;
	mendilian["HG2"] = 0;
	mendilian["HG3"] = 0;
	mendilian["HG4"] = 0;

	vector<int> methods;
	methods.push_back(0); //assembly
	methods.push_back(0); //mapping

	std::map<std::string, int> method;
	method["allpass"] = 0;
	method["assemblyticsfalcon"] = 0;
	method["assemblyticsPBcR"] = 0;
	method["bkscan"] = 0;
	method["breakseq"] = 0;
	method["cnvnatorlowcov"] = 0;
	method["commonlaw"] = 0;
	method["Cortex"] = 0;
	method["delly"] = 0;
	method["FB"] = 0;
	method["fermikitraw"] = 0;
	method["fermikitsv"] = 0;
	method["GATKHC"] = 0;
	method["GATKHCSBGrefine"] = 0;
	method["hap1kb"] = 0;
	method["HySA"] = 0;
	method["Krunchall"] = 0;
	method["lumpy"] = 0;
	method["manta"] = 0;
	method["MetaSV"] = 0;
	method["MSPacMonboth"] = 0;
	method["parliament"] = 0;
	method["PB10Xdip"] = 0;
	method["PBHoneyBaylor"] = 0;
	method["pbsv"] = 0;
	method["pindelpass"] = 0;
	method["scalpel"] = 0;
	method["smrtsvdip"] = 0;
	method["snifflesngmlr"] = 0;
	method["Spiral"] = 0;
	method["SpiralSDKrefine"] = 0;
	method["svaba"] = 0;
	method["SvEvents"] = 0;
	method["SVrefine10Xhap1"] = 0;
	method["SVrefine10Xhap12"] = 0;
	method["SVrefine10Xhap2"] = 0;
	method["SVrefineDISCOVAR"] = 0;
	method["SVrefineDISCOVARDovetail"] = 0;
	method["SVrefineFalcon1"] = 0;
	method["SVrefineFalcon1Dovetail"] = 0;
	method["SVrefineFalcon2"] = 0;
	method["SVrefineFalcon2Bionano"] = 0;
	method["SVrefinePBcR"] = 0;
	method["SVrefinePBcRDovetail"] = 0;
	method["tardis"] = 0;
	method["tnscope"] = 0;
	method["vcfBeta"] = 0;

	int count = 0;
	std::string sample = "";
	std::string tech_str = "";
	std::string caller = "";
	for (size_t i = 0; i < res.size(); i++) {
		if (count == 0 && res[i] != '_') {
			sample += res[i];
		}
		if (count == 1 && res[i - 1] == '_') {
			tech_str = res[i];
			tech_str += res[i + 1];
			if (res[i + 2] != '_') {
				tech_str += res[i + 2];
			}
		}
		if (count == 2 && res[i] != ';') {
			caller += res[i];
		}

		if (res[i] == '_') {
			count++;
		}
		if (res[i] == ';') {
			mendilian[sample]++;
			tech[tech_str]++;
			method[caller]++;
			if (strcmp(caller.c_str(), "assemblyticsfalcon") == 0 || strcmp(caller.c_str(), "assemblyticsPBcR") == 0 || strcmp(caller.c_str(), "SVrefine") == 0) {
				methods[0]++;
			} else {
				methods[1]++;
			}
			sample = "";
			tech_str = "";
			caller = "";
			count = 0;
		}
	}

	std::stringstream ss;
	ss << ";TECH=";
	ss << "tot:";
	int num = 0;
	for (std::map<std::string, int>::iterator i = tech.begin(); i != tech.end(); i++) {
		if ((*i).second > 0) {
			num++;
		}
	}
	ss << num;
	for (std::map<std::string, int>::iterator i = tech.begin(); i != tech.end(); i++) {
		ss << ",";
		ss << (*i).first;
		ss << ":";
		ss << (*i).second;
	}
	ss << ";MEND=";
	ss << "tot:";
	num = 0;
	for (std::map<std::string, int>::iterator i = mendilian.begin(); i != mendilian.end(); i++) {
		if ((*i).second > 0) {
			num++;
		}
	}
	ss << num;
	for (std::map<std::string, int>::iterator i = mendilian.begin(); i != mendilian.end(); i++) {
		ss << ",";
		ss << (*i).first;
		ss << ":";
		ss << (*i).second;
	}
	ss << ";Callers=";
	ss << "tot:";
	num = 0;
	for (std::map<std::string, int>::iterator i = method.begin(); i != method.end(); i++) {
		if ((*i).second > 0) {
			num++;
		}
	}
	ss << num;
	for (std::map<std::string, int>::iterator i = method.begin(); i != method.end(); i++) {
		ss << ",";
		ss << (*i).first;
		ss << ":";
		ss << (*i).second;

	}
	ss << ";METHODS=";
	ss << "Assembly:";
	ss << methods[0];
	ss << ",Mapping:";
	ss << methods[1];

	return ss.str();

}
void summary_giab(std::string filename, std::string output) {

	//parse header to read support vector correctly
	size_t buffer_size = 200000;
	char*buffer = new char[buffer_size];
	ifstream myfile;

	myfile.open(filename.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Annotation Parser: could not open file: " << filename.c_str() << endl;
		exit(0);
	}

	FILE *file;
	file = fopen(output.c_str(), "w");

	myfile.getline(buffer, buffer_size);
	while (!myfile.eof() && buffer[1] == '#') { //avoid header
		myfile.getline(buffer, buffer_size);
	}

	int count = 0;
	std::vector<std::string> names;
	std::string name;
	for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
		if (count > 8 && buffer[i] != '\t') {
			name += buffer[i];
		}
		if (buffer[i] == '\t') {
			count++;
			if (count > 9) {
				names.push_back(name);
				name.clear();
			}
		}
	}
	names.push_back(name);
	/*for (size_t i = 0; i < names.size(); i++) {
	 std::cout << names[i] << std::endl;
	 }*/
	myfile.getline(buffer, buffer_size);
	std::string tmp = "";

	while (!myfile.eof()) {
		if (buffer[0] != '#') { //just in case:
			int count = 0;

			for (size_t i = 0; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
				if (count > 6 && strncmp("SUPP_VEC=", &buffer[i], 9) == 0) {
					vector<int> caller = parse_vec(&buffer[i + 9]);
					//cout<<"compose:"<<endl;
					for (size_t j = 0; j < caller.size(); j++) { //prepare
						tmp += names[caller[j]];
						tmp += ";";
						//	cout<< names[caller[j]]<<endl;
					}
					std::string entry = string(buffer);
					size_t found = entry.find_first_of(';');
					if (found != std::string::npos) {
						//print first part,
						fprintf(file, "%s", entry.substr(0, found).c_str());
						//print result;
						fprintf(file, "%s", parse_results(tmp).c_str());
						//cout<<parse_results(tmp)<<endl;
						//print rest;
						fprintf(file, "%s", entry.substr(found).c_str());
						fprintf(file, "%c", '\n');
						tmp = "";
						break;
					}
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
		}
		myfile.getline(buffer, buffer_size);
	}
}
void summary_giab_old(std::string filename, std::string output) {
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
		if (type == -1) {
			type = 5;
		}
		support_sum[support]++;
		if (svs[i].size < 50) {
			if (type >= 0 && type < 6) {
				svs_summary[support][type][0]++;
			}
			svs_support_size[0][support]++;
		} else if (svs[i].size < 100) {
			if (type >= 0 && type < 6) {
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

