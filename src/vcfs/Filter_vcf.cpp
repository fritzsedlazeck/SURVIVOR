/*
 * Filter_vcf.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */
#include "Filter_vcf.h"
//parses all genotype calls and summarizes them.
void translate(const char * entry, int min_alternative_pairs, float min_alt_ref_ratio, int max_genotype, strentry & result) {
//	GT:GL:GQ:FT:RC:DR:DV:RR:RV
	size_t i = 0;
	int count = 0;
	float dv = 0; //"# high-quality variant pairs"
	float dr = 0; //"# high-quality reference pairs"
	int gq = 0; //"Genotype Quality"
	int rr = 0;
	int rv = 0;
	int lumpy = 0;
	//std::cout<<entry<<std::endl;
	while (entry[i] != '\t' && entry[i] != '\n' && entry[i] != '\0') {
		if (count == 2 && entry[i - 1] == ':') {
			gq = atoi(&entry[i]);
		}
		if (count == 5 && entry[i - 1] == ':') {
			dr = atoi(&entry[i]);
		}
		if (count == 6 && entry[i - 1] == ':') {
			dv = atoi(&entry[i]);
		}
		if (count == 7 && entry[i - 1] == ':') {
			rr = atoi(&entry[i]);
		}
		if (count == 8 && entry[i - 1] == ':') {
			rv = atoi(&entry[i]);
		}
		if (count == 9 && entry[i - 1] == ':') {
			lumpy = atoi(&entry[i]);
		}
		if (entry[i] == ':') {
			count++;
		}
		i++;
	}
//	std::cout<<dv<<" "<<dv/dr<<" "<<gq<<std::endl;
	if (dv == 0 && dr == 0) {
		result.not_covered++;
	} else if (dv + rv < min_alternative_pairs) {
		result.not_valid++;
	} else if (dr > 0 && (dv + rv) / (dr + rr) < min_alt_ref_ratio) {
		result.not_valid++;
	} else if (gq > max_genotype) {
		result.not_valid++;
	} else {
		//std::cout<<"valid"<<std::endl;
		result.valid++;
	}
}

//check if SV should be kept:
bool pass_filter(strvcfentry region, std::vector<strregion> ignore_regions) {
	//TODO: compare to bed file to filter out regions!
	for (size_t i = 0; i < ignore_regions.size(); i++) {
		if (strcmp(region.start.chr.c_str(), ignore_regions[i].start.chr.c_str()) == 0) {
			if (region.start.pos > ignore_regions[i].start.pos && region.start.pos < ignore_regions[i].stop.pos) {
				return false;
			}
		}

		if (strcmp(region.stop.chr.c_str(), ignore_regions[i].start.chr.c_str()) == 0) {
			if (region.stop.pos > ignore_regions[i].start.pos && region.stop.pos < ignore_regions[i].stop.pos) {
				return false;
			}
		}

	}

	return true;
}

//parse the bed file that defines regions that should be ignored:
std::vector<strregion> parse_bed(std::string filename) {
	std::vector<strregion> ignore_regions;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "BED Parser: could not open file: " << filename.c_str() << std::endl;
		return ignore_regions;
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		int count = 0;
		strregion tmp;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				tmp.start.chr += buffer[i];
				tmp.stop.chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				tmp.start.pos = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				tmp.stop.pos = atoi(&buffer[i]);
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		ignore_regions.push_back(tmp);
		myfile.getline(buffer, buffer_size);
	}

	myfile.close();
	return ignore_regions;
}

void filter_vcf(std::string vcf_file, std::string genomic_regions,int min_size, int max_size,double min_AF,int min_reads, std::string outputvcf) {

	if(max_size!=-1){
		std::cerr<<"Warning: Max size threshold set, TRA wont be reported as their size cannot be assesst."<<std::endl;
	}
	std::vector<strregion> ignore_regions;
	if(strncmp(genomic_regions.c_str(),"NA",2)!=0){
		ignore_regions = parse_bed(genomic_regions);
	}
	//parse vcf file and write vcf file in one go:

	std::vector<std::string> names;
	std::string buffer;
	std::ifstream myfile;

	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}
	FILE *file;
	file = fopen(outputvcf.c_str(), "w");

	getline(myfile,buffer);
	int deleted = 0;
	while (!myfile.eof()) {
		if (buffer[0] == '#') { //write header info:
			fprintf(file, "%s", buffer.c_str());
			fprintf(file, "%c", '\n');
		} else {
			strvcfentry sv= parse_vcf_entry(buffer);
			int size=(min_size)+1;

			if(sv.type!=3 && sv.type!=5 && sv.type!=-1){
				size=sv.sv_len;
			}
			//std::cout<<sv.sv_len<<" "<<sv.num_reads.second<<std::endl;
			//if((size>min_size && (size< max_size || max_size==-1))){
			//	std::cout<<"size_pass: "<<size<<" "<<min_size<<" "<<max_size<<std::endl;
			//}

			if (((sv.af ==-1 || (sv.af>min_AF) ) && (ignore_regions.empty() || pass_filter(sv, ignore_regions) )) && ( (size>min_size && (size< max_size || max_size==-1)) && (sv.num_reads.second == -1 || sv.num_reads.second >min_reads)) ) {
				fprintf(file, "%s", buffer.c_str());
				fprintf(file, "%c", '\n');
			} else {
				deleted++;
			}
		}
		getline(myfile,buffer);
	}
	std::cout << "SVs ignored: " << deleted << std::endl;
	myfile.close();
	fclose(file);
}

void filter_vcf_sniffles(std::string vcf_file, int min_lenght, std::string outputvcf) {

	std::vector<std::string> names;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}
	FILE *file;
	file = fopen(outputvcf.c_str(), "w");

	myfile.getline(buffer, buffer_size);
	int deleted = 0;
	while (!myfile.eof()) {
		if (buffer[0] == '#') { //write header info:
			fprintf(file, "%s", buffer);
			fprintf(file, "%c", '\n');
		} else {
			int count = 0;
			int len = 0;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 7 && strncmp(&buffer[i], "SVLEN=NA", 8) == 0) {

					len = min_lenght * 2;
					break;
				} else if (count == 7 && strncmp(&buffer[i], "SVLEN=", 6) == 0) {
					len = atoi(&buffer[i + 6]);
					//std::cout<<"len "<<len<<std::endl;
					break;
				}

				if (buffer[i] == '\t') {
					count++;
				}
			}

			if (len >= min_lenght) {
				fprintf(file, "%s", buffer);
				fprintf(file, "%c", '\n');
			} else {
				deleted++;
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	std::cout << "WE deleted: " << deleted << std::endl;
	myfile.close();
	fclose(file);
}
struct sv_trio {
	bool son;
	bool mother;
	bool father;
};
double report_norm(double e1, double e2) {
	if (e2 == 0) {
		return 0;
	}
	return e1;//round((e1 / e2)* 10000)/100 ;
}
void summarize_paper_gaib(std::string venn_file) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(venn_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: " << venn_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	myfile.getline(buffer, buffer_size);

	std::map<std::string, std::vector<double> > trio_summary;
	std::vector<double> tmp;
	tmp.resize(8, 0);
	trio_summary["DEL"] = tmp;
	trio_summary["DUP"] = tmp;
	trio_summary["INS"] = tmp;
	trio_summary["INV"] = tmp;
	trio_summary["TRA"] = tmp;

	while (!myfile.eof()) {
		int count = 0;

		std::string key = "";
		sv_trio trio;
		trio.father = false;
		trio.mother = false;
		trio.son = false;

		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && i < 3) {
				key += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				trio.son = (bool) (buffer[i] == '1');
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				trio.father = (bool) (buffer[i] == '1');
			}
			if (count == 3 && buffer[i - 1] == '\t') {
				trio.mother = (bool) (buffer[i] == '1');
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		if (trio_summary.find(key) == trio_summary.end()) {
			trio_summary[key] = tmp;
		}
		trio_summary[key][0]++; // total
		if (trio.son) {
			trio_summary[key][1]++; // son
		}
		if (trio.father) {
			trio_summary[key][2]++; // father
		}
		if (trio.mother) {
			trio_summary[key][3]++; // mother
		}
		if (trio.son && (trio.mother && trio.father)) {
			trio_summary[key][4]++;
		}
		if (trio.son && (trio.mother || trio.father)) {
			trio_summary[key][5]++;
		}
		if (!trio.son && (trio.mother && trio.father)) {
			trio_summary[key][6]++;
		}
		if (trio.son && (!trio.mother && !trio.father)) {
			trio_summary[key][7]++;
		}

		myfile.getline(buffer, buffer_size);
	}
	std::cout << "DEL/DUP/INS/INV/TRA" << std::endl;
	std::cout << trio_summary["DEL"][0] << '/' << trio_summary["DUP"][0] << '/' << trio_summary["INS"][0] << '/' << trio_summary["INV"][0] << '/' << trio_summary["TRA"][0] << '\t';
	std::cout << trio_summary["DEL"][1] << '/' << trio_summary["DUP"][1] << '/' << trio_summary["INS"][1] << '/' << trio_summary["INV"][1] << '/' << trio_summary["TRA"][1] << '\t';
	std::cout << trio_summary["DEL"][2] << '/' << trio_summary["DUP"][2] << '/' << trio_summary["INS"][2] << '/' << trio_summary["INV"][2] << '/' << trio_summary["TRA"][2] << '\t';
	std::cout << trio_summary["DEL"][3] << '/' << trio_summary["DUP"][3] << '/' << trio_summary["INS"][3] << '/' << trio_summary["INV"][3] << '/' << trio_summary["TRA"][3] << '\t';

	std::cout << report_norm(trio_summary["DEL"][4], trio_summary["DEL"][0]) << '/' << report_norm(trio_summary["DUP"][4], trio_summary["DUP"][0]) << '/' << report_norm(trio_summary["INS"][4], trio_summary["INS"][0]) << '/' << report_norm(trio_summary["INV"][4], trio_summary["INV"][0]) << '/'<< report_norm(trio_summary["TRA"][4], trio_summary["TRA"][0]) << '\t';
	std::cout << report_norm(trio_summary["DEL"][5], trio_summary["DEL"][0]) << '/' << report_norm(trio_summary["DUP"][5], trio_summary["DUP"][0]) << '/' << report_norm(trio_summary["INS"][5], trio_summary["INS"][0]) << '/' << report_norm(trio_summary["INV"][5], trio_summary["INV"][0]) << '/'<< report_norm(trio_summary["TRA"][5], trio_summary["TRA"][0]) << '\t';
	std::cout << report_norm(trio_summary["DEL"][6], trio_summary["DEL"][0]) << '/' << report_norm(trio_summary["DUP"][6], trio_summary["DUP"][0]) << '/' << report_norm(trio_summary["INS"][6], trio_summary["INS"][0]) << '/' << report_norm(trio_summary["INV"][6], trio_summary["INV"][0]) << '/'<< report_norm(trio_summary["TRA"][6], trio_summary["TRA"][0]) << '\t';
	std::cout << report_norm(trio_summary["DEL"][7], trio_summary["DEL"][0]) << '/' << report_norm(trio_summary["DUP"][7], trio_summary["DUP"][0]) << '/' << report_norm(trio_summary["INS"][7], trio_summary["INS"][0]) << '/' << report_norm(trio_summary["INV"][7], trio_summary["INV"][0]) << '/'<< report_norm(trio_summary["TRA"][7], trio_summary["TRA"][0]) << std::endl;
}

