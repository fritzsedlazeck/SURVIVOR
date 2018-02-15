/*
 * process_Lumpy.cpp
 *
 *  Created on: Feb 24, 2015
 *      Author: fsedlaze
 */

#include "Process_Lumpy.h"

short get_type(const char * type) {
	//std::cout<<type<<std::endl;
	if (strncmp(type, "DELETION", 8) == 0) {
		return 0;
	} else if (strncmp(type, "DUPLICATION", 11) == 0) {
		return 1;
	} else if (strncmp(type, "INVERSION", 8) == 0) {
		return 2;
	} else if (strncmp(type, "INTERCHROM", 9) == 0) {
		return 3;
	} else {
		std::cerr << "Unknown type!" << std::endl;
	}
	return -1;
}

int get_support(const char * line) {
//++,11;--,12
	size_t i = 0;
	int support = 0;
	while (line[i] != '\t') {
		if (line[i - 1] == ',') {
			support += atoi(&line[i]);
		}
		i++;
	}
	return support;
}

strregion get_coords(const char * buffer) {
//MAX:III:213395;III:227239
	size_t i = 0;
	int count = 1;
	strregion region;
	bool flag = false;
	while (buffer[i] != '\t') {

		if (count == 1 && buffer[i] != ':') {
			region.start.chr += buffer[i];
		}
		if (count == 2 && buffer[i - 1] == ':') {
			region.start.pos = atoi(&buffer[i]);
		}

		if (flag && buffer[i] != ':') {
			region.stop.chr += buffer[i];
		}
		if (count == 3 && buffer[i - 1] == ':') {
			region.stop.pos = atoi(&buffer[i]);
		}
		if (buffer[i] == ';') {
			flag = true;
		}
		if (buffer[i] == ':') {
			flag = false;
			count++;
		}
		i++;
	}
	//std::cout<<region.start.chr<<" "<<region.stop.chr<<std::endl;
	return region;
}
bool equal_region(strcoordinate c1, strcoordinate c2, int max_dist) {
	return (bool) (strcmp(c1.chr.c_str(), c2.chr.c_str()) == 0 && abs(c1.pos - c2.pos) < max_dist);
}

int get_entry(strregion region, short type, std::vector<strvcfentry> & entries, int max_dist) {
	for (size_t i = 0; i < entries.size(); i++) {

		if (entries[i].type == type) {
			//start - start
			//std::cout<<"\t"<<entries[i].start.chr<<" "<<entries[i].start.pos<<" "<<entries[i].stop.chr<<" "<<entries[i].stop.pos<<" "<<entries[i].type<<std::endl;
			if (equal_region(entries[i].start, region.start, max_dist) && equal_region(entries[i].stop, region.stop, max_dist)) {
				return i;
			} else if (equal_region(entries[i].start, region.stop, max_dist) && equal_region(entries[i].stop, region.start, max_dist)) {
				return i;
			}
		}
	}
	return -1;
}
strvcfentry create_entry(strregion region, double eval, int support, short type, int id) {
	strvcfentry tmp;
//	III     5104    DEL00000002     N       <DEL>   .       LowQual IMPRECISE;CIEND=-305,305;CIPOS=-305,305;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.5.9;CHR2=III;END=15991;SVLEN=10887;CT=3to5;PE=2;MAPQ=60        GT:GL:GQ:FT:RC:DR:DV:RR:RV      1/1:-12,-0.602059,0:6:LowQual:816:0:2:0:0
	tmp.start = region.start;
	tmp.stop = region.stop;
	tmp.type = type;
	tmp.sup_lumpy = support;

	std::ostringstream convert;   // stream used for the conversion
	convert << region.start.chr;
	convert << "\t";
	convert << region.start.pos;      // insert the textual representation of 'Number' in the characters in the stream
	convert << "\t";
	convert << trans_type(type);
	convert << "00";
	convert << id;
	convert << "LUM\tN\t<";
	convert << trans_type(type);
	if (tmp.sup_lumpy < 4) {
		convert << ">\t.\tLowQual\tIMPRECISE;SVTYPE=";
	} else {
		convert << ">\t.\tPASS\tIMPRECISE;SVTYPE=";
	}
	convert << trans_type(type);
	convert << ";SVMETHOD=LUMPYv0.2.9;CHR2=";
	convert << region.stop.chr;
	convert << ";END=";
	convert << region.stop.pos;
	convert << ";EVAL=";
	convert << eval;

	if (tmp.type == 3) {
		convert << ";SVLEN=0;PE=";
	} else {
		convert << ";SVLEN=";
		convert << region.stop.pos - region.start.pos;
		convert << ";PE=";
	}
	convert << support;
	convert << "\tGT:GL:GQ:FT:RC:DR:DV:RR:RV\t";
	tmp.header = convert.str();
	std::stringstream s;
	s << "1/1:0,0,0:0:PASS:0:0:";
	s << tmp.sup_lumpy;
	s << ":0:0";
	tmp.calls["lumpy"] = s.str();
	return tmp;
}
void parse_lumpy(std::string lumpy_bede, std::vector<strvcfentry> & entries, int min_number_supporting, double max_eval) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(lumpy_bede.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Lumpy Parser: could not open file: " << lumpy_bede.c_str() << std::endl;
		exit(0);
	}
	int call_id = entries.size();
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		int count = 0;
		short type = -1;
		double eval = 99;

		strregion region;
		int support = 0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 7 && buffer[i - 1] == '\t') {
				eval = atof(&buffer[i]);
			}

			if (count > 9 && strncmp(&buffer[i], "TYPE:", 5) == 0) {
				//get type;
				type = get_type(&buffer[i + 5]);
			}

			if (count > 10 && strncmp(&buffer[i], "STRANDS", 7) == 0) {
				//get support val;
				support = get_support(&buffer[i + 7]);
			}

			if (count > 10 && strncmp(&buffer[i], "MAX:", 4) == 0) {
				//get positions;
				region = get_coords(&buffer[i + 4]);
			}

			if (buffer[i] == '\t') {
				count++;
			}
		}

		//filter the parsed SV:
		if (support > min_number_supporting && eval < max_eval) {
			//std::cout<<eval<<" "<<type<<" "<<region.start.pos<<" "<<region.stop.pos<<std::endl;
			//detect overlap with delly:
			//if no overlap construct vcf entry for Lumpy:
			entries.push_back(create_entry(region, eval, support, type, call_id));
			call_id++;
		}
		myfile.getline(buffer, buffer_size);
	}

	myfile.close();
}

void print_header(std::string name, std::string output) {
	FILE *file;
	file = fopen(output.c_str(), "w");

	fprintf(file, "%s", "##fileformat=VCFv4.1\n");
	fprintf(file, "%s", "##fileDate=20150217\n");
	fprintf(file, "%s", "##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(file, "%s", "##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(file, "%s", "##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(file, "%s", "##ALT=<ID=TRA,Description=\"Translocation\">\n");
	fprintf(file, "%s", "##ALT=<ID=INS,Description=\"Insertion\">\n");

	fprintf(file, "%s", "##FILTER=<ID=LowQual,Description=\"PE support below 3.\">\n");
	fprintf(file, "%s", "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n");
	fprintf(file, "%s", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n");

	fprintf(file, "%s", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Normalized high-quality read count for the SV\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">\n");
	fprintf(file, "%s", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
	fprintf(file, "%s", name.c_str());
	fprintf(file, "%c", '\n');

	fclose(file);
}

void print_entries(std::string output, std::vector<strvcfentry>& entries) {
	FILE *file;
	file = fopen(output.c_str(), "a");
	for (size_t i = 0; i < entries.size(); i++) {
		fprintf(file, "%s", entries[i].header.c_str());

		for (std::map<std::string, std::string>::iterator j = entries[i].calls.begin(); j != entries[i].calls.end(); j++) {
			fprintf(file, "%s", (*j).second.c_str());
		}
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}

void process_Lumpy(std::string lumpy_bede, int min_number_supporting, float max_eval, std::string output) {
	std::vector<strvcfentry> entries;			//= parse_vcf(delly_vcf); //get delly calls

	parse_lumpy(lumpy_bede, entries, min_number_supporting, pow(10, max_eval));

	print_header(lumpy_bede, output);
	print_entries(output, entries);

}

std::string parse_line(std::string buffer) {

	int count = 0;
	std::stringstream ss;
	strvcfentry entry = parse_vcf_entry(buffer);
	for (size_t i = 0; i < buffer.size() && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
		if (count == 4 &&entry.type==3) {//tra
			if (buffer[i - 1] == '\t') {
				if (!entry.strands.first && entry.strands.second) {
					ss << "]";
					ss << entry.stop.chr.c_str();
					ss << ':';
					ss << entry.stop.pos;
					ss << "]N";
				} else {
					ss << "N[";
					ss << entry.stop.chr.c_str();
					ss << ':';
					ss << entry.stop.pos;
					ss << '[';
				}
				ss << '\t';
			}
		} else if (count == 7) {
			if (buffer[i - 1] == '\t') {
				if(entry.type==3) {
					ss << "SVLEN=100000;SVTYPE=BND;SVMETHOD=SURVIVORv2;CIPOS=-500,500;CIEND=-500,500;STRANDS=";
				}else{
					ss << "SVLEN=";
					ss<< entry.sv_len;
					ss<<";SVTYPE=";
					ss<<trans_type(entry.type);
					ss<<";SVMETHOD=SURVIVORv2;CIPOS=-500,500;CIEND=-500,500;STRANDS=";
				}
				if (entry.strands.first) {
					ss << "+";
				} else {
					ss << "-";
				}
				if (entry.strands.second) {
					ss << "+";
				} else {
					ss << "-";
				}
				ss << ";\t";
			}
		} else {
			ss << buffer[i];
		}

		if (buffer[i] == '\t') {
			count++;
		}
	}
	return ss.str();

}
void trans_vcf(std::string in_vcf, std::string out_vcf) {

	FILE *file;
	file = fopen(out_vcf.c_str(), "w");

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(in_vcf.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << in_vcf.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if (buffer[0] == '#') {
			fprintf(file, "%s", buffer);
			fprintf(file, "%c", '\n');
		} else {
			fprintf(file, "%s",parse_line(std::string(buffer)).c_str());
			fprintf(file, "%c", '\n');
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
}
