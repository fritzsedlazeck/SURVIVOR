/*
 * Extract_Seq.cpp
 *
 *  Created on: Mar 18, 2015
 *      Author: fsedlaze
 */

#include "Extract_Seq.h"
std::map<std::string, std::string> parse_ref(std::string reference_file) {
	std::map<std::string, std::string> genome;

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(reference_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: "
				<< reference_file.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	std::string seq;
	std::string name;
	while (!myfile.eof()) {
		if (buffer[0] == '>') {
			if (!seq.empty()) {
				genome[name] = seq;
				std::cout<<name<<" "<<seq.size()<<std::endl;
				name.clear();
				seq.clear();
			}
			for (size_t i = 1;
					(i < buffer_size && buffer[i] != ' ')
							&& (buffer[i] != '\0' && buffer[i] != '\n'); i++) {
				name += buffer[i];
			}
		} else {
			seq += std::string(buffer);
		}
		myfile.getline(buffer, buffer_size);
	}
	genome[name] = seq;
	name.clear();
	seq.clear();
	myfile.close();
	return genome;
}
void check_coords(int & start, int len, int max_len) {
	if (start < 0) {
		start = 0;
	}
	if (start + len > max_len) {
		start = (max_len - len) - 1;
	}
}
std::string create_name(strvcfentry entry, std::string prefix) {
	std::ostringstream tmp;
	tmp << ">";
	tmp << prefix;
	tmp << "_";
	tmp << trans_type(entry.type);
	tmp << "_";
	tmp << entry.start.chr;
	tmp << "_";
	tmp << entry.start.pos;
	tmp << "_";
	tmp << entry.stop.chr;
	tmp << "_";
	tmp << entry.stop.pos;
	return tmp.str();
}
std::string rev(std::string seq) {
	std::string tmp;
	for (std::string::reverse_iterator i = seq.rbegin(); i != seq.rend(); i++) {
		tmp += (*i);
	}
	return tmp;
}
void extract_DEL(strvcfentry entry, std::map<std::string, std::string> &ref,
		int len, FILE *&file) {
	/*deletions
	 ref ....1234....
	 alt ....    ....
	 #test for ref and alt: for short DELs (up to ~500nt)
	 #expect shorter band for the alt
	 ref	.....1234.....
	 ...>......<...
	 */
	//TODO hard coded length constrain:
	if (entry.stop.pos - entry.start.pos < 500) {
		fprintf(file, "%s", create_name(entry, "reg_1").c_str());
		fprintf(file, "%c", '\n');
		int start = entry.start.pos - len;
		check_coords(start, len, ref[entry.start.chr].size());
		fprintf(file, "%s", ref[entry.start.chr].substr(start, len).c_str());
		fprintf(file, "%c", '\n');
		fprintf(file, "%s", create_name(entry, "reg_2").c_str());
		fprintf(file, "%c", '\n');
		start = entry.stop.pos;
		check_coords(start, len, ref[entry.stop.chr].size());
		fprintf(file, "%s", ref[entry.stop.chr].substr(start, len).c_str());
		fprintf(file, "%c", '\n');
	}
}

void extract_INV(strvcfentry entry, std::map<std::string, std::string> &ref,
		int len, FILE *&file) {
	/*	ref ....1234....
	 alt ....4321....

	 #test for ref:
	 ref ....1234....
	 ...>.<......

	 #test for alt:
	 ref ....4321....
	 ...>.<......
	 */
	int region = (len / 2);

	fprintf(file, "%s", create_name(entry, "ref").c_str());
	fprintf(file, "%c", '\n');

	int start = entry.start.pos - region;
	check_coords(start, len, ref[entry.start.chr].size());
	fprintf(file, "%s", ref[entry.start.chr].substr(start, len).c_str());
	fprintf(file, "%c", '\n');
	fprintf(file, "%s", create_name(entry, "alt").c_str());
	fprintf(file, "%c", '\n');

	start = entry.start.pos - region;
	check_coords(start, region, ref[entry.start.chr].size());
	fprintf(file, "%s", ref[entry.start.chr].substr(start, region).c_str());

	start = entry.stop.pos - region;
	check_coords(start, region, ref[entry.stop.chr].size());
	fprintf(file, "%s", rev(ref[entry.stop.chr].substr(start, region)).c_str());
	fprintf(file, "%c", '\n');
}

void extract_TRA(strvcfentry entry, std::map<std::string, std::string> &ref,
		int len, FILE *&file) {
	/*
	 #test for ref:
	 ref	....1234.... [another place] ,,,,,,,,
	 ...>..<.....

	 #test for alt:
	 alt ....12,,,,,, [another place] ,,,,34....
	 ...>..<,,,,,
	 */

	int region = (len / 2);

	fprintf(file, "%s", create_name(entry, "ref").c_str());
	fprintf(file, "%c", '\n');
	int start = entry.start.pos - region;
	check_coords(start, len, ref[entry.start.chr].size());
	fprintf(file, "%s", ref[entry.start.chr].substr(start, len).c_str());
	fprintf(file, "%c", '\n');
	fprintf(file, "%s", create_name(entry, "alt").c_str());
	fprintf(file, "%c", '\n');

	start = entry.start.pos - region;
	check_coords(start, len, ref[entry.start.chr].size());
	fprintf(file, "%s", ref[entry.start.chr].substr(start, region).c_str());

	start = entry.stop.pos - region;
	check_coords(start, len, ref[entry.stop.chr].size());
	fprintf(file, "%s", ref[entry.stop.chr].substr(start, region).c_str());
	fprintf(file, "%c", '\n');
}

void extract_DUP(strvcfentry entry, std::map<std::string, std::string> &ref,
		int len, FILE *&file) {
	std::cout << entry.start.chr << " " << entry.start.pos << " "
			<< entry.stop.chr << " " << entry.stop.pos << std::endl;

	if (entry.stop.pos - entry.start.pos < 500) {
		/*
		 * #duplications
		 ref ....1234....
		 alt ....12341234....

		 #test for ref and alt: for short DUPs (up to ~500nt)
		 ref	.....1234.....
		 ...>......<...
		 */

		//TODO think about the region vs. two seq for breakpoints?
		fprintf(file, "%s", create_name(entry, "Dup_small_1").c_str());
		fprintf(file, "%c", '\n');

		int start = entry.start.pos - len;
		check_coords(start, len, ref[entry.start.chr].size());
		fprintf(file, "%s", ref[entry.start.chr].substr(start, len).c_str());
		fprintf(file, "%c", '\n');
		fprintf(file, "%s", create_name(entry, "Dup_small_2").c_str());
		fprintf(file, "%c", '\n');

		start = entry.stop.pos;
		check_coords(start, len, ref[entry.stop.chr].size());
		fprintf(file, "%s", ref[entry.stop.chr].substr(start, len).c_str());
		fprintf(file, "%c", '\n');

	} else {
		/*
		 #for for longer dups (> 500nt)
		 test for alt:
		 alt ....12341234....
		 var ......>..<......

		 #in the ref this should not make a product:
		 ref	.....1234.....
		 ......<>......
		 *
		 */
		//half into the Dup:
		int pos = entry.start.pos + ((entry.stop.pos - entry.start.pos) / 2);
		std::cout << pos << " " << len << std::endl;
		fprintf(file, "%s", create_name(entry, "Dup_large_1").c_str());
		fprintf(file, "%c", '\n');
		int start = pos - len;
		check_coords(start, len, ref[entry.start.chr].size());
		fprintf(file, "%s", ref[entry.start.chr].substr(start, len).c_str());
		fprintf(file, "%c", '\n');

		fprintf(file, "%s", create_name(entry, "Dup_large_2").c_str());
		fprintf(file, "%c", '\n');
		start = pos;
		check_coords(start, len, ref[entry.stop.chr].size());
		fprintf(file, "%s", ref[entry.stop.chr].substr(start, len).c_str());
		fprintf(file, "%c", '\n');

	}

}
void extract_breakpoint_seq(std::string vcf_file, std::string reference_file,
		int len, std::string outputfile) {

	std::vector<strvcfentry> svs = parse_vcf(vcf_file,0);
	std::map<std::string, std::string> ref = parse_ref(reference_file);
	std::cout<<"REF "<<ref.size()<<std::endl;
	FILE *file;
	file = fopen(outputfile.c_str(), "w");

	for (size_t i = 0; i < svs.size(); i++) {
		switch (svs[i].type) {
		case 0: //DEL
			std::cout << "DEL" << std::endl;
			extract_DEL(svs[i], ref, len, file);
			break;
		case 1: //DUP
			std::cout << "DUP" << std::endl;
			extract_DUP(svs[i], ref, len, file);
			break;
		case 2: //INV
			std::cout << "INV" << std::endl;
			extract_INV(svs[i], ref, len, file);
			break;
		case 3: // TRA
			std::cout << "TRA" << std::endl;
			extract_TRA(svs[i], ref, len, file);
			break;
		default:
			std::cerr << "Unkown: " << svs[i].type << std::endl;
			break;
		}
	}
}
