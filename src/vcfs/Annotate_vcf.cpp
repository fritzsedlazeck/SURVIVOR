/*
 * Annotate_vcf.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */
#include "Annotate_vcf.h"
std::string parse_genename(const char * buffer) {
	size_t i = 5;
	std::string name;
	while (buffer[i] != ';' && buffer[i] != ':') {
		name += buffer[i];
		i++;
	}
	return name;
}

std::vector<strgene> parse_annotation(std::string gtf_file) {
	std::vector<strgene> genes;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(gtf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "GTF Parser: could not open file: " << gtf_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		int count = 0;
		strgene tmp;
		tmp.count = 0;
		tmp.region.start.pos = -1; //used as indication that line was not parsed;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				tmp.region.start.chr += buffer[i];
				tmp.region.stop.chr += buffer[i];
			}

			if (count == 2 && buffer[i - 1] == '\t') {
				if (strncmp(&buffer[i], "gene", 4) != 0) { // gene???
					break;
				}
			}
			if (count == 3 && buffer[i - 1] == '\t') {
				tmp.region.start.pos = atoi(&buffer[i]);
			}

			if (count == 4 && buffer[i - 1] == '\t') {
				tmp.region.stop.pos = atoi(&buffer[i]);
			}
			if (count == 8 && strncmp(&buffer[i], "gene:", 5) == 0) { //gene:
				tmp.gene_name = parse_genename(&buffer[i]);
				break;
			}

			if (buffer[i] == '\t') {
				count++;
			}
		}
		if (tmp.region.start.pos != -1 && strcmp(tmp.region.start.chr.c_str(), "MT") != 0) {
			genes.push_back(tmp);
		}
		myfile.getline(buffer, buffer_size);
	}

	myfile.close();
	return genes;
}

std::vector<strgene> parse_annotation_bed(std::string gtf_file) {
	std::vector<strgene> genes;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(gtf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "BED Parser: could not open file: " << gtf_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		int count = 0;
		strgene tmp;
		tmp.count = 0;
		tmp.region.start.pos = -1; //used as indication that line was not parsed;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				tmp.region.start.chr += buffer[i];
				tmp.region.stop.chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				tmp.region.start.pos = atoi(&buffer[i]);
			}

			if (count == 2 && buffer[i - 1] == '\t') {
				tmp.region.stop.pos = atoi(&buffer[i]);
			}
			if (count == 3 && buffer[i] != '\t') { //gene:
				tmp.gene_name += buffer[i];
			}

			if (buffer[i] == '\t') {
				count++;
			}
		}
		if (tmp.region.start.pos != -1 && strcmp(tmp.region.start.chr.c_str(), "MT") != 0) {
			genes.push_back(tmp);
		}
		myfile.getline(buffer, buffer_size);
	}

	myfile.close();
	return genes;
}

bool is_within(strvcfentry entry, strgene gene, int max_distance) {
	//if coordinate of SV is within the gene:
	if (strcmp(entry.start.chr.c_str(), gene.region.start.chr.c_str()) == 0) {
		return (bool) (entry.start.pos > gene.region.start.pos - max_distance && entry.start.pos < gene.region.stop.pos + max_distance);
	}
	if (strcmp(entry.stop.chr.c_str(), gene.region.start.chr.c_str()) == 0) {
		return (bool) (entry.stop.pos > gene.region.start.pos - max_distance && entry.stop.pos < gene.region.stop.pos + max_distance);
	}
	return false;
}

bool is_overlapping(strvcfentry entry, strgene gene, int max_distance) {
	if (is_within(entry, gene, max_distance)) {
		return true;
	}

	//start and stop chr should be on the same chr
	if (strcmp(entry.start.chr.c_str(), gene.region.start.chr.c_str()) == 0) {
		return (bool) (entry.start.pos - max_distance < gene.region.start.pos && entry.stop.pos + max_distance > gene.region.stop.pos);
	}

	return false;
}

int get_num_strains(strvcfentry entry) { //run over the calls and count number of LowQual and PASS
	//parse call
	int num_support = 0;

	for (std::map<std::string, std::string>::iterator i = entry.calls.begin(); i != entry.calls.end(); i++) {
		//it should always be just one entry;
		const char *line = (*i).second.c_str();

		for (size_t pos = 0; pos < (*i).second.size(); pos++) {
			if (strncmp("LowQual", &line[pos], 7) == 0 || strncmp("PASS", &line[pos], 4) == 0) {
				num_support++;
			}
		}
	}
	return num_support;
}
void print_result(std::string genename, FILE *& file) {
	fprintf(file, "%s", genename.c_str());

	fprintf(file, "%c", '\n');
}
void print_result(std::string genename, int type, FILE *& file) {
	fprintf(file, "%s", genename.c_str());
	fprintf(file, "%c", '\t');
	fprintf(file, "%s", trans_type(type).c_str());
	fprintf(file, "%c", '\n');
}
void comp_overlap(std::vector<strvcfentry> entries, std::vector<strgene> genes, int max_distance, int min_num_occurance, int max_num_occurance, int type, std::string outputfile) {

	FILE *file;
	file = fopen(outputfile.c_str(), "w");
	std::map<std::string, bool> reported_genes;
	for (size_t i = 0; i < entries.size(); i++) {

		if (type == -1 || entries[i].type == type) {
			int support = get_num_strains(entries[i]);
			if ((min_num_occurance == -1 || support > min_num_occurance) && (max_num_occurance == -1 || support < max_num_occurance)) {
				if (entries[i].type == 0 || entries[i].type == 1) { // DUP and DEL (whole gene is influenced)

					for (size_t j = 0; j < genes.size(); j++) {
						if (is_overlapping(entries[i], genes[j], max_distance)) {
							if (reported_genes.find(genes[j].gene_name) == reported_genes.end()) {
								reported_genes[genes[j].gene_name] = true;
								print_result(genes[j].gene_name, file);
							}
						}
					}
				} else {
					for (size_t j = 0; j < genes.size(); j++) {
						if (is_within(entries[i], genes[j], max_distance)) {
							if (reported_genes.find(genes[j].gene_name) == reported_genes.end()) {
								reported_genes[genes[j].gene_name] = true;
								print_result(genes[j].gene_name, file);
							}
						}
					}
				}
			}
		}
	}
	reported_genes.clear();
	fclose(file);
}

void overlap_gtf(std::string vcf_file, std::string gtf_file, int max_distance, int min_num_occurance, int max_num_occurance, int type, std::string output) {
	std::vector<strvcfentry> entries = parse_vcf(vcf_file);
	std::vector<strgene> genes = parse_annotation_bed(gtf_file);
	std::cout << "parsed genes: " << genes.size() << std::endl;
	//std::cout<<"gene[0] "<< genes[0].region.start<<std::endl;
	comp_overlap(entries, genes, max_distance, min_num_occurance, max_num_occurance, type, output);
}
/////////////////////////////////////////////////////////////////////////////// new ///////////////////////////////////////////////////////////////
void change(std::string & chr) {
	//std::cout<<"CHAGE: "<<chr<<std::endl;
	if (strcmp(chr.c_str(), "chr1") == 0) {
		chr = "III";
	} else if (strcmp(chr.c_str(), "chr2") == 0) {
		chr = "II";
	} else if (strcmp(chr.c_str(), "chr3") == 0) {
		chr = "I";
	} else {
		std::cout << "ERROR in Annotation!" << std::endl;
	}
}
bool is_within(SV_reg entry, strgene gene, int max_distance) {
	//if coordinate of SV is within the gene:
	if (strcmp(entry.start.chr.c_str(), gene.region.start.chr.c_str()) == 0) {
		return (bool) (entry.start.pos > gene.region.start.pos - max_distance && entry.start.pos < gene.region.stop.pos + max_distance);
	}
	if (strcmp(entry.stop.chr.c_str(), gene.region.start.chr.c_str()) == 0) {
		return (bool) (entry.stop.pos > gene.region.start.pos - max_distance && entry.stop.pos < gene.region.stop.pos + max_distance);
	}
	return false;
}

bool is_overlapping(SV_reg entry, strgene gene, int max_distance) {
	if (is_within(entry, gene, max_distance)) {
		return true;
	}

	//start and stop chr should be on the same chr
	if (strcmp(entry.start.chr.c_str(), gene.region.start.chr.c_str()) == 0) {
		return (bool) (entry.start.pos - max_distance < gene.region.start.pos && entry.stop.pos + max_distance > gene.region.stop.pos);
	}

	return false;
}

std::map<std::string, std::vector<LTR_reg> > parse_LTR(std::string LTR_file) {
	std::map<std::string, std::vector<LTR_reg> > ltrs;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(LTR_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "GTF Parser: could not open file: " << LTR_file.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	int count = 0;
	std::string name;
	std::vector<std::string> names;
	for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
		if (count > 2 && buffer[i] != '\t') {
			name += buffer[i];
		}
		if (buffer[i] == '\t' && !name.empty()) {
			names.push_back(name);
			name.clear();
		}
		if (buffer[i] == '\t') {
			count++;
		}
	}
	if (!name.empty()) {
		names.push_back(name);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		count = 0;
		LTR_reg tmp;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				tmp.chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				change(tmp.chr);
				tmp.start = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				tmp.stop = atoi(&buffer[i]);
			}

			if (count > 2 && buffer[i - 1] == '\t') {
				if (atoi(&buffer[i]) == 1) {
					ltrs[names[count - 3]].push_back(tmp);
				}
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	return ltrs;
}

bool is_called(const char * buffer) {
	//std::cout << buffer[0] << buffer[1] << buffer[2] << std::endl;
	return (strncmp(buffer, "1/1", 3) == 0 || strncmp(buffer, "0/1", 3) == 0);
}
bool min_distance(int coord, int start, int stop, int dist) {
	return ((abs(coord - start) < dist || abs(coord - stop) < dist) || (start > coord && coord < stop));
}
bool is_in(SV_reg tmp, std::vector<LTR_reg> ltrs, int min_dist) {

	for (size_t i = 0; i < ltrs.size(); i++) {
		if (strcmp(tmp.start.chr.c_str(), ltrs[i].chr.c_str()) == 0) {

			if (min_distance(tmp.start.pos, ltrs[i].start, ltrs[i].stop, min_dist)) {
				std::cout << "HiT" << std::endl;
				return true;
			}
		}
		if (strcmp(tmp.stop.chr.c_str(), ltrs[i].chr.c_str()) == 0) {
			std::cout << "stop CHR" << std::endl;
			if (min_distance(tmp.stop.pos, ltrs[i].start, ltrs[i].stop, min_dist)) {
				std::cout << "HiT" << std::endl;
				return true;
			}
		}
	}
	return false;
}
std::map<std::string, std::vector<SV_reg> > parse_SV(std::string vcf_file, std::map<std::string, std::vector<LTR_reg> > ltrs, int min_dist_LTR) {
	std::map<std::string, std::vector<SV_reg> > sv_entries;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "GTF Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	while (!myfile.eof() && buffer[1] == '#') { //avoid header
		myfile.getline(buffer, buffer_size);
	}

	//parse names:
	int count = 0;
	std::string name;
	std::vector<std::string> names;
	for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
		if (count > 8 && buffer[i] != '\t') {
			name += buffer[i];
		}
		if (buffer[i] == '\t' && !name.empty()) {
			//std::cout<<name<<std::endl;
			names.push_back(name);
			name.clear();
		}
		if (buffer[i] == '\t') {
			count++;
		}
	}
	if (!name.empty()) {
		names.push_back(name);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		count = 0;
		SV_reg tmp;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				tmp.start.chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				tmp.start.pos = atoi(&buffer[i]);
			}
			if (count == 4 && buffer[i - 1] == '<') {
				tmp.type = get_type(std::string(&buffer[i]));
			}
			if (count == 7 && strncmp(&buffer[i], "CHR2=", 5) == 0) {
				tmp.stop = parse_stop(&buffer[i]);
			}
			if (count < 8) {
				tmp.header += buffer[i];
			}
			if (count > 8 && buffer[i - 1] == '\t') {
				if (is_called(&buffer[i]) && (ltrs.find(names[count - 9]) != ltrs.end() && !is_in(tmp, ltrs[names[count - 9]], min_dist_LTR))) {
					std::cout << "hit" << std::endl;
					sv_entries[names[count - 9]].push_back(tmp);
				}
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		myfile.getline(buffer, buffer_size);
	}

	myfile.close();
	return sv_entries;
}

void gene_overlap(std::string vcf_file, std::string LTR_file, std::string gtf_file, int min_dist_LTR, int min_dist_gene, std::string output) {

	std::map<std::string, std::vector<LTR_reg> > ltrs = parse_LTR(LTR_file);
	std::cout << "LTR:" << ltrs.size() << std::endl;
	std::map<std::string, std::vector<SV_reg> > sv_entries = parse_SV(vcf_file, ltrs, min_dist_LTR);
	std::cout << "SV:" << sv_entries.size() << std::endl;
	ltrs.clear();
	std::vector<strgene> genes = parse_annotation(gtf_file);
	for (std::map<std::string, std::vector<SV_reg> >::iterator i = sv_entries.begin(); i != sv_entries.end(); i++) {

		std::vector<SV_reg> entries = (*i).second;

		FILE *file;
		std::string out = output;
		out += (*i).first;
		file = fopen(out.c_str(), "w");

		std::map<std::string, bool> reported_genes;
		for (size_t pos = 0; pos < entries.size(); pos++) {
			if (entries[pos].type == 0 || entries[pos].type == 1) { // DUP and DEL (whole gene is influenced)
				for (size_t j = 0; j < genes.size(); j++) {
					if (is_overlapping(entries[pos], genes[j], min_dist_gene)) {
						if (reported_genes.find(genes[j].gene_name) == reported_genes.end()) { //to avoid printing the same name twice
							genes[j].count++;
							reported_genes[genes[j].gene_name] = true;
							print_result(genes[j].gene_name, entries[pos].type, file);
						}
					}
				}
			} else {
				for (size_t j = 0; j < genes.size(); j++) {
					if (is_within(entries[pos], genes[j], min_dist_gene)) {
						if (reported_genes.find(genes[j].gene_name) == reported_genes.end()) { //to avoid printing the same name twice
							genes[j].count++;
							reported_genes[genes[j].gene_name] = true;
							print_result(genes[j].gene_name, entries[pos].type, file);
						}
					}
				}
			}
		}
		fclose(file);
		reported_genes.clear();
	}
	FILE *file;
	std::string out = output;
	out += "summary";
	file = fopen(out.c_str(), "w");
	for (size_t j = 0; j < genes.size(); j++) {
		fprintf(file, "%s", genes[j].gene_name.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", genes[j].count);
		fprintf(file, "%c", '\n');
	}
	fclose(file);

}

struct sv_str {
	std::string name;
	short type;
	int start;
	int stop;
	std::string start_chr;
	std::string stop_chr;
};

std::string parse_chr(char * buffer) {
	size_t i = 0;
	std::string chr;
	while (buffer[i] != ';' && buffer[i] != '\t') {
		chr += buffer[i];
		i++;
	}
	return chr;
}

std::string trans(int type) {
	switch (type) {
	case 0:
		return "DEL";
		break;
	case 1:
		return "DUP";
		break;
	case 2:
		return "TRA";
		break;
	case 3:
		return "INV";
		break;
	}
	return "NA";
}
short get_type(char * buffer) {
	size_t i = 0;
	std::string type;
	while (buffer[i] != '>' && buffer[i] != '\t') {
		type += buffer[i];
		i++;
	}
	if (strcmp(type.c_str(), "DEL") == 0) {
		return 0;
	} else if (strcmp(type.c_str(), "DUP") == 0) {
		return 1;
	} else if (strcmp(type.c_str(), "TRA") == 0) {
		return 2;
	} else if (strcmp(type.c_str(), "INV") == 0) {
		return 3;
	}
	return -1;
}

int is_overlapping(std::string chr, int start, int stop, std::vector<sv_str> svs_elements, int max_distance) {

	for (size_t i = 0; i < svs_elements.size(); i++) {
		//breakpoint is within the genes:
		if (strcmp(svs_elements[i].start_chr.c_str(), chr.c_str()) == 0) {
			if (svs_elements[i].start > start - max_distance && svs_elements[i].start < stop + max_distance) {
				return i;
			}
		} else if (strcmp(svs_elements[i].stop_chr.c_str(), chr.c_str()) == 0) {
			if (svs_elements[i].stop > start - max_distance && svs_elements[i].stop < stop + max_distance) {
				return i;
			}
		}
		//gene is within the SV:
		if (svs_elements[i].type < 2 && strcmp(svs_elements[i].start_chr.c_str(), chr.c_str()) == 0) { //for del and dup!
			if (svs_elements[i].start - max_distance < start && svs_elements[i].stop + max_distance > stop) {
				return i;
			}
		}
	}
	return -1;
}
void generate_gene_list(std::string vcf_file, std::string annotation, int max_distance, std::string output) {

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "GTF Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}
	std::vector<sv_str> svs_elements;
	myfile.getline(buffer, buffer_size);
	//parse vcf:
	while (!myfile.eof()) { //avoid header
		if (buffer[0] != '#') {
			int count = 0;
			sv_str svs;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 0 && buffer[i] != '\t') {
					svs.start_chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					svs.start = atoi(&buffer[i]);
				}
				if (count == 4 && buffer[i - 1] == '<') {
					svs.type = get_type(&buffer[i]);
				}
				if (count == 7 && strncmp(&buffer[i], ";END=", 5) == 0) {
					svs.stop = atoi(&buffer[i + 5]);
				}
				if (count == 7 && strncmp(&buffer[i], ";CHR2=", 6) == 0) {
					svs.stop_chr = parse_chr(&buffer[i + 6]);
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			if (svs.stop_chr.empty() && (svs.type < 2 || svs.type == 3)) {
				svs.stop_chr = svs.start_chr;
			}
			std::cout << "SV: " << svs.start_chr << " " << svs.start << " " << svs.stop_chr << " " << svs.stop << std::endl;
			svs_elements.push_back(svs);
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();

//	std::cout<<"SV found "<<svs_elements.size()<<std::endl;
	//genes:
	myfile.open(annotation.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "GTF Parser: could not open file: " << annotation.c_str() << std::endl;
		exit(0);
	}


	std::map<std::string, bool> gene_names;
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		int count = 0;
		int start;
		int stop;
		std::string chr;
		int overlaps = -1;
		std::string gene_name;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				chr += buffer[i];
			}
			if (count == 3 && buffer[i - 1] == '\t') {
				start = atoi(&buffer[i]);
			}

			if (count == 4 && buffer[i - 1] == '\t') {
				stop = atoi(&buffer[i]);
				//std::cout<<chr<<" "<<start<<" "<<stop<<std::endl;
				overlaps = is_overlapping(chr, start, stop, svs_elements, max_distance);
				if (overlaps == -1) {
					break;
				} else {
					std::cout << "hit" << std::endl;
				}
			}

			if (count == 8 && buffer[i] != '\t') {
				gene_name += buffer[i];
			}
			if (count > 8) {
				gene_names[gene_name] = true;
				/*		fprintf(file, "%c", '\t');
				 fprintf(file, "%s", trans(svs_elements[overlaps].type).c_str());
				 fprintf(file, "%c", '\t');
				 fprintf(file, "%s", svs_elements[overlaps].start_chr.c_str());
				 fprintf(file, "%c", '\t');
				 fprintf(file, "%i", svs_elements[overlaps].start);
				 fprintf(file, "%c", '\t');
				 fprintf(file, "%s", svs_elements[overlaps].stop_chr.c_str());
				 fprintf(file, "%c", '\t');
				 fprintf(file, "%i", svs_elements[overlaps].stop);
				 fprintf(file, "%c", '\n');
				 */
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	FILE *file;
	file = fopen(output.c_str(), "w");
	for (std::map<std::string, bool>::iterator i = gene_names.begin(); i != gene_names.end(); i++) {
		fprintf(file, "%s", (*i).first.c_str());
		fprintf(file, "%c", '\n');
	}
	myfile.close();
	fclose(file);
}
