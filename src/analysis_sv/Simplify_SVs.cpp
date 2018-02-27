/*
 * Simplify_SVs.cpp
 *
 *  Created on: Nov 28, 2017
 *      Author: sedlazec
 */

#include "Simplify_SVs.h"

std::vector<std::string> parse_support(char * vec, vector<std::string> names) {
	std::vector<std::string> accessions;
	size_t i = 0;
	while (vec[i] != ';') {
		if (vec[i] == '1') {
			accessions.push_back(names[i]);
		}
		i++;
	}
	return accessions;
}
std::string parse_seq(char * vec) {
	size_t i = 0;
	std::string chr = "";
	while (vec[i] != ';' && vec[i] != '\t') {
		chr += vec[i];
		i++;
	}
	return chr;
}
std::string parse_gene_name(char * vec) {
	size_t i = 0;
	std::string name = "";
	bool parse = false;
	while (vec[i] != ';' && vec[i] != '\t') {
		if (strncmp("Name=", &vec[i], 5) == 0) {
			i = i + 5;
			parse = true;
		}
		if (parse) {
			name += vec[i];
		}
		if (vec[i] == ',') {
			parse = false;
		}
		i++;
	}
	return name.substr(0, name.size() - 1); //chop of the last comma
}

sv_simple_str parse_line_sv(char * buffer, int buffer_size, vector<std::string> names, std::string & gene_name) {
	int count = 0;

	sv_simple_str sv;
	for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
		if (count == 0 && buffer[i] != '\t') {
			sv.start.chr += buffer[i];
		}
		if (count == 1 && buffer[i - 1] == '\t') {
			sv.start.pos = atoi(&buffer[i]);
		}
		if (count == 4 && (buffer[i] != '\t' && (buffer[i] != '<' && buffer[i] != '>'))) {
			sv.svtype += buffer[i];
		}
		if (count == 7 && strncmp(&buffer[i], "SUPP_VEC=", 9) == 0) {
			sv.accessions = parse_support(&buffer[i + 9], names);
		}
		if (count == 7 && strncmp(&buffer[i], "END=", 4) == 0) {
			sv.stop.pos = atoi(&buffer[i + 4]);
		}
		if (count == 7 && strncmp(&buffer[i], "CHR2=", 5) == 0) {
			sv.stop.chr = parse_seq(&buffer[i + 5]);
		}
		if (count == 7 && strncmp(&buffer[i], "STRANDS=", 8) == 0) {
			sv.strands.first = buffer[i + 8] == '+';
			sv.strands.second = buffer[i + 9] == '+';
		}
		if (count == 7 && strncmp(&buffer[i], ";gene_id=", 9) == 0) {
			gene_name = parse_gene_name(&buffer[i + 9]);
			break;
		}

		if (buffer[i] == '\t') {
			count++;
		}
	}
	return sv;
}

void print_gene_sv(std::string gene_name, sv_simple_str entry, FILE *&file2, int pop_size, std::map<std::string, std::vector<std::string> > populations) {

	fprintf(file2, "%s", gene_name.c_str());
	fprintf(file2, "%c", '\t');
	fprintf(file2, "%s", entry.svtype.c_str());
	fprintf(file2, "%c", '\t');
	fprintf(file2, "%s", entry.start.chr.c_str());
	fprintf(file2, "%c", ':');
	fprintf(file2, "%i", entry.start.pos);
	fprintf(file2, "%c", '-');
	fprintf(file2, "%s", entry.stop.chr.c_str());
	fprintf(file2, "%c", ':');
	fprintf(file2, "%i", entry.stop.pos);
	fprintf(file2, "%c", '\t');
	if (entry.strands.first == true) {
		fprintf(file2, "%c", '+');
	} else {
		fprintf(file2, "%c", '-');
	}

	if (entry.strands.second == true) {
		fprintf(file2, "%c", '+');
	} else {
		fprintf(file2, "%c", '-');
	}
	fprintf(file2, "%c", '\t');
	for (size_t t = 0; t < entry.accessions.size(); t++) {
		fprintf(file2, "%s", entry.accessions[t].c_str());
		if (t + 1 < entry.accessions.size()) {
			fprintf(file2, "%c", ',');
		}
	}
	fprintf(file2, "%c", '\t');
	fprintf(file2, "%i", (int) entry.accessions.size() );
	fprintf(file2, "%c", '\t');
	fprintf(file2, "%f", (double) entry.accessions.size() / (double) pop_size);
	//compute AF!
	for (std::map<std::string, std::vector<std::string> >::iterator j = populations.begin(); j != populations.end(); j++) {
		int count = 0;
		for (size_t z = 0; z < entry.accessions.size(); z++) {
			for (size_t t = 0; t < (*j).second.size(); t++) {
				if (strcmp(entry.accessions[z].c_str(), (*j).second[t].c_str()) == 0) {
					count++;
				}
			}
		}
		fprintf(file2, "%c", '\t');
		fprintf(file2, "%i", count);
		fprintf(file2, "%c", '\t');
		fprintf(file2, "%f", (double) count /(double) (*j).second.size());

	}
	fprintf(file2, "%c", '\n');
}
std::map<std::string, std::vector<std::string> > parse_populations(std::string pop_file) {
	std::map<std::string, std::vector<std::string> > pop;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	ifstream myfile;
	myfile.open(pop_file.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Pop Parser: could not open file: " << pop_file.c_str() << endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);

	while (!myfile.eof()) {
		std::string sample = "";
		std::string pop_id = "";
		int count = 0;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				sample += buffer[i];
			}
			if (count > 0 && buffer[i] != '\t') {
				if (buffer[i] == ' ') {
					pop_id += "_";
				} else {
					pop_id += buffer[i];
				}
			}
			if (buffer[i] == '\t') {
				if (!pop_id.empty()) {
					pop[pop_id].push_back(sample);
					pop_id = "";
				}
				count++;
			}
		}
		if (!pop_id.empty()) {
			pop[pop_id].push_back(sample);
			pop_id = "";
		}

		myfile.getline(buffer, buffer_size);
	}

	return pop;
}
void simplify_svs(std::string filename, std::string pop_file, int min_size, std::string output) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	ifstream myfile;

	std::map<std::string, std::vector<std::string> > population = parse_populations(pop_file);

	myfile.open(filename.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Annotation Parser: could not open file: " << filename.c_str() << endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);

	vector<std::string> names;
	std::map<std::string, std::vector<sv_simple_str> > svs;

	FILE *file2;
	file2 = fopen(output.c_str(), "w");
	fprintf(file2, "%s", "Genes\tSVtype\tPositions\tstrands\taccessions\ttotal_num\ttotal_AF");
	for (std::map<std::string, std::vector<std::string> >::iterator j = population.begin(); j != population.end(); j++) {
		fprintf(file2, "%s", "\t");
		fprintf(file2, "%s", (*j).first.c_str());
		fprintf(file2, "%s", "_num");

		fprintf(file2, "%s", "\t");
		fprintf(file2, "%s", (*j).first.c_str());
		fprintf(file2, "%s", "_AF");
	}
	fprintf(file2, "%s", "\n");

	while (!myfile.eof()) {
		if (buffer[0] == '#' && buffer[1] == 'C') {
			int count = 0;
			std::string name = "";
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count > 8 && buffer[i] != '\t') {
					name += buffer[i];
				}
				if (buffer[i] == '\t') {
					if (!name.empty()) {
						names.push_back(name);
						name.clear();
					}
					count++;
				}
			}
			if (!name.empty()) {
				names.push_back(name);
				name.clear();
			}
			cout << "Names: " << names.size() << endl;
		} else if (buffer[0] != '#') {
			std::string gene_name = "NA";
			sv_simple_str tmp = parse_line_sv(buffer, buffer_size, names, gene_name);
			print_gene_sv(gene_name, tmp, file2, names.size(), population);
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	fclose(file2);
}
