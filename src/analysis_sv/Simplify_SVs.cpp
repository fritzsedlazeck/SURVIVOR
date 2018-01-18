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
sv_simple_str parse_line_sv(char * buffer, int buffer_size, vector<std::string> names, std::string & gene_name) {
	int count = 0;
	gene_name = "";
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
			gene_name = parse_seq(&buffer[i + 9]);
			break;
		}

		if (buffer[i] == '\t') {
			count++;
		}
	}
	return sv;
}

void simplify_svs(std::string filename, int min_size, std::string output) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	ifstream myfile;

	myfile.open(filename.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Annotation Parser: could not open file: " << filename.c_str() << endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);

	vector<std::string> names;
	std::map<std::string, std::vector<sv_simple_str> > svs;

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
			std::string gene_name;
			sv_simple_str tmp = parse_line_sv(buffer, buffer_size, names, gene_name);
			svs[gene_name].push_back(tmp);
		}
		myfile.getline(buffer, buffer_size);
	}

	myfile.close();

	FILE *file2;
	file2 = fopen(output.c_str(), "w");
	fprintf(file2, "%s", "Genes\tSVtype\tPositions\tstrands\taccessions\n");
	for (std::map<std::string, std::vector<sv_simple_str> >::iterator i = svs.begin(); i != svs.end(); i++) {
		for (size_t j = 0; j < (*i).second.size(); j++) {
			fprintf(file2, "%s", (*i).first.c_str());
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%s", (*i).second[j].svtype.c_str());
			fprintf(file2, "%c", '\t');
			fprintf(file2, "%s", (*i).second[j].start.chr.c_str());
			fprintf(file2, "%c", ':');
			fprintf(file2, "%i", (*i).second[j].start.pos);
			fprintf(file2, "%c", '-');
			fprintf(file2, "%s", (*i).second[j].stop.chr.c_str());
			fprintf(file2, "%c", ':');
			fprintf(file2, "%i", (*i).second[j].stop.pos);
			fprintf(file2, "%c", '\t');
			if ((*i).second[j].strands.first == true) {
				fprintf(file2, "%c", '+');
			} else {
				fprintf(file2, "%c", '-');
			}

			if ((*i).second[j].strands.second == true) {
				fprintf(file2, "%c", '+');
			} else {
				fprintf(file2, "%c", '-');
			}

			for (size_t t = 0; t < (*i).second[j].accessions.size(); t++) {
				fprintf(file2, "%c", '\t');
				fprintf(file2, "%s", (*i).second[j].accessions[t].c_str());
			}
			fprintf(file2, "%c", '\n');
		}
	}
	fclose(file2);

}
