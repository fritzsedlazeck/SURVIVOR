/*

 * Convert_VCF_to_BED.cpp
 *
 *  Created on: Mar 3, 2015
 *      Author: fsedlaze
 */

#include "Convert_VCF_to_BED.h"

std::vector<int> get_strains(strvcfentry entry) {
	std::cout << "start" << std::endl;
	std::vector<int> names;
	for (std::map<std::string, std::string>::iterator i = entry.calls.begin(); i != entry.calls.end(); i++) {
		//it should always be just one entry;
		const char *line = (*i).second.c_str();
		int id = 0;
		std::cout << (*i).second << std::endl;
		for (size_t pos = 0; pos < (*i).second.size(); pos++) {
//			if (strncmp("1/1", &line[pos], 3)==0){//(strncmp("LowQual", &line[pos], 7) == 0|| strncmp("PASS", &line[pos], 4) == 0)) {
			//names.push_back(id);
			//	id++;
			//}else
			if (strncmp("0/0", &line[pos], 3) == 0 || strncmp("./.", &line[pos], 3) == 0) {
				id++;
			} else if (pos == 0 || (*i).second[pos - 1] == '\t') {
				names.push_back(id);
				id++;
			}
		}
	}
	std::cout << "end" << std::endl;
	return names;
}
std::string get_names(const char * buffer) {
	size_t i = 0;
	std::string name;
	// /seq/schatz/fritz/UCL/real//JB1108_dist1000_overlap.vcf_dist1000_overlap.vcf.regionfilt2_20reads.vcf

	while (buffer[i] != '\t' && buffer[i] != '\0' && buffer[i] != '\n') {

		name += buffer[i];
		i++;
	}
	//std::cout<<name<<std::endl;
	//return name;
///seq/schatz/fritz/UCL/real//JB1113_dist1000_overlap.vcf_dist1000_overlap.vcf.regionfilt2_20reads.vcf
//	int start = name.find_last_of("/") + 1;
//	int stop = name.find_first_of("_");
	//std::cout<<start<<" "<<stop<<" "<<name<<std::endl;
	return name;	//.substr(start,stop-start);
}
std::vector<std::string> parse_names(std::string vcf_file) {
	std::vector<std::string> names;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(vcf_file.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "GTF Parser: could not open file: " << vcf_file.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if (buffer[0] == '#' && buffer[1] == 'C') {

			int count = 0;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count > 9 && buffer[i - 1] == '\t') {
					//std::cout<<"Hit"<<std::endl;
					names.push_back(get_names(&buffer[i]));
					//	std::cout<<"names "<<names.size()<<std::endl;
				}

				if (buffer[i] == '\t') {
					count++;
				}
			}
			myfile.close();
			return names;
		}
		myfile.getline(buffer, buffer_size);

	}
	return names;
}

void convert_vcf(std::string vcf_file, std::string output) {
	FILE *file;
	file = fopen(output.c_str(), "w");

	std::vector<strvcfentry> entries = parse_vcf(vcf_file, 0);

	int id = 0;
	for (size_t i = 0; i < entries.size(); i++) {
		if (strcmp(entries[i].start.chr.c_str(), entries[i].stop.chr.c_str()) == 0) {
			fprintf(file, "%s", entries[i].start.chr.c_str());
			fprintf(file, "%c", '\t');
			if (entries[i].start.pos > entries[i].stop.pos) {
				int start = entries[i].start.pos;
				entries[i].start.pos = entries[i].stop.pos;
				entries[i].stop.pos = start;
			}
			fprintf(file, "%i", entries[i].start.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].stop.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%s", trans_type(entries[i].type).c_str());
			fprintf(file, "%c", '\n');
		}

	}

	fclose(file);
	exit(0);

	//fprintf(file, "%s", "chr\tStart\tStop\tType\tID\tnames\n");
	for (size_t i = 0; i < entries.size(); i++) {
		bool flag = false;
		fprintf(file, "%s", entries[i].start.chr.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", entries[i].start.pos);
		fprintf(file, "%c", '\t');
		if (strcmp(entries[i].start.chr.c_str(), entries[i].stop.chr.c_str()) == 0) {
			fprintf(file, "%i", entries[i].stop.pos);
		} else {
			flag = true;
			fprintf(file, "%i", entries[i].start.pos + 1);
		}
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", trans_type(entries[i].type).c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", id);
		fprintf(file, "%c", '\t');
		std::vector<int> ids = get_strains(entries[i]);
		fprintf(file, "%i", (int) ids.size());
		/*		fprintf(file, "%c", '\t');
		 for (size_t j = 0; j < ids.size(); j++) {
		 fprintf(file, "%s", names[ids[j]].c_str());
		 fprintf(file, "%c", '\t');
		 }
		 */
		fprintf(file, "%c", '\n');

		if (flag) { //print the ending pos
			fprintf(file, "%s", entries[i].stop.chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].stop.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].stop.pos + 1);
			fprintf(file, "%c", '\t');
			fprintf(file, "%s", trans_type(entries[i].type).c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", id);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", (int) ids.size());
			fprintf(file, "%c", '\t');
			for (size_t j = 0; j < ids.size(); j++) {
				//	fprintf(file, "%s", names[ids[j]].c_str());
				fprintf(file, "%c", '\t');
			}
			fprintf(file, "%c", '\n');
		}
		id++;
	}
}
std::vector<strvcfentry> parse_vcf_slim(std::string filename) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}
	std::vector<strvcfentry> calls;
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			strvcfentry tmp;
			tmp.sup_lumpy = 0;
			tmp.caller_id = 0;
			//std::cout<<buffer<<std::endl;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {

				if (count == 0 && buffer[i] != '\t') {
					tmp.start.chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					tmp.start.pos = atoi(&buffer[i]);
					//std::cout<<tmp.start.pos<<std::endl;
				}
				if (count == 7 && buffer[i - 1] == '\t') {
					tmp.stop = parse_stop(&buffer[i]);
					if (tmp.stop.chr.empty()) {
						tmp.stop.chr = tmp.start.chr;
					}
					tmp.strands = parse_strands(&buffer[i]);
					//std::cout<<tmp.stop.chr<<std::endl;
				}
				if (count > 9 && buffer[i - 1] == '\t') {
					if (buffer[i + 2] == '1') {
						tmp.sup_lumpy++;
					}
					tmp.caller_id++;
				}
				if (count == 4 && buffer[i - 1] == '<') {

					tmp.type = get_type(std::string(&buffer[i]));
				}
				if (count >= 0 && count < 9) {
					tmp.header += buffer[i];
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			calls.push_back(tmp);
			tmp.calls.clear();
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	std::cout << calls.size() << std::endl;
	return calls;
}
void convert_vcf_bede(std::string vcffile, int min_length, std::string output) {
	std::vector<strvcfentry> entries = parse_vcf_slim(vcffile);
	std::cout << "ENT: " << entries.size() << std::endl;
	FILE *file;
	file = fopen(output.c_str(), "w");
	int id = 1;
	fprintf(file, "%s", "Chr1\tStart\tStop\tChr2\tStart\tStop\tID\teval\tstrand1\tstrand2\ttype\tNumReadsSupport\tAlleleFreq\n");
	for (size_t i = 0; i < entries.size(); i++) {
		if (strcmp(entries[i].start.chr.c_str(), entries[i].stop.chr.c_str()) != 0 || abs(entries[i].start.pos - entries[i].stop.pos) > min_length) {
			//	bool flag = false;
			//int score=strcmp(entries[i].start.chr.c_str(),entries[i].stop.chr.c_str());
			//if(score<0 || (score==0 && entries[i].start.pos<entries[i].stop.pos)){
			fprintf(file, "%s", entries[i].start.chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].start.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].start.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%s", entries[i].stop.chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].stop.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].stop.pos);
			//}else{
			/*	fprintf(file, "%s", entries[i].stop.chr.c_str());
			 fprintf(file, "%c", '\t');
			 fprintf(file, "%i", entries[i].stop.pos);
			 fprintf(file, "%c", '\t');
			 fprintf(file, "%i", entries[i].stop.pos);
			 fprintf(file, "%c", '\t');
			 fprintf(file, "%s", entries[i].start.chr.c_str());
			 fprintf(file, "%c", '\t');
			 fprintf(file, "%i", entries[i].start.pos);
			 fprintf(file, "%c", '\t');
			 fprintf(file, "%i", entries[i].start.pos);
			 }*/
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", id);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", -1);
			fprintf(file, "%c", '\t');
			if (!entries[i].strands.first) {
				fprintf(file, "%c", '+');
			} else {
				fprintf(file, "%c", '-');
			}
			fprintf(file, "%c", '\t');
			if (!entries[i].strands.second) {
				fprintf(file, "%c", '+');
			} else {
				fprintf(file, "%c", '-');
			}
			fprintf(file, "%c", '\t');
			fprintf(file, "%s", trans_type(entries[i].type).c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%f", (double) entries[i].sup_lumpy / (double) entries[i].caller_id);

			fprintf(file, "%c", '\n');
		}
		id++;
	}
}

short get_type_bed(std::string type) {
	if (strncmp(type.c_str(), "DEL", 3) == 0) {
		return 0;
	} else if (strncmp(type.c_str(), "DUP", 3) == 0) {
		return 1;
	} else if (strncmp(type.c_str(), "INV", 3) == 0) {
		return 2;
	} else if (strncmp(type.c_str(), "INS", 3) == 0) {
		return 4;
	} else {
		std::cerr << "Unknown type! " << type << std::endl;
	}
	return -1;

}

std::string print_entry_bed(strvcfentry & region) {

//	III     5104    DEL00000002     N       <DEL>   .       LowQual IMPRECISE;CIEND=-305,305;CIPOS=-305,305;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.5.9;CHR2=III;END=15991;SVLEN=10887;CT=3to5;PE=2;MAPQ=60        GT:GL:GQ:FT:RC:DR:DV:RR:RV      1/1:-12,-0.602059,0:6:LowQual:816:0:2:0:0

	std::ostringstream convert;   // stream used for the conversion
	convert << region.start.chr;
	convert << "\t";
	convert << region.start.pos;      // insert the textual representation of 'Number' in the characters in the stream
	convert << "\t";
	convert << trans_type(region.type);
	convert << "00";
	convert << "BED\tN\t<";
	convert << trans_type(region.type);
	convert << ">\t.\tLowQual\tIMPRECISE;SVTYPE=";
	convert << trans_type(region.type);
	convert << ";SVMETHOD=BEDFILE;CHR2=";
	convert << region.stop.chr;
	convert << ";END=";
	convert << region.stop.pos;
	convert << ";SVLEN=";
	convert << region.stop.pos - region.start.pos;
	convert << ";PE=";
	convert << 1;
	convert << "\tGT:GL:GQ:FT:RC:DR:DV:RR:RV\t";
	std::stringstream s;
	s << "1/1:0,0,0:0:PASS:0:0:";
	s << 1;
	s << ":0:0";
	//std::cout<<convert.str()<<std::endl;
	return convert.str();
}

void process_bed_file(std::string bedfile, std::string type, std::string output) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(bedfile.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Bed Parser: could not open file: " << bedfile.c_str() << std::endl;
		exit(0);
	}
	FILE *file;
	file = fopen(output.c_str(), "w");

	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		int count = 0;
		strvcfentry region;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				region.start.chr += buffer[i];
				region.stop.chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				region.start.pos = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				region.stop.pos = atoi(&buffer[i]);
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		region.type = get_type_bed(type);
		fprintf(file, "%s", print_entry_bed(region).c_str());
		fprintf(file, "%c", '\n');
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	fclose(file);
}

char trans_strands(bool strand) {
	if (strand) {
		return '+';
	}
	return '-';
}

void parse_VCF_to_bed(std::string vcffile, int min_size, int max_size, std::string output) {
	std::vector<strvcfentry> entries = parse_vcf(vcffile, min_size);
	FILE *file;
	file = fopen(output.c_str(), "w");

	for (size_t i = 0; i < entries.size(); i++) {
		if (entries[i].sv_len < max_size || max_size == -1) {
			//	if ((entries[i].type == 0 && entries[i].sv_len < 3000) || (entries[i].type == 4 && entries[i].sv_len < 311)) {
			fprintf(file, "%s", entries[i].start.chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].start.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].start.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%s", entries[i].stop.chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].stop.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", entries[i].stop.pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", (int) i);
			fprintf(file, "%c", '\t');
			fprintf(file, "%c", ',');
			fprintf(file, "%c", '\t');
			fprintf(file, "%c", trans_strands(entries[i].strands.first));
			fprintf(file, "%c", '\t');
			fprintf(file, "%c", trans_strands(entries[i].strands.second));
			fprintf(file, "%c", '\t');
			fprintf(file, "%s", trans_type(entries[i].type).c_str());

			fprintf(file, "%c", '\n');

			//	fprintf(file, "%s", entries[i].stop.chr.c_str());
			//	fprintf(file, "%c", '\t');
			//	fprintf(file, "%i", entries[i].stop.pos);
			//	fprintf(file, "%c", '\t');
			//	fprintf(file, "%i", (int) i);
			//	fprintf(file, "%c", '\n');
		}
	}
	fclose(file);
}

void change_insert_pos(std::string vcffile, std::string output) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(vcffile.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Bed Parser: could not open file: " << vcffile.c_str() << std::endl;
		exit(0);
	}
	FILE *file;
	file = fopen(output.c_str(), "w");

	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {

		if (buffer[0] == '#') {
			fprintf(file, "%s", buffer);
			fprintf(file, "%c", '\n');
		} else {
			int start = 0;
			bool is_insert = false;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (buffer[i] == '\t' && start == 0) {
					start = atoi(&buffer[i + 1]);
				}
				if (strncmp("<INS>", &buffer[i], 5) == 0) {
					is_insert = true;
				}
				if (is_insert && strncmp("END=", &buffer[i], 4) == 0) {
					while (buffer[i] != ';') {
						i++;
					}
					fprintf(file, "%s", "END=");
					fprintf(file, "%i", start + 1);
				}
				fprintf(file, "%c", buffer[i]);
			}
			fprintf(file, "%c", '\n');

		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	fclose(file);
}

void prepare_positions(std::string vcffile, std::string output) {
	FILE *file;
	file = fopen(output.c_str(), "w");

	std::vector<strvcfentry> entries = parse_vcf(vcffile, 0);
	for (size_t i = 0; i < entries.size(); i++) {
		fprintf(file, "%s", entries[i].start.chr.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", entries[i].start.pos);
		fprintf(file, "%c", '\t');

		fprintf(file, "%c", ' ');
		fprintf(file, "%s", entries[i].stop.chr.c_str());
		fprintf(file, "%c", ' ');
		fprintf(file, "%i", entries[i].stop.pos);

	}
}

void prepare_svviz(std::string vcffile, std::string bam, std::string ref, std::string output) {

	FILE *file;
	file = fopen(output.c_str(), "w");

	std::vector<strvcfentry> entries = parse_vcf(vcffile, 0);
	for (size_t i = 0; i < entries.size(); i++) {
		if (entries[i].type == 3) {
			fprintf(file, "%s", "svviz -t bkend -b ");
			fprintf(file, "%s", bam.c_str());
			fprintf(file, "%c", ' ');
			fprintf(file, "%s", ref.c_str());
			fprintf(file, "%c", ' ');
			fprintf(file, "%s", entries[i].start.chr.c_str());
			fprintf(file, "%c", ' ');
			fprintf(file, "%i", entries[i].start.pos);
			fprintf(file, "%c", ' ');
			if (entries[i].strands.first) {
				fprintf(file, "%c", '+');
			} else {
				fprintf(file, "%c", '-');
			}
			fprintf(file, "%c", ' ');
			fprintf(file, "%s", entries[i].stop.chr.c_str());
			fprintf(file, "%c", ' ');
			fprintf(file, "%i", entries[i].stop.pos);

			if (entries[i].strands.second) {
				fprintf(file, "%s", " + >>res_svviz\n");
			} else {
				fprintf(file, "%s", " - >>res_svviz\n");
			}
		}
	}
	fclose(file);
}
