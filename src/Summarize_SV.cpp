/*
 * Summarize_SV.cpp
 *
 *  Created on: Nov 18, 2015
 *      Author: fsedlaze
 */

#include "Summarize_SV.h"

void adjust(std::vector<int> & vec, int dist) {
	while (vec.size() < dist + 1) {
		vec.push_back(0);
	}
}

/* Rcode:
 if (!require("RColorBrewer")) {
 #install.packages("RColorBrewer")
 library(RColorBrewer)
 }
 cols=(brewer.pal(5,"Set1"))
 pdf('sniffels01_s5_len.pdf')
 t=read.table('sniffels01_s5_05_05_calls.summary',header=T)
 plot(log10(t[,1]),t[,2],xlab="log10(length(bp))",ylab="# of SVs",type='l',col=cols[1])
 lines(log10(t[,1]),t[,3],col=cols[2])
 lines(log10(t[,1]),t[,4],col=cols[3])
 lines(log10(t[,1]),t[,5],col=cols[4])
 legend('topright',legend=c('DEL','DUP','INV','INS'),lwd=2,col=cols)
 dev.off()

 pdf('sniffels01_s5_chr.pdf')
 t=read.table('sniffels01_s5_05_05_calls.summary_CHR',header=T)
 plot(c(1:length(t[,1])),t[,2],ylim=c(0,max(t[,c(2:5)])),ylab="# of SVs",col=cols[1],xlab="chromosome")
 points(c(1:length(t[,1]))-0.01,t[,3],col=cols[2])
 points(c(1:length(t[,1]))+0.01,t[,4],col=cols[3])
 points(c(1:length(t[,1]))+0.05,t[,5],col=cols[4])
 legend('topright',legend=c('DEL','DUP','INV','INS'),lwd=2,col=cols)
 dev.off()
 *
 *
 */
int get_support(vector<int> support) {
	int count = 0;
	for (size_t i = 0; i < support.size(); i++) {
		if (support[i] > 0) {
			count++;
		}
	}
	return count;
}
void summary_SV(std::string vcf_file, std::string output) {
	std::vector<int> len_Del;
	std::vector<int> len_Dup;
	std::vector<int> len_Inv;
	std::vector<int> len_Ins;
	std::vector<int> support;
	int TRA = 0;
	int step = 1000;
	std::map<std::string, std::map<int, int> > SV_chrs;

	std::vector<strvcfentry> entries = parse_vcf(vcf_file, 0);

	for (size_t i = 0; i < entries.size(); i++) {
		//summarize the support:
		int id = get_support(entries[i].caller_supports);
		if (id != 0) {
			while (id >= support.size()) {
				support.push_back(0);
			}
			if (id == 1) {
				std::cout << entries[i].start.pos << std::endl;
			}
			support[id]++;
		}
		int dist = abs(entries[i].stop.pos - entries[i].start.pos) / step;
		if (entries[i].type == 0) {
			adjust(len_Del, dist);
			len_Del[dist]++;
		} else if (entries[i].type == 1) {
			adjust(len_Inv, dist);
			len_Inv[dist]++;
		} else if (entries[i].type == 2) {
			adjust(len_Dup, dist);
			len_Dup[dist]++;
		} else if (entries[i].type == 4) {
			adjust(len_Ins, dist);
			len_Ins[dist]++;
		} else if (entries[i].type == 3) {
			TRA++;
		}
		std::string chr = entries[i].start.chr;
		if (SV_chrs.find(chr) != SV_chrs.end() || SV_chrs[chr].find(entries[i].type) != SV_chrs[chr].end()) {
			SV_chrs[chr][entries[i].type]++;
		} else {
			SV_chrs[chr][entries[i].type] = 1;
		}
	}
	std::cout << "Parsing done: " << TRA << endl;
	FILE * file;
	std::string out = output;
	out += "support";
	file = fopen(out.c_str(), "w");
	for (size_t i = 0; i < support.size(); i++) {
		fprintf(file, "%i", (int) (i));
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", support[i]);
		fprintf(file, "%c", '\n');
	}
	fclose(file);
	int maxim = std::max(std::max((int) len_Del.size(), (int) len_Dup.size()), (int) len_Inv.size());
	file = fopen(output.c_str(), "w");
	fprintf(file, "%s", "Len(max)\tDel(1kb)\tDup(1kb)\tInv(1kb)\tINS(1kb)\tTRA(1kb)\n");
	for (size_t i = 0; i < maxim + 1; i++) {
		fprintf(file, "%i", (int) (i + 1) * step);
		fprintf(file, "%c", '\t');
		if (i < len_Del.size()) {
			fprintf(file, "%i", len_Del[i]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if (i < len_Dup.size()) {
			fprintf(file, "%i", len_Dup[i]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if (i < len_Inv.size()) {
			fprintf(file, "%i", len_Inv[i]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if (i < len_Ins.size()) {
			fprintf(file, "%i", len_Ins[i]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if (i == 0) {
			fprintf(file, "%i", TRA);
		} else {
			fprintf(file, "%i", 0);
		}

		fprintf(file, "%c", '\n');
	}
	fclose(file);

	out = output;
	out += "_CHR";
	file = fopen(out.c_str(), "w");
	bool flag = true;
	for (std::map<std::string, std::map<int, int> >::iterator i = SV_chrs.begin(); i != SV_chrs.end(); i++) {

		if (flag) { //print the header:
			fprintf(file, "%s", "Chr\tDEL\tDUP\tINV\tINS\tTRA\n");
			flag = false;
		}

		fprintf(file, "%s", (*i).first.c_str());
		fprintf(file, "%c", '\t');
		if ((*i).second.find(0) != (*i).second.end()) {
			fprintf(file, "%i", (*i).second[0]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if ((*i).second.find(1) != (*i).second.end()) {
			fprintf(file, "%i", (*i).second[1]);
		} else {
			fprintf(file, "%i", 0);
		}

		fprintf(file, "%c", '\t');
		if ((*i).second.find(2) != (*i).second.end()) {
			fprintf(file, "%i", (*i).second[2]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if ((*i).second.find(4) != (*i).second.end()) {
			fprintf(file, "%i", (*i).second[4]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if ((*i).second.find(3) != (*i).second.end()) {
			fprintf(file, "%i", (*i).second[3]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}

void summary_venn(std::string filename, std::string output) {
	std::vector<std::vector<int> > mat;
	std::vector<std::string> names;

	size_t buffer_size = 200000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);

	int num = 0;
	while (!myfile.eof()) {
		if (num == 0) {
			std::string name;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (buffer[i] == '\t') {
					names.push_back(name);
					name.clear();
				} else {
					name += buffer[i];
				}
			}
			if (!name.empty()) {
				names.push_back(name);
			}
			std::vector<int> tmp;
			tmp.assign(names.size(), 0);
			mat.assign(names.size(), tmp);
			num++;
		} else {
			int count = 0;
			std::vector<int> ids;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (buffer[i - 1] == '\t' && buffer[i] == '1') {
					ids.push_back(count - 1); //-1 since we have the ID in the first column.
					//	std::cout<<count<<std::endl;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}

			for (size_t i = 0; i < ids.size(); i++) {
				for (size_t j = 0; j < ids.size(); j++) {
					if (ids[i] <= ids[j]) {
						//if (i != j || ids.size() == 1) {
							mat[ids[i]][ids[j]]++;
						//}
					}
				}
			}
			num++;
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	FILE * file = fopen(output.c_str(), "w");
	//fprintf(file, "%c", '\t');
	for (size_t i = 0; i < names.size(); i++) {
		//fprintf(file, "%i", (int) i);
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", names[i].c_str());
		//fprintf(file, "%c", '\n');
	}
	fprintf(file, "%c", '\t');
	/*for (size_t j = 0; j < names.size(); j++) {
		fprintf(file, "%i", (int) j);
		fprintf(file, "%c", '\t');
	}
	fprintf(file, "%c", '\n');
*/
	for (size_t i = 0; i < mat.size(); i++) {
		fprintf(file, "%i", (int) i);
		fprintf(file, "%c", '\t');
		for (size_t j = 0; j < mat.size(); j++) {
			fprintf(file, "%i", mat[i][j]);
			fprintf(file, "%c", '\t');
		}
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}
