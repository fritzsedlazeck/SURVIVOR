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
void summary_SV(std::string filename, std::string output) {
	std::vector<int> len_Del;
	std::vector<int> len_Dup;
	std::vector<int> len_Inv;
	std::vector<int> len_Ins;
	int TRA = 0;
	int step = 1000;
	std::map<std::string, std::map<std::string, int> > SV_chrs;

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			std::string chr;

			int dist = 0;
			std::string type;
			//cout<<buffer<<endl;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 0 && buffer[i] != '\t') {
					chr += buffer[i];
				}
				if (count == 4 && (buffer[i] != '\t' && (buffer[i] != '>' && buffer[i] != '<'))) {
					type += buffer[i];
				}
				if (count == 7 && strncmp("SVLEN=", &buffer[i], 6) == 0) {
					dist = atoi(&buffer[i + 6]);
					break;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}


			dist = dist / step;
			//cout<<dist<<endl;
			if (SV_chrs.find(chr) != SV_chrs.end() || SV_chrs[chr].find(type) != SV_chrs[chr].end()) {
				SV_chrs[chr][type]++;
			} else {
				SV_chrs[chr][type] = 1;
			}

			if (strcmp(type.c_str(), "DEL") == 0) {
				adjust(len_Del, dist);
				len_Del[dist]++;
			} else if (strcmp(type.c_str(), "INV") == 0) {
				adjust(len_Inv, dist);
				len_Inv[dist]++;
			} else if (strcmp(type.c_str(), "DUP") == 0) {
				adjust(len_Dup, dist);
				len_Dup[dist]++;
			} else if (strcmp(type.c_str(), "INS") == 0) {
				adjust(len_Ins, dist);
				len_Ins[dist]++;
			} else if (strcmp(type.c_str(), "TRA") == 0) {
				TRA++;
			}
		}
		myfile.getline(buffer, buffer_size);
	}

	std::cout << "Parsing done: " << TRA << endl;
	myfile.close();

	int maxim = std::max(std::max((int) len_Del.size(), (int) len_Dup.size()), (int) len_Inv.size());
	FILE * file;
	file = fopen(output.c_str(), "w");
	fprintf(file, "%s", "Len(max)\tDel(1kb)\tDup(1kb)\tInv(1kb)\tINS(1kb)\tTRA(1kb)\n");
	for (size_t i = 0; i < maxim + 1; i++) {
		fprintf(file, "%i",(int)(i+1)*step);
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

	std::string out = output;
	out += "_CHR";
	file = fopen(out.c_str(), "w");
	bool flag = true;
	for (std::map<std::string, std::map<std::string, int> >::iterator i = SV_chrs.begin(); i != SV_chrs.end(); i++) {

		if (flag) { //print the header:
			fprintf(file, "%s", "Chr\tDEL\tDUP\tINV\tINS\tTRA\n");
			flag = false;
		}

		fprintf(file, "%s", (*i).first.c_str());
		fprintf(file, "%c", '\t');
		if ((*i).second.find("DEL") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["DEL"]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if ((*i).second.find("DUP") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["DUP"]);
		} else {
			fprintf(file, "%i", 0);
		}

		fprintf(file, "%c", '\t');
		if ((*i).second.find("INV") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["INV"]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if ((*i).second.find("INS") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["INS"]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\t');
		if ((*i).second.find("TRA") != (*i).second.end()) {
			fprintf(file, "%i", (*i).second["TRA"]);
		} else {
			fprintf(file, "%i", 0);
		}
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}
