/*
 * test_cov.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: fsedlaze
 */
#include "test_cov.h"

struct SV_sim {
	int pos;
	int times;
};
std::vector<SV_sim> get_SV_pos(int num_SV, int genome) {
	std::vector<SV_sim> svs;
	SV_sim tmp;
	tmp.times=0;

	svs.resize(num_SV,tmp);
	for(size_t i=0;i<svs.size();i++){
		svs[i].pos=rand() % genome;
	}
	return svs;
}

int sim_readlength(int read_length) {
	if (rand() % 100 < 50) {
		return read_length - rand() % (read_length / 10);
	}
	return read_length + rand() % (read_length / 10);
}

// boxplot(t(t[,c(2:11)]),xlab="coverage",ylab="min 5 reads supporting")

void est_cov(int read_length, int num_SV, int min_overlap, int min_support,int coverage) {
	srand(time(NULL));
	long genome=300000000;
	/* int min_overlap=50;
	 int read_length=400;
	 int num_Sv=2;
	 int min_support=5;
	 int coverage=10;
	 */
//	std::cout<<"Start"<<std::endl;
//	for (size_t cov = 1; cov < 30; cov++) { //x axis!
		//for (size_t read_length=100; read_length<20000;read_length+=100){
	//	int coverage = cov;
		cout << coverage << "\t";
		for (size_t times = 0; times < 20; times++) {
			int covered = 0;
			std::vector<SV_sim> svs = get_SV_pos(num_SV, genome); // gives X random pos.
			long num_reads = (genome / (long) read_length) * (long) coverage;
			for (long i = 0; i < num_reads; i++) {
				long pos = rand() % genome;
				int length = read_length;//sim_readlength(read_length);

				//detect overlap taking into account the min overlap + read length
				for (size_t j = 0; j < svs.size(); j++) {
					if (pos + min_overlap < svs[j].pos && (pos + length) - min_overlap > svs[j].pos) {
						svs[j].times++;
						if(svs[j].times==min_support+1){
							covered++;
						}
					}
				}
			}
			std::cout << covered << "\t";
			svs.clear();
		}
		std::cout << endl;
//	}
}
void count_valid_reads(double allowed_n_ratio) {
	int count = 0;
	while (!cin.eof()) {
		string line;
		getline(cin, line);

		if (!cin.fail()) {
			if (line[0] != '>') {
				//error
				double len = (double) line.size();
				double ns = 0;
				for (size_t i = 0; i < line.size(); i++) {
					if (line[i] == 'N' || line[i] == 'n') {
						ns++;
					}
				}
				if (ns / len < allowed_n_ratio) {
					count++;
				}
			}
		} else {
			break;
		}
	}
	cout << "Number of valid reads: " << count << std::endl;
}
