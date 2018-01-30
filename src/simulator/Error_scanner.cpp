/*
 * Error_scanner.cpp
 *
 *  Created on: Jun 30, 2017
 *      Author: sedlazec
 */

#include "Error_scanner.h"

void add_event(int pos, size_t & i, vector<differences_str> & events) {
	//insert sorted into vector:
	while (i < events.size() && pos > events[i].position) {
		i++;
	}
	differences_str ev;
	ev.position = pos;
	ev.type = 0; //mismatch
	events.insert(events.begin() + i, ev);
}

int get_index(char nuc) {
	switch (nuc) {
	case 'A':
		return 0;
		break;
	case 'C':
		return 1;
		break;
	case 'G':
		return 2;
		break;
	case 'T':
		return 3;
		break;
	case '-':
		return 4;
		break;
	default:
		return 5;
		break;
	}
	return 5;
}

int find_pos(vector<differences_str> events, int pos) {
	for (size_t i = 0; i < events.size(); i++) {
		if (events[i].position == pos) {
			return i;
		}
	}
	return -1;
}

int get_id(char ref) {
	switch (ref) {
	case 'A':
		return 0;
		break;
	case 'C':
		return 1;
		break;
	case 'G':
		return 2;
		break;
	case 'T':
		return 3;
		break;
	case '-':
		return 4;
		break;
	default:
		break;
	}

	cerr << "Unknonw base: "<<ref << endl;
	exit(1);
}

void computeAlignment(std::vector<CigarOp> cigar_data, char * md, std::string read_seq, std::vector<vector<vector<int> > > & error_mat,std::string dir) {
	int pos = 0;

	std::string ref = read_seq;
	std::string read = read_seq;

	for (size_t i = 0; i < cigar_data.size(); i++) {
		if (cigar_data[i].Type == 'I') {
			//cout << "I: " << pos << endl;
			for (int t = 0; t < cigar_data[i].Length; t++) {
				ref[pos] = '-';
				pos++;
			}
		} else if (cigar_data[i].Type == 'D') {
			for (int t = 0; t < cigar_data[i].Length; t++) {
				read.insert(pos, "-");
				ref.insert(pos, "D");
				pos++;
			}
		} else if (cigar_data[i].Type == 'S') {
			read.erase(pos, cigar_data[i].Length);
			ref.erase(pos, cigar_data[i].Length);
		} else if (cigar_data[i].Type == 'M') {
			pos += cigar_data[i].Length;
		} else if (cigar_data[i].Type == 'H') {

		}
	}

	pos = 0; //operates on the read level
	bool match = false;
	bool gap = false;
	//comp_aln = clock();
	size_t i = 0;
	size_t j = 0;
	int aln_pos = 0;
	while (md[i] != '\t') {
		if (md[i] == '^') { //deletion!
			gap = true;
		}
		if ((atoi(&md[i]) == 0 && md[i] != '0') && md[i] != '^') { //is not a number
		//i: pos on md sting
		//pos: pos on ref?
		//	cout << 67 + pos + 1 << " " << md[i] << endl;

			while (j < read.size() && aln_pos != pos) {
				if (ref[j] != '-') {
					aln_pos++;
				}
				j++;
			}
			if (!gap) { // only mismatches are stored. We should have the rest from CIGAR
				ref[j] = md[i];
				//	cout << "M: " << pos << " " << md[i] << " " << read_seq.substr(read_pos - 3, 7) << endl; //[read_pos - 1] << " " << read_seq[read_pos] << " " << read_seq[read_pos + 1] << endl;
				//store in sorted order:
				//		error_mat[get_index(md[i])][read_seq[pos]];//check that out!
				pos++; //just the pos on ref!
			} else {
				if(ref[j]=='D'){
					ref[j] = md[i];
				}else{ //TODO improve this: sometimes they differ by +1 - +2 positions ^^
					//cout<<"Error :" <<j << " "<<pos;
					int t=j-10;
					while(ref[t]!='D'){
						t++;
					}
					//cout<<" correct: "<<t<<endl;
					ref[t]= md[i];
					j=t;
				}

				//		cout << "D: " << pos << " " << md[i] << " " << ref.substr(j - 1, 3) << endl;
				pos++;
			}
			match = false;
		} else if (!match) {
			match = true;
			pos += atoi(&md[i]);
			gap = false;
		}
		i++;
	}

	cout << " ref: "<<dir<<" "<< ref << endl;
	cout << "read: "<<dir<<" "<< read << endl;
	cout << endl;
	pos = 0;

	vector<int> tmp1;
	tmp1.assign(5, 0);
	vector<vector<int> > tmp2;
	tmp2.assign(5, tmp1);
	for (size_t i = 0; i < ref.size(); i++) {
		while (pos + 1 > error_mat.size()) {
			error_mat.push_back(tmp2);
		}
		error_mat[pos][get_id(ref[i])][get_id(read[i])]++;
		if (read[i] != '-') {
			pos++;
		}
	}

}

vector<differences_str> summarizeAlignment(std::vector<CigarOp> cigar_data, char * md) { //,std::string read_seq, vector<vector<int > >error_mat) {
	//clock_t comp_aln = clock();
	vector<differences_str> events;
	int pos = 0; //this->getPosition();
	differences_str ev;

	for (size_t i = 0; i < cigar_data.size(); i++) {
		if (cigar_data[i].Type == 'D') {
			// in MD!
			ev.position = pos;
			ev.type = cigar_data[i].Length; //deletion
			events.push_back(ev);
			pos += cigar_data[i].Length;
		} else if (cigar_data[i].Type == 'I') {
			ev.position = pos;
			ev.type = cigar_data[i].Length * -1; //insertion
			//run pos - pos+len -> X vs. -
			//	error_mat['-'][get_index(read_seq[pos])]; //ref read
			events.push_back(ev);
		} else if (cigar_data[i].Type == 'M') {
			//in MD
			pos += cigar_data[i].Length;
		} else if (cigar_data[i].Type == 'N') {
			pos += cigar_data[i].Length;
		} else if (cigar_data[i].Type == 'S' && pos != 0) { /// Used for reads ranging into an inser
			string sa;
			ev.position = pos;
			ev.type = pos + cigar_data[i].Length; //such that we know when the read really ended.
			events.push_back(ev);
		}
	}
	/*
	 std::cout << "FIRST:" << std::endl;
	 for (size_t i = 0; i < events.size(); i++) {
	 if (abs(events[i].type) > 200) {
	 cout << events[i].position << " " << events[i].type << endl;
	 }
	 }
	 cout << endl;
	 */
	//set ref length requ. later on:pos
	//cout<<" comp len: "<<this->ref_len<<" "<<pos<<" "<<this->getPosition()<<endl;
	pos = 0;
	bool match = false;
	bool gap = false;
	int ref_pos = 0;
	size_t pos_events = 0;
	//comp_aln = clock();
	size_t i = 0;
	while (md[i] != '\t') {
		if (md[i] == '^') { //deletion!
			gap = true;
		}
		if ((atoi(&md[i]) == 0 && md[i] != '0')) { //is not a number
			if (!gap) { // only mismatches are stored. We should have the rest from CIGAR
				//correct for shift in position with respect to the ref:
				//
				while (ref_pos < (int) events.size() && pos > events[ref_pos].position) {
					if (events[ref_pos].type > 0) {
						pos += events[ref_pos].type;
					}
					ref_pos++;
				}
				//store in sorted order:
				add_event(pos, pos_events, events);
				//		error_mat[get_index(md[i])][read_seq[pos]];//check that out!
				pos++; //just the pos on ref!
			}
			match = false;
		} else if (!match) {
			match = true;
			pos += atoi(&md[i]);
			gap = false;
		}
		i++;
	}

	/*std::cout << "SECOND:" << std::endl;
	 for (size_t i = 0; i < events.size(); i++) {
	 if (abs(events[i].type) > 200) {
	 cout << events[i].position << " " << events[i].type << endl;
	 }
	 }
	 cout << endl;
	 */
	return events;
}

size_t get_length(std::vector<CigarOp> CigarData) {
	size_t len = 0; //orig_length;
	for (size_t i = 0; i < CigarData.size(); i++) {
		//cout<<CigarData[i].Length<<":"<<CigarData[i].Type<<";";
		if (CigarData[i].Type == 'D' || CigarData[i].Type == 'M' || CigarData[i].Type == 'N') {
			len += CigarData[i].Length;
		}
	}
	return len;
}
std::vector<CigarOp> parse_cigar(char * buffer) {
	size_t i = 0;
//	2S2M1I4M1D27M1D6M1D3M117S
	std::vector<CigarOp> cigar;
	CigarOp tmp;
	tmp.Length = -1;
	while (buffer[i] != '\t') {
		if ((atoi(&buffer[i]) == 0 && buffer[i] != '0')) { //not a number:
			tmp.Type = buffer[i];
			cigar.push_back(tmp);
			tmp.Length = -1;
		} else {
			if (tmp.Length == -1) {
				tmp.Length = atoi(&buffer[i]);
			}
		}
		i++;
	}
	return cigar;
}

void store_diffs(vector<differences_str> diffs, std::vector<read_position> & error_profile) {

	read_position tmp;
	tmp.del = 0;
	tmp.ins = 0;
	tmp.match = 0;
	tmp.mismatch = 0;
//	cout<<diffs.size() - 1<<endl;
	size_t size = diffs[diffs.size() - 1].position;
	//cout<<"size: "<<size<< " "<< error_profile.size()<<endl;
	while (size > error_profile.size()) { //check if the length is ok.
		error_profile.push_back(tmp);
	}

	int pos = 0;
	for (size_t i = 0; i + 1 < diffs.size(); i++) { //last position in that array is basically a flag!
		int diff = diffs[i].position - pos;
		//	cout << "Diff:" << diff << " " << pos << " " << error_profile.size() << endl;
		for (int j = 0; j < diff - 1; j++) {
			error_profile[pos].match++;
			error_profile[pos].total++;
			pos++;
		}
		int j = 0;
		int border = abs(diffs[i].type);
		if (border == 0) {
			border = 1;
		}
		while (j < border) {
			if (diffs[i].type < 0) { //ins
				error_profile[pos].ins++;
			} else if (diffs[i].type > 0) { //del
				error_profile[pos].del++;
			} else { //subst:
				error_profile[pos].mismatch++;
			}
			error_profile[pos].total++;
			pos++;
			j++;
		}

	}
	while (pos < (int) size) {
		error_profile[pos].match++;
		error_profile[pos].total++;
		pos++;
	}

}
void generate_error_profile(int min_length, bool comp_error_mat, std::string output) {

	std::vector<read_position> error_profile;

	std::vector<std::vector<vector<int> > > error_mat;

	//we are using cin to get the data.
	double num = 0;

	while (!cin.eof()) {
		string line;
		getline(cin, line);
		if (!cin.fail()) {
			if (line[0] != '@') {
				int count = 0;
				std::string seq;
				std::vector<CigarOp> cigar_data;
				std::string dir="";
				for (size_t i = 0; i < line.size(); i++) {
					if(count==1 && line[i]!='\t' ){
					//	if (comp_error_mat) {
							dir+=line[i];
					//	}
					}
					if (count == 5 && line[i - 1] == '\t') {
						cigar_data = parse_cigar(&line[i]); //TODO
						if ((int) get_length(cigar_data) < min_length) {
							break;
						}
					}
					if (count == 9 && line[i] != '\t' && comp_error_mat) {
						seq += line[i];
					}
					if (count > 11) {
						if (strncmp(&line[i], "MD:Z:", 5) == 0) {
							//cout << "MD!" << endl;
							num++;
							vector<differences_str> diffs = summarizeAlignment(cigar_data, &line[i + 5]);
							store_diffs(diffs, error_profile);
							if (comp_error_mat) {
								computeAlignment(cigar_data, &line[i + 5], seq, error_mat,dir);
							}
							if ((int) num % 10000 == 0) {
								cout << "Scanned: " << num << endl;
							}
						}
					}
					if (line[i] == '\t') {
						count++;
					}
				}

			}

		} else {
			break;
		}
	}

	FILE *file;
	file = fopen(output.c_str(), "w");
	fprintf(file, "%s", "Pos\tP(stop)\tP(match)\tP(mismatch)\tP(ins)]\tP(del)\n");
	for (size_t i = 0; i < error_profile.size(); i++) {
		fprintf(file, "%i", (int) i);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", 1 - (error_profile[i].total / num));
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", error_profile[i].match / error_profile[i].total);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", error_profile[i].mismatch / error_profile[i].total);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", error_profile[i].ins / error_profile[i].total);
		fprintf(file, "%c", '\t');
		fprintf(file, "%f", error_profile[i].del / error_profile[i].total);

		//cout << i << "\t" << 1 - (error_profile[i].total / num) << "\t" << error_profile[i].match << "\t" << error_profile[i].mismatch << "\t" << error_profile[i].del << "\t" << error_profile[i].ins << std::endl;
		fprintf(file, "%c", '\n');
	}
	fclose(file);

	if (comp_error_mat) {
		std::string out = output;
		out += "_errormat.txt";
		file = fopen(out.c_str(), "w");
		fprintf(file, "%s", "Pos\tRefallele\tRead(A)\tRead(C)\tRead(G)\tRead(T)\tRead(-)\n");

		for (size_t i = 0; i < error_mat.size(); i++) {

			for (size_t j = 0; j < error_mat[i].size(); j++) {
				fprintf(file, "%i", (int) i);
				fprintf(file, "%c", '\t');
				switch (j) {
				case 0:
					fprintf(file, "%c", 'A');
					break;
				case 1:
					fprintf(file, "%c", 'C');
					break;
				case 2:
					fprintf(file, "%c", 'G');
					break;

				case 3:
					fprintf(file, "%c", 'T');
					break;
				case 4:
					fprintf(file, "%c", '-');
					break;
				}

				for (size_t t = 0; t < error_mat[i][j].size(); t++) {
					fprintf(file, "%c", '\t');
					fprintf(file, "%i", error_mat[i][j][t]);
				}
				fprintf(file, "%c", '\n');
			}

		}
		fclose(file);
	}

	cout << "Number of valid reads: " << num << std::endl;

}
