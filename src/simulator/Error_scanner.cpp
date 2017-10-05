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

int get_index(char nuc){
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
vector<differences_str> summarizeAlignment(std::vector<CigarOp> cigar_data, char * md){//,std::string read_seq, vector<vector<int > >error_mat) {
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
	bool gap=false;
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
				while (ref_pos < (int)events.size() && pos > events[ref_pos].position) {
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
	while (pos < (int)size) {
		error_profile[pos].match++;
		error_profile[pos].total++;
		pos++;
	}

}
void generate_error_profile(int min_length, std::string output) {

	std::vector<read_position> error_profile;

	std::vector< vector<int> > error_mat;
	vector<int> tmp ;
	tmp.assign(6,0);
	error_mat.assign(6,tmp);

	//we are using cin to get the data.
	double num = 0;
	std::string seq;
	while (!cin.eof()) {
		string line;
		getline(cin, line);
		if (!cin.fail()) {
			if (line[0] != '@') {
				int count = 0;
				std::vector<CigarOp> cigar_data;

				for (size_t i = 0; i < line.size(); i++) {
					if (count == 5 && line[i - 1] == '\t') {

						cigar_data = parse_cigar(&line[i]); //TODO
						if ((int)get_length(cigar_data) < min_length) {
							break;
						}
					}
					if(count==10 && line[i]!='\t'){
						seq+=line[i];
					}
					if (count > 11) {
						if (strncmp(&line[i], "MD:Z:", 5) == 0) {
							//cout << "MD!" << endl;
							num++;
							vector<differences_str> diffs= summarizeAlignment(cigar_data, &line[i + 5]);
							store_diffs(diffs, error_profile);
							if((int)num%10000==0){
								cout<<"Scanned: "<<num<<endl;
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

	cout << "Number of valid reads: " << num << std::endl;

}
