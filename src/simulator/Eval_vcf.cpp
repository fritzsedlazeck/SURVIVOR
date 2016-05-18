/*
 * Eval_vcf.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */

#include "Eval_vcf.h"
std::string trans_type(short type) {
	//0=DEL,1=DUP,2=INV,3=TRA, 4=INS
	switch (type) {
	case 0:
		return "DEL";
		break;
	case 1:
		return "DUP";
		break;
	case 2:
		return "INV";
		break;
	case 3:
		return "TRA";
		break;
	case 4:
		return "INS";
		break;
	default:
		return "NA";
		break;
	}
}
std::vector<strsimul> parse_bed_simul(std::string filename) {
	std::vector<strsimul> simulated;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "BED Parser: could not open file: " << filename.c_str()
				<< std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		int count = 0;
		strsimul tmp;
		tmp.identified = false;
		for (size_t i = 0;
				i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n';
				i++) {
			if (count == 0 && buffer[i] != '\t') {
				tmp.start.chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				tmp.start.pos = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i] != '\t') {
				tmp.stop.chr += buffer[i];
			}
			if (count == 3 && buffer[i - 1] == '\t') {
				tmp.stop.pos = atoi(&buffer[i]);
			}
			if (count == 4 && buffer[i - 1] == '\t') {
				tmp.type = get_type(std::string(&buffer[i]));
				break;
			}

			if (buffer[i] == '\t') {
				count++;
			}
		}
		simulated.push_back(tmp);
		myfile.getline(buffer, buffer_size);
	}

	myfile.close();
	return simulated;
}
bool match_coords(strsimul c1, strvcfentry c2, int max_allowed_dist) {

	if ((strcmp(c1.start.chr.c_str(), c2.start.chr.c_str()) == 0
			&& abs(c1.start.pos - c2.start.pos) < max_allowed_dist)) {
		if(c1.type==4){
			return true;
		}
		return (strcmp(c1.stop.chr.c_str(), c2.stop.chr.c_str()) == 0
				&& abs(c1.stop.pos - c2.stop.pos) < max_allowed_dist);

	} else if ((strcmp(c1.stop.chr.c_str(), c2.start.chr.c_str()) == 0
			&& abs(c1.stop.pos - c2.start.pos) < max_allowed_dist)) {
		if(c1.type==4){
				return true;
			}
		return (strcmp(c1.start.chr.c_str(), c2.stop.chr.c_str()) == 0
				&& abs(c1.start.pos - c2.stop.pos) < max_allowed_dist);

	}
	return false;

}
void reprotvals(const char * entry) {
	size_t i = 0;
	int count = 0;
	float dv = 0; //"# high-quality variant pairs"
	float dr = 0; //"# high-quality reference pairs"
	int gq = 0; //"Genotype Quality"
	int rr = 0;
	int rv = 0;
	//std::cout<<entry<<std::endl;
	while (entry[i] != '\t' && entry[i] != '\n' && entry[i] != '\0') {
		if (count == 2 && entry[i - 1] == ':') {
			gq = atoi(&entry[i]);
		}
		if (count == 5 && entry[i - 1] == ':') {
			dr = atoi(&entry[i]);
		}
		if (count == 6 && entry[i - 1] == ':') {
			dv = atoi(&entry[i]);
		}
		if (count == 7 && entry[i - 1] == ':') {
			rr = atoi(&entry[i]);
		}
		if (count == 8 && entry[i - 1] == ':') {
			rv = atoi(&entry[i]);
		}
		if (entry[i] == ':') {
			count++;
		}
		i++;
	}
	// std::cout<<dv<<" "<<rv<<" "<<dr<<" "<<rr<<" "<<gq<<std::endl;
}
strreport init_report() {
	strreport tmp;
	tmp.del = 0;
	tmp.dup = 0;
	tmp.inv = 0;
	tmp.tra = 0;
	tmp.other = 0;
	tmp.ins=0;
	return tmp;
}
std::string print_report(strreport report) {
	//0=DEL,1=DUP,2=INV,3=TRA
	std::stringstream ss; //create a stringstream
	ss << report.del; //add number to the stream
	ss << "/";
	ss << report.dup;
	ss << "/";
	ss << report.inv;
	ss << "/";
	ss << report.tra;
	ss << "/";
	//ss << report.ins;
	//ss << "/";
	ss << report.other+ report.ins;
	return ss.str();

}
void add_to_report(strreport & report, short type) {
	//0=DEL,1=DUP,2=INV,3=TRA
	switch (type) {
	case 0:
		report.del++;
		break;
	case 1:
		report.dup++;
		break;
	case 2:
		report.inv++;
		break;
	case 3:
		report.tra++;
		break;
	case 4:
		report.ins++;
		break;
	default:
		report.other++;
		break;
	}
}
void eval_calls(std::vector<strvcfentry> entries, std::vector<strsimul> simul,
		int max_allowed_dist, std::string output) {

	strreport notfound = init_report();
	strreport additional = init_report();
	FILE * right;
	FILE * addition;
	std::string out=output;
	out+="_right.vcf";
	right = fopen(out.c_str(), "w");

	out=output;
	out+="addition.vcf";
	addition=fopen(out.c_str(), "w");

	//for (size_t i = 0; i < entries.size(); i++) {
	//	std::cout<<entries[i].start.chr<<" "<<entries[i].start.pos<<" "<<entries[i].stop.chr<<" "<<entries[i].stop.pos<<" "<<entries[i].type<<std::endl;
	//}

	// std::cout<<"type dv"<<" "<<"rv"<<" "<<"dr"<<" "<<"rr"<<" "<<"gq"<<std::endl;
	for (size_t i = 0; i < entries.size(); i++) {
		bool found = false;
		for (size_t j = 0; j < simul.size(); j++) {
			if (simul[j].type == entries[i].type) {
				if (match_coords(simul[j], entries[i], max_allowed_dist)) { //check if order is perserved!
					simul[j].identified = true;
					found = true;
				}
			}
		}
		if (!found) {
			fprintf(addition, "%s", entries[i].header.c_str());
			for(std::map<std::string, std::string >::iterator tz=entries[i].calls.begin();tz!=entries[i].calls.end();tz++){
				fprintf(addition, "%s",(*tz).second.c_str());
			}
			fprintf(addition, "%c", '\n');
			add_to_report(additional, entries[i].type);
			//	std::cout<<"additional found: "<<entries[i].type<<" "<<entries[i].start.chr<<" "<<entries[i].start.pos<<" "<<entries[i].stop.chr<<" "<<entries[i].stop.pos<<std::endl;
		} else {
			fprintf(right, "%s", entries[i].header.c_str());
			for (std::map<std::string, std::string>::iterator tz =
					entries[i].calls.begin(); tz != entries[i].calls.end();
					tz++) {
				fprintf(right, "%s", (*tz).second.c_str());
			}
			fprintf(right, "%c", '\n');
		}
	}

	int refound = 0;
	for (size_t j = 0; j < simul.size(); j++) {
		if (!simul[j].identified) {
			add_to_report(notfound, simul[j].type);
		} else {
			refound++;
		}
	}
	std::cout << " Overall: " << simul.size() << " " << refound
			<< " " << print_report(notfound) << " "
			<< print_report(additional) << std::endl;
}
void summarize_simul(std::vector<strsimul> simul){
	std::vector<double> svs;
	std::vector<double> svs_length;
	for(size_t i=0;i<6;i++){
		svs.push_back(0);
		svs_length.push_back(0);
	}

	for(size_t i =0;i<simul.size();i++){
		svs[simul[i].type]++;
		if(simul[i].type!=3){ //No TRA
			svs_length[simul[i].type]+=simul[i].stop.pos-simul[i].start.pos;
		}
	}

	std::cout<<"Simulated:"<<std::endl;
	std::cout<<"type\t#events\tavg. length\n"<<std::endl;
	std::cout<<"DEL\t"<<svs[0]<<"\t"<<svs_length[0]/svs[0]<<std::endl;
	std::cout<<"DUP\t"<<svs[1]<<"\t"<<svs_length[1]/svs[1]<<std::endl;
	std::cout<<"INV\t"<<svs[2]<<"\t"<<svs_length[2]/svs[2]<<std::endl;
	std::cout<<"TRA\t"<<svs[3]<<"\t"<<svs_length[3]/svs[3]<<std::endl;
	std::cout<<"INS\t"<<svs[4]<<"\t"<<svs_length[4]/svs[4]<<std::endl;

}

void eval_vcf(std::string vcf_file, std::string bed_file, int max_allowed_dist,
		std::string output) {
	std::vector<strvcfentry> entries = parse_vcf(vcf_file);
	//prase simulated
	std::vector<strsimul> simul = parse_bed_simul(bed_file);
	summarize_simul(simul);

	//compare overlap
	eval_calls(entries, simul, max_allowed_dist, output);
}
