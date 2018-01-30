/*
 * Convert_Pindel.cpp
 *
 *  Created on: Mar 3, 2015
 *      Author: fsedlaze
 */

#include "Convert_Pindel.h"

short get_type_pind(const char * type){
	if (strncmp(type, "DEL", 3) == 0 || strncmp(type, "RPL", 3) == 0) {
		return 0;
	} else if (strncmp(type, "DUP", 3) == 0) {
		return 1;
	} else if (strncmp(type, "INV", 3) == 0) {
		return 2;
		//no interchrom!
	} else if (strncmp(type, "INS", 3) == 0) {
		return 4;
	} else {
		std::cerr << "Unknown type! "<<type << std::endl;
	}
	return -1;

}
strvcfentry create_entry(strregion region, int support, short type,int id) {
	strvcfentry tmp;
//	III     5104    DEL00000002     N       <DEL>   .       LowQual IMPRECISE;CIEND=-305,305;CIPOS=-305,305;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.5.9;CHR2=III;END=15991;SVLEN=10887;CT=3to5;PE=2;MAPQ=60        GT:GL:GQ:FT:RC:DR:DV:RR:RV      1/1:-12,-0.602059,0:6:LowQual:816:0:2:0:0
	tmp.start = region.start;
	tmp.stop= region.stop;
	tmp.type=type;
	tmp.sup_lumpy=support;

	std::ostringstream convert;   // stream used for the conversion
	convert << region.start.chr;
	convert <<  "\t";
	convert << region.start.pos;      // insert the textual representation of 'Number' in the characters in the stream
	convert <<  "\t";
	convert << trans_type(type);
	convert <<  "00";
	convert << id;
	convert << "PIN\tN\t<";
	convert << trans_type(type) ;
	if(tmp.sup_lumpy<4){
		convert << ">\t.\tLowQual\tIMPRECISE;SVTYPE=";
	}else{
		convert << ">\t.\tPASS\tIMPRECISE;SVTYPE=";
	}
	convert << trans_type(type);
	convert << ";SVMETHOD=PINDELv0.2.5a8;CHR2=";
	convert << region.stop.chr;
	convert << ";END=";
	convert << region.stop.pos;

	if(tmp.type==3){
		convert << ";SVLEN=0;PE=";
	}else{
		convert << ";SVLEN=";
		convert << region.stop.pos-region.start.pos;
		convert << ";PE=";
	}
	convert <<support;
	convert <<"\tGT:GL:GQ:FT:RC:DR:DV:RR:RV\t";
	tmp.header=convert.str();
	std::stringstream s;
	s<<"1/1:0,0,0:0:PASS:0:0:";
	s<<support;
	s<<":0:0";
	tmp.calls["lumpy"]=s.str();
	return tmp;
}

void parse_pindel(std::string pindel_vcf, std::vector<strvcfentry> & entries,
		int min_number_supporting, int min_length) {

	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;
	myfile.open(pindel_vcf.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Pindel Parser: could not open file: " << pindel_vcf.c_str()
				<< std::endl;
		exit(0);
	}
	int call_id = entries.size();
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if(buffer[0]!='#'){
			int count=0;
			strregion region;
			int support = 0;
			short type=-2;
			bool flag=false;
			for (size_t i = 0;i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n';i++) {
				if(count==0 && buffer[i]!='\t'){
					region.start.chr+=buffer[i];
					region.stop.chr+=buffer[i];
				}
				if(count==1 && buffer[i-1]=='\t'){
					region.start.pos=atoi(&buffer[i]);
				}
				if(count==7 && strncmp(&buffer[i],"END=",4)==0){
					region.stop.pos=atoi(&buffer[i+4]);
				}
				if(count==7 && strncmp(&buffer[i],"SVLEN=0",7)==0){
					flag=true;
				}
				if(count==7 && strncmp(&buffer[i],"SVTYPE=",7)==0){
					type=get_type_pind(&buffer[i+7]);
				}
				if(count==9 && (buffer[i-1]==',' || (flag && buffer[i-1]==':'))){
					support=atoi(&buffer[i]);

				}
				if(buffer[i]=='\t'){
					count++;
				}
			}
			//std::cout<<support<<" "<<region.stop.pos-region.start.pos<<std::endl;
			if(support>min_number_supporting && (region.stop.pos-region.start.pos > min_length || flag) ){
				entries.push_back(create_entry(region,support,type,call_id));
				call_id++;
			}

		}
		myfile.getline(buffer, buffer_size);
	}
}


void process_Pindel(std::string pindel_vcf, int min_number_supporting,
		int min_length, std::string output) {
	std::vector<strvcfentry> entries; //= parse_vcf(delly_vcf); //get delly calls

	parse_pindel(pindel_vcf, entries, min_number_supporting, min_length);

	print_header(pindel_vcf, output);
	print_entries(output, entries);
}
