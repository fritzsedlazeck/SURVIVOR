//============================================================================
// Name        : Arnes_little_helper.cpp
// Author      : Fritz Sedlazeck
// Version     :
// Copyright   : artistic license
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simulator/SV_Simulator.h"
#include "simulator/Pac_Simulator.h"
#include "simulator/Eval_vcf.h"
#include "vcfs/Combine_3_VCF.h"
#include "vcfs/Annotate_vcf.h"
#include "vcfs/Filter_vcf.h"
#include "convert/Process_Lumpy.h"
#include "convert/Convert_Pindel.h"
#include "convert/ConvertMQ0Bed.h"
int main(int argc, char *argv[]) {
	if (argc > 1) {
		switch (atoi(argv[1])) {
		case 1:
			if (argc == 6) {
				bool coordinates = bool(atoi(argv[4])==0);
				simulate_SV(std::string(argv[2]), std::string(argv[3]), coordinates, std::string(argv[5]));
				std::cout << "SV simulated" << std::endl;
			} else if (argc == 3) {
				generate_parameter_file(std::string(argv[2]));
				std::cout << "Parameter file generated" << std::endl;
			} else {
				std::cerr << "No parameters provided:" << std::endl;
				std::cerr << "To generate a example parameter file call the option again and specify the file name:" << std::endl;
				std::cerr << "To simulate SV:" << std::endl;
				std::cerr << "1: Reference fasta file" << std::endl;
				std::cerr << "2: Parameter file" << std::endl;
				std::cerr << "3: 0/1 indicating if coordinates should be with respect to reference or new reads" << std::endl;
				std::cerr << "4: output prefix" << std::endl;
			}
			break;
		case 2:
			if (argc == 4) {
				simulate_pac(std::string(argv[2]), std::string(argv[3]));
				std::cout << "SV simulated" << std::endl;
			} else {
				std::cerr << "No parameters provided:" << std::endl;
				std::cerr << "1: Reference fasta file" << std::endl;
				std::cerr << "2: output prefix" << std::endl;
			}
			break;
		case 3:
			if (argc == 6) {
				//eval VCF calls for SV
				eval_vcf(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "VCF file" << std::endl;
				std::cerr << "BED file from simulations" << std::endl;
				std::cerr << "Max allowed distance to simulated" << std::endl;
				std::cerr << "Output eval file" << std::endl;
			}
			break;
		case 4:
			if (argc == 5) {
				//merge VCF calls for SV
				merge_vcf(std::string(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "list of files to process" << std::endl;
				std::cerr << "Max allowed distance to be counted the same" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			break;
		case 5:
			if (argc == 7) {
				//merge 3 SV calls from the same strain
				combine_calls(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), atoi(argv[5]), std::string(argv[6]));
			} else {
				std::cerr << "VCF delly" << std::endl;
				std::cerr << "VCF lumpy" << std::endl;
				std::cerr << "VCF pindel" << std::endl;
				std::cerr << "max dist" << std::endl;
				std::cerr << "Output prefix" << std::endl;
			}
			break;
		case 6:
			if (argc == 6) {
				generate_gene_list(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "VCF file" << std::endl;
				std::cerr << "filtered gene file" << std::endl;
				std::cerr << "overlap" << std::endl;
				std::cerr << "outputfile" << std::endl;
			}
			break;
		case 7:
			if (argc == 8) {
				//filter SV calls delly
				filter_vcf(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), atof(argv[5]), atoi(argv[6]), std::string(argv[7]));
			} else {
				std::cerr << "VCF file to filter" << std::endl;
				std::cerr << "BED file with regions to ignore" << std::endl;
				std::cerr << "Min number of supporting pairs" << std::endl;
				std::cerr << "min ratio of support vs. ref pairs " << std::endl;
				std::cerr << "Max Genotype Quality" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			break;
		case 8:
			if (argc == 6) {
				//convert Lumpy to VCF
				process_Lumpy(std::string(argv[2]), atoi(argv[3]), atof(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "Bede file from Lumpy" << std::endl;
				std::cerr << "Min number of supporting reads" << std::endl;
				std::cerr << "Max eval" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			break;
		case 9:
			if (argc == 6) {
				//convert Pindel to VCF
				process_Pindel(std::string(argv[2]), atoi(argv[3]), atof(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "Bede file from Lumpy" << std::endl;
				std::cerr << "Min number of supporting reads" << std::endl;
				std::cerr << "Min length" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			break;

		case 10:
			if (argc == 5) {
				//Convert a MQ0 coverage file to bed file for filtering SV
				comp_mq0bed(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]));
			} else {
				std::cerr << "cov MQ0 file" << std::endl;
				std::cerr << "border size " << std::endl;
				std::cerr << "min coverage to be considerd " << std::endl;
			}
			break;
		default:
			break;
		}

	} else {
		std::cerr << "Possible options" << std::endl;
		std::cerr << "1: Simulate SV on genome" << std::endl;
		std::cerr << "2: Simulate PacBio reads" << std::endl;
		std::cerr << "3: Evaluate SV calling" << std::endl;
		std::cerr << "4: Merge SV calls (vcf) " << std::endl;
		std::cerr << "5: Consensus call from 2/3 callers" << std::endl;
		std::cerr << "6: Extract genes influenced by SVs" << std::endl;
		std::cerr << "7: Filter and convert SV calls from Delly" << std::endl;
		std::cerr << "8: Filter and convert SV calls from Lumpy" << std::endl;
		std::cerr << "9: Filter and convert SV calls from Pindel" << std::endl;
		std::cerr << "10: Summarize MQ 0 coverage to bed file" << std::endl;
	}
	return 0;
}
