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
#include "Summarize_SV.h"
#include "convert/Convert_Honey_tails.h"
#include "convert/Convert_Assemblytics.h"
#include "convert/Convert_Bionano.h"
#include "convert/Convert_VCF_to_BED.h"
#include "merge_vcf/Paramer.h"
#include "merge_vcf/combine_svs.h"
#include "convert/Process_Coverage.h"
#include "analysis_sv/GIAB_summary.h"
#include "vcfs/Detect_nested.h"
#include "simulator/Error_scanner.h"
#include "simulator/Sim_reads.h"



Parameter* Parameter::m_pInstance = NULL;
int main(int argc, char *argv[]) {
	if (argc > 1) {
		switch (atoi(argv[1])) {
		case 1:
			if (argc == 6) {
				bool coordinates = bool(atoi(argv[4]) == 0);
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
				std::cerr << "3: 0= simulated reads; 1= real reads " << std::endl;
				std::cerr << "4: output prefix" << std::endl;
			}
			break;
		case 2:
			if (argc > 2) {
				if (strcmp(argv[2], "scan") == 0) {
					if (argc == 5) {
						generate_error_profile(atoi(argv[3]), std::string(argv[4]));
					} else {
						std::cerr << "Required parameters:" << std::endl;
						std::cerr << "How to run: samtools view your_file.bam | ./SURVIVOR 2 scan 10000 error.txt"<<std::endl;
						std::cerr << "1: Min read length" << std::endl;
						std::cerr << "2: output " << std::endl;
					}
				} else if (strcmp(argv[2], "simul") == 0) {
					if (argc == 7) {
						simulate_reads(std::string(argv[3]), std::string(argv[4]),atoi(argv[5]), std::string(argv[6]));

					} else {
						std::cerr << "No parameters provided:" << std::endl;
						std::cerr << "1: Reference fasta file" << std::endl;
						std::cerr << "2: error profile file (see scan)" << std::endl;
						std::cerr << "3: Coverage" << std::endl;
						std::cerr << "4: output prefix" << std::endl;
					}


				} else {
					std::cerr << "Unkown option!" << std::endl;
				}
			} else {
				std::cerr << "Choose options:" << std::endl;
				std::cerr << "\'scan\': generate error profile" << std::endl;
				std::cerr << "\'simul\': simulate reads (Pacbio / Nanopore)" << std::endl;
				std::cerr << "We included a error profile for Nanopore and Pacbio data in the SURVIVOR folder."<<std::endl;
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
			if (argc == 6) {
				//merge VCF calls for SV
				merge_vcf(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "list of files to process" << std::endl;
				std::cerr << "Max allowed distance to be counted the same" << std::endl;
				std::cerr << "Min observed" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			break;
			/*case 5:
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
			 break;*/
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
			if (argc == 5) {
				//convert Pindel to VCF
				process_bed_file(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]));
				//process_Pindel(std::string(argv[2]), atoi(argv[3]), atof(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "Bed file" << std::endl;
				std::cerr << "Type" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			break;
		case 10:
			if (argc == 5) {
				//merge 3 SV calls from the same strain
				process_Honey(std::string(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "Honey tails file" << std::endl;
				std::cerr << "min size" << std::endl;
				std::cerr << "Output" << std::endl;
			}
			break;
		case 11:
			if (argc == 5) {
				//convert Assemblytics to VCF
				process_Assemblytics(std::string(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "Bed file from Assemblytics" << std::endl;
				std::cerr << "Min size to keep" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			break;
		case 12:
			if (argc == 5) {
				//Convert a MQ0 coverage file to bed file for filtering SV
				comp_mq0bed(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]));
			} else {
				std::cerr << "cov MQ0 file" << std::endl;
				std::cerr << "border size " << std::endl;
				std::cerr << "min coverage to be considerd " << std::endl;
			}
			break;
		case 13:
			if (argc == 6) {
				summary_SV(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), std::string(argv[5]));
				std::cout << "You can find an R script in the src/R-scripts/ to create plots given the summary output files." << std::endl;
			} else {
				std::cerr << "vcf file" << std::endl;
				std::cerr << "Min SV size (disable: -1)" << std::endl;
				std::cerr << "Max SV size (disable: -1)" << std::endl;
				std::cerr << "output summary file" << std::endl;
			}
			break;

		case 5: //prev 14!
			if (argc == 10) {
				//merge 3 SV calls from the same strain
				//	combine_calls_new(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), std::string(argv[5]));
				combine_calls_svs(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), std::string(argv[9]));
			} else {
				std::cerr << "Tab file with names" << std::endl;
				std::cerr << "max distance between breakpoints " << std::endl;
				std::cerr << "Minimum number of supporting caller" << std::endl;
				std::cerr << "Take the type into account (1==yes, else no)" << std::endl;
				std::cerr << "Take the strands of SVs into account (1==yes, else no)" << std::endl;
				std::cerr << "Estimate distance based on the size of SV (1==yes, else no)." << std::endl;
				std::cerr << "Minimum size of SVs to be taken into account." << std::endl;
				std::cerr << "Output prefix" << std::endl;
			}
			break;
		case 15:
			if (argc == 5) {
				//eval calls paper
				eval_paper(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]));
			} else {
				std::cerr << "vcf file" << std::endl;
				std::cerr << "Tab file with names" << std::endl;
				std::cerr << "max dist" << std::endl;
			}
			break;
		case 16:
			if (argc == 4) {
				//eval calls paper
				process_Bionano(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "*.smap file" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			break;
		case 17:
			if (argc == 4) {
				//eval calls paper
				process_CG(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "*.smap file" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			break;
		case 18:
			if (argc == 6) {
				parse_VCF_to_bed(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "vcf file" << std::endl;
				std::cerr << "min size" << std::endl;
				std::cerr << "max size" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			break;

		case 19:
			if (argc == 4) {
				//summarize venn
				summary_venn(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "vcf venn file" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			break;
		case 20:
			if (argc == 5) {
				filter_vcf_sniffles(std::string(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "vcf input file" << std::endl;
				std::cerr << "min length" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			break;
		case 21:
			if (argc == 3) {
				summarize_paper_gaib(std::string(argv[2]));
			} else {
				std::cerr << "vcf input file" << std::endl;
			}
			break;
		case 22:
			if (argc == 4) {
				change_insert_pos(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "vcf input file" << std::endl;
				std::cerr << "corrected vcf output file" << std::endl;
			}
			break;
		case 23:
			if (argc == 6) {
				summarize_badcoverage(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "sambada depth input file" << std::endl;
				std::cerr << "window size" << std::endl;
				std::cerr << "min coverage " << std::endl;
				std::cerr << "corrected vcf output file" << std::endl;
			}
			break;
		case 24:
			if (argc == 5) {
				summarize_VCF_files(std::string(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "input file: list of vcf" << std::endl;
				std::cerr << "min size " << std::endl;
				std::cerr << "output file" << std::endl;
			}
			break;
		case 25:
			if (argc == 4) {
				summary_giab(std::string(argv[2]), std::string(argv[3]));

			} else {
				std::cerr << "input file _venn" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			break;
		case 26:
			if (argc == 4) {
				convert_vcf(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "input vcf file" << std::endl;
				std::cerr << "output bed file" << std::endl;
			}
			break;
		case 27:
			if (argc == 4) {
				detect_nested(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "input vcf file" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			break;
		case 28:
			if (argc == 6) {
				prepare_svviz(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "input vcf file" << std::endl;
				std::cerr << "input bam file" << std::endl;
				std::cerr << "input ref file" << std::endl;
				std::cerr << "output svviz file" << std::endl;
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
		std::cerr << "9: Filter and convert SV calls from BED files" << std::endl;
		std::cerr << "10: Convert SV calls from PBHoney (tails)" << std::endl;
		std::cerr << "11: Convert SV calls from Assemblytics" << std::endl;
		std::cerr << "12: Summarize MQ 0 coverage to bed file" << std::endl;
		std::cerr << "13: Summarize SVs events in VCF file" << std::endl;
		//	std::cerr << "14: Combine calls from different vcf files" << std::endl;
	}
	return 0;
}
