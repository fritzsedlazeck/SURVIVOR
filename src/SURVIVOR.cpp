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
#include "analysis_sv/Summ_mat.h"
#include "vcfs/Generate_distMat.h"
#include "snp_overlap/Overlap_snps.h"
#include "convert/Convert_MUMmer.h"
#include "simulator/test_cov.h"
#include "analysis_sv/Simplify_SVs.h"
#include "analysis_sv/MUMmer_overlap.h"
#include "analysis_sv/Select_samples.h"
#include "convert/Convert_hapcut2.h"
#include "convert/Update_bam_pacbio.h"

Parameter* Parameter::m_pInstance = NULL;

void official_interface(int argc, char *argv[]) {
	if (argc > 1) {
		if (strcmp(argv[1], "simSV") == 0) {
			if (argc == 7) {
				bool coordinates = bool(atoi(argv[5]) == 0);
				simulate_SV(std::string(argv[2]), std::string(argv[3]), atof(argv[4]), coordinates, std::string(argv[6]));
				std::cout << "Done: SV+SNP simulated" << std::endl;
			} else if (argc == 3) {
				generate_parameter_file(std::string(argv[2]));
				std::cout << "Parameter file generated" << std::endl;
			} else {

				std::cerr << "No parameters provided:" << std::endl;
				std::cerr << "To generate a example parameter file call the option again and specify the file name:" << std::endl;
				std::cerr << std::endl;
				std::cerr << "To simulate SV:" << std::endl;
				std::cerr << "1: Reference fasta file" << std::endl;
				std::cerr << "2: Parameter file" << std::endl;
				std::cerr << "3: SNP mutations frequency (0-1)" << std::endl;
				std::cerr << "4: 0= simulated reads; 1= real reads " << std::endl;
				std::cerr << "5: output prefix" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "scanreads") == 0) {

			if (argc == 4) {
				//	bool error_mat = bool(atoi(argv[3]) == 1);
				generate_error_profile(atoi(argv[2]), false, std::string(argv[3]));
			} else {
				std::cerr << "Required parameters:" << std::endl;
				std::cerr << "How to run: samtools view your_file.bam | ./SURVIVOR scanreads 1000 error.txt" << std::endl;
				std::cerr << "1: Min read length" << std::endl;
				//	std::cerr << "2: Comp error mat (1, else not)" << std::endl;
				std::cerr << "2: output " << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "simreads") == 0) {

			if (argc == 6) {
				simulate_reads(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "No parameters provided:" << std::endl;
				std::cerr << "1: Reference fasta file" << std::endl;
				std::cerr << "2: error profile file (see scan)" << std::endl;
				std::cerr << "3: Coverage" << std::endl;
				std::cerr << "4: output prefix" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "eval") == 0) {
			if (argc == 6) {
				//eval VCF calls for SV
				eval_vcf(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "VCF file" << std::endl;
				std::cerr << "BED file from simulations" << std::endl;
				std::cerr << "Max allowed distance to simulated" << std::endl;
				std::cerr << "Output eval file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "smaptovcf") == 0) {
			if (argc == 4) {
				process_Bionano(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "*.smap file" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "merge") == 0) {
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
			exit(0);
		} else if (strcmp(argv[1], "filter") == 0) {
			if (argc == 9) {
				//filter SV calls delly
				filter_vcf(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]), atoi(argv[7]), std::string(argv[8]));
			} else {
				std::cerr << "VCF file to filter" << std::endl;
				std::cerr << "BED file with regions to ignore (NA to disable)" << std::endl;
				std::cerr << "Min SV size (-1 to disable)" << std::endl;
				std::cerr << "Max SV size (-1 to disable)" << std::endl;
				std::cerr << "Min allele frequency (0-1)" << std::endl;
				std::cerr << "Min number of reads support: RE flag (-1 to disable)" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "vcftobed") == 0) {
			if (argc == 6) {
				parse_VCF_to_bed(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "vcf file" << std::endl;
				std::cerr << "min size" << std::endl;
				std::cerr << "max size" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "bedtovcf") == 0) {
			if (argc == 5) {
				//convert Pindel to VCF
				process_bed_file(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]));
				//process_Pindel(std::string(argv[2]), atoi(argv[3]), atof(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "Bed file" << std::endl;
				std::cerr << "Type" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "bincov") == 0) {
			if (argc == 5) {
				//Convert a MQ0 coverage file to bed file for filtering SV
				comp_mq0bed(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]));
			} else {
				std::cerr << "cov MQ0 file" << std::endl;
				std::cerr << "border size " << std::endl;
				std::cerr << "min coverage to be considerd " << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "stats") == 0) {
			if (argc == 7) {
				summary_SV(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), std::string(argv[6]));
				std::cout << "You can find an R script in the src/R-scripts/ to create plots given the summary output files." << std::endl;
			} else if (argc == 5) {
				summary_SV_stream(atoi(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "vcf file" << std::endl;
				std::cerr << "Min SV size (disable: -1)" << std::endl;
				std::cerr << "Max SV size (disable: -1)" << std::endl;
				std::cerr << "Min number read support (disable: -1)" << std::endl;
				std::cerr << "output summary file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "compMUMMer") == 0) {
			if (argc == 6) {
				overlapp_mummer(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "SVs VCF file" << std::endl;
				std::cerr << "Output of MUMMer Show-diff" << std::endl;
				std::cerr << "Allowed distance (e.g. 100bp)" << std::endl;
				std::cerr << "Output: annotated vcf file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "select_sample") == 0) {
			if (argc == 4) {
				select_greedy(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "SVs VCF file" << std::endl;
				std::cerr << "Output: ranked file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "bedpetovcf") == 0) {
			if (argc == 4) {
				process_Lumpy(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "Bedpe file" << std::endl;
				std::cerr << "Output: vcf file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "hapcuttovcf") == 0) {
			if (argc == 5) {
				process_hapcut(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "original SNP file" << std::endl;
				std::cerr << "Hapcut2 final file" << std::endl;
				std::cerr << "Output: vcf file" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "convertAssemblytics") == 0) {
			if (argc == 5) {
				process_Assemblytics(std::string(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "Bed file from Assemblytics" << std::endl;
				std::cerr << "Min size to keep" << std::endl;
				std::cerr << "Output vcf file" << std::endl;
			}
			exit(0);
		}
		/*else if (strcmp(argv[1], "updateBamfile") == 0) {
		 if (argc == 5) {
		 process_sam_forpacbio(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]));
		 } else {
		 std::cerr << "original SNP file" << std::endl;
		 std::cerr << "Hapcut2 final file" << std::endl;
		 std::cerr << "Output: vcf file" << std::endl;
		 }
		 exit(0);

		 }*/

	}
	std::cerr << "Program: SURVIVOR (Tools for Structural Variations in the VCF format)" << std::endl;
	std::cerr << "Version: " << Parameter::Instance()->version << std::endl;
	std::cerr << std::endl;
	std::cerr << "Usage: SURVIVOR <command> [options]" << std::endl;
	std::cerr << std::endl;

	std::cerr << "Commands:" << std::endl;
	std::cerr << "-- Simulation/ Evaluation" << std::endl;
	std::cerr << "\tsimSV\tSimulates SVs and SNPs on a reference genome." << std::endl;
	std::cerr << "\tscanreads\tObtain error profiles form mapped reads for simulation." << std::endl;
	std::cerr << "\tsimreads\tSimulates long reads (Pacio or ONT)." << std::endl;
	std::cerr << "\teval\tEvaluates a VCF file after SV calling over simulated data." << std::endl;
	std::cerr << std::endl;

	std::cerr << "-- Comparison/filtering" << std::endl;
	std::cerr << "\tmerge\tCompare or merge VCF files to generate a consensus or multi sample vcf files." << std::endl;
	std::cerr << "\tfilter\tFilter a vcf file based on size and/or regions to ignore" << std::endl;
	std::cerr << "\tstats\tReport multipe stats over a VCF file" << std::endl;
	std::cerr << "\tcompMUMMer\tAnnotates a VCF file with the breakpoints found with MUMMer (Show-diff)." << std::endl;
	std::cerr << std::endl;

	std::cerr << "-- Conversion" << std::endl;
	std::cerr << "\tbincov\tBins coverage vector to a bed file to filter SVs in low MQ regions" << std::endl;
	std::cerr << "\tvcftobed\tConverts a VCF file to a bed file" << std::endl;
	std::cerr << "\tbedtovcf\tConverts a bed file to a VCF file " << std::endl;
	std::cerr << "\tsmaptovcf\tConverts the smap file to a VCF file (beta version)" << std::endl;
	std::cerr << "\tbedpetovcf\tConverts a bedpe file ot a VCF file (beta version)" << std::endl;
	std::cerr << "\thapcuttovcf\tConverts the Hapcut2 final file to a VCF file using the original SNP file provided to Hapcut2" << std::endl;
	std::cerr << "\tconvertAssemblytics\tConverts Assemblytics to a VCF file" <<std::endl;

	exit(0);
}

//todo output: pb 0 -> 1

int main(int argc, char *argv[]) {

	official_interface(argc, argv);
	if (argc > 1) {
		switch (atoi(argv[1])) {
		case 1:
			if (argc == 7) {
				bool coordinates = bool(atoi(argv[5]) == 0);
				simulate_SV(std::string(argv[2]), std::string(argv[3]), atof(argv[4]), coordinates, std::string(argv[6]));
				std::cout << "Done: SV+SNP simulated" << std::endl;
			} else if (argc == 3) {
				generate_parameter_file(std::string(argv[2]));
				std::cout << "Parameter file generated" << std::endl;
			} else {
				std::cerr << "No parameters provided:" << std::endl;
				std::cerr << "To generate a example parameter file call the option again and specify the file name:" << std::endl;
				std::cerr << "To simulate SV:" << std::endl;
				std::cerr << "1: Reference fasta file" << std::endl;
				std::cerr << "2: Parameter file" << std::endl;
				std::cerr << "3: SNP mutations frequency (0-1)" << std::endl;
				std::cerr << "4: 0= simulated reads; 1= real reads " << std::endl;
				std::cerr << "5: output prefix" << std::endl;
			}
			break;
		case 2:
			if (argc > 2) {
				if (strcmp(argv[2], "scan") == 0) {
					if (argc == 6) {
						bool error_mat = bool(atoi(argv[4]) == 1);
						generate_error_profile(atoi(argv[3]), error_mat, std::string(argv[5]));
					} else {
						std::cerr << "Required parameters:" << std::endl;
						std::cerr << "How to run: samtools view your_file.bam | ./SURVIVOR 2 scan 10000 error.txt" << std::endl;
						std::cerr << "1: Min read length" << std::endl;
						std::cerr << "2: Comp error mat (1, else not)" << std::endl;
						std::cerr << "3: output " << std::endl;
					}
				} else if (strcmp(argv[2], "simul") == 0) {
					if (argc == 7) {
						simulate_reads(std::string(argv[3]), std::string(argv[4]), atoi(argv[5]), std::string(argv[6]));
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
				std::cerr << "We included a error profile for Nanopore and Pacbio data in the SURVIVOR folder." << std::endl;
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
			if (argc == 4) {
				process_Bionano(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "*.smap file" << std::endl;
				std::cerr << "output file" << std::endl;
			}
			break;
		case 5:
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
		case 6: //SURVIVOR_ant??
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
			/*if (argc == 8) {
			 //filter SV calls delly
			 //	filter_vcf(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), std::string(argv[7]));
			 } else {
			 std::cerr << "VCF file to filter" << std::endl;
			 std::cerr << "BED file with regions to ignore (NA to disable)" << std::endl;
			 std::cerr << "Min SV size (-1 to disable)" << std::endl;
			 std::cerr << "Max SV size (-1 to disable)" << std::endl;
			 std::cerr << "Min num of reads (-1 to disable)" << std::endl;
			 std::cerr << "Output vcf file" << std::endl;
			 }*/
			break;
		case 8:
			if (argc == 6) {
				parse_VCF_to_bed(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "vcf file" << std::endl;
				std::cerr << "min size" << std::endl;
				std::cerr << "max size" << std::endl;
				std::cerr << "output file" << std::endl;
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
			if (argc == 7) {
				summary_SV(std::string(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), std::string(argv[6]));
				std::cout << "You can find an R script in the src/R-scripts/ to create plots given the summary output files." << std::endl;
			} else if (argc == 5) {
				summary_SV_stream(atoi(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "vcf file" << std::endl;
				std::cerr << "Min SV size (disable: -1)" << std::endl;
				std::cerr << "Max SV size (disable: -1)" << std::endl;
				std::cerr << "Min number read support (disable: -1)" << std::endl;
				std::cerr << "output summary file" << std::endl;
			}
			break;
		case 14:
			if (argc == 7) {
				//		summary_MT(std::string(argv[2]), std::string(argv[6]));
			} else {
				std::cerr << "Parsed sam file" << std::endl;
				std::cerr << "output summary file" << std::endl;
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
				std::cerr << "input file merged vcf file" << std::endl;
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
		case 28: //not needed
			if (argc == 6) {
				prepare_svviz(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "input vcf file" << std::endl;
				std::cerr << "input bam file" << std::endl;
				std::cerr << "input ref file" << std::endl;
				std::cerr << "output svviz file" << std::endl;
			}
			break;
		case 29:
			// Make matrix for Y: sample X SV (varianten)
			//Run through window (X achse) count # times same vector is observed.
			//a1=#unique. a2=#2 times observed,.... jede zeile = 1 window.
			if (argc == 4) {
				//summarize_svs_table_window(std::string(argv[2]), atoi(argv[3]), std::string(argv[4]));
				summarize_svs_table_window_stream(atoi(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "input vcf sumary file" << std::endl;
				std::cerr << "window size" << std::endl;
				std::cerr << "output  file" << std::endl;
			}
			break;
		case 30:
			if (argc == 6) {
				generate_dist_mat(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "input vcf SVs file" << std::endl;
				std::cerr << "input vcf SNP file" << std::endl;
				std::cerr << "input weighted file" << std::endl;
				std::cerr << "output  file" << std::endl;
			}
			break;

		case 31:
			if (argc == 8) {
				overlap_snpsGWASDB(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), std::string(argv[7]));
				//overlap_snps(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), atoi(argv[5]),atoi(argv[6]), std::string(argv[7]));
			} else {
				std::cerr << "input vcf SVs file" << std::endl;
				std::cerr << "input vcf SNP file" << std::endl;
				std::cerr << "max distance" << std::endl;
				std::cerr << "min SV length" << std::endl;
				std::cerr << "min AF (0-100)" << std::endl;
				std::cerr << "output  file" << std::endl;
			}
			break;
		case 32:
			if (argc == 7) {
				overlap_snps_gwas(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), atoi(argv[5]), std::string(argv[6]));
			} else {
				std::cerr << "input vcf SVs file" << std::endl;
				std::cerr << "input random vcf SVs file" << std::endl;
				std::cerr << "max distance" << std::endl;
				std::cerr << "min SV length" << std::endl;
				std::cerr << "output  file" << std::endl;
			}
			break;
		case 33:
			if (argc == 5) {
				convert_mummer_svs(std::string(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "input SVs MUMmer file (Show-diff)" << std::endl;
				std::cerr << "min SV length" << std::endl;
				std::cerr << "output vcf file" << std::endl;
			}
			break;
		case 34:
			//if(argc == 4){
			if (argc == 7) {
				est_cov(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
			} else {
				std::cerr << "Read length" << std::endl;
				std::cerr << "Number of SVs" << std::endl;
				std::cerr << "Minimum overlap of read (bp)" << std::endl;
				std::cerr << "Minimum required overlapp of reads (support)" << std::endl;
				std::cerr << "Targeted genome coverage" << std::endl;
			}
			//	}
			break;
		case 35:
			if (argc == 6) {
				simplify_svs(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "input SURVIVOR_ant vcf file" << std::endl;
				std::cerr << "input population file" << std::endl;
				std::cerr << "min SV length" << std::endl;
				std::cerr << "output table file" << std::endl;
			}
			break;
		case 36:
			if (argc == 6) {
				overlapp_mummer(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "SVs VCF file" << std::endl;
				std::cerr << "list of mummer files" << std::endl;
				std::cerr << "Allowed distance" << std::endl;
				std::cerr << "output vcf file" << std::endl;
			}
			break;
		case 37:
			if (argc == 6) {
				generate_random_regions(std::string(argv[2]), std::string(argv[3]), atoi(argv[4]), std::string(argv[5]));
			} else {
				std::cerr << "Genome fasta file" << std::endl;
				std::cerr << "SVs VCF file" << std::endl;
				std::cerr << "min SVs size" << std::endl;
				std::cerr << "output bed file" << std::endl;
			}
			break;
		case 38:
			if (argc == 4) {
				trans_vcf(std::string(argv[2]), std::string(argv[3]));
			} else {
				std::cerr << "SVs VCF file" << std::endl;
				std::cerr << "SVs VCF file" << std::endl;
			}
			break;
		default:
			break;

		}
	} else {

		//split in : manipulating VCF, Simulation, merging,
		std::cerr << "Possible options" << std::endl;
		std::cerr << "1: Simulate SV on genome" << std::endl;
		std::cerr << "2: Simulate PacBio/ONT reads" << std::endl;
		std::cerr << "3: Evaluate SV calling based on simulation option 1" << std::endl;
		std::cerr << "4: Convert Bionano smap to vcf file " << std::endl;
		std::cerr << "5: Consensus call from multipe SV vcf files" << std::endl;
		std::cerr << "6: Extract genes influenced by SVs" << std::endl;
		std::cerr << "7: Filter SV calls" << std::endl;
		std::cerr << "8: Converts vcf to bedpe" << std::endl;
		std::cerr << "9: Filter and convert SV calls from BED files" << std::endl;
		std::cerr << "10: Convert SV calls from PBHoney (tails)" << std::endl;
		std::cerr << "11: Convert SV calls from Assemblytics" << std::endl;
		std::cerr << "12: Summarize MQ 0 coverage to bed file" << std::endl;
		std::cerr << "13: Summarize SVs events in VCF file" << std::endl;
		/*	std::cerr << "15: Fast eval for Sniffles paper" << std::endl;
		 std::cerr << "16: Convert Bionano (.smap) to vcf" << std::endl;
		 std::cerr << "17: Convert Bionano (.smap) to vcf" << std::endl;
		 std::cerr << "18: Convert vcf SV calls to bed" << std::endl;
		 std::cerr << "19: Pairwise overlap of large SV merged calls (option 5)" << std::endl;
		 std::cerr << "20: Summarize MQ 0 coverage to bed file" << std::endl;
		 std::cerr << "21: NA" << std::endl;
		 std::cerr << "22: NA" << std::endl;
		 std::cerr << "23: Summarize sambada coverage file for filtering" << std::endl;
		 std::cerr << "24: Summarize GiaB input call sets" << std::endl;
		 std::cerr << "25: Summarize GiaB merged calls" << std::endl;
		 std::cerr << "26: VCF -> BED" << std::endl;
		 std::cerr << "27: Detect potential nested/adjecent SVs based on Sniffles" << std::endl;
		 std::cerr << "28: NA" << std::endl;
		 std::cerr << "29: Summarize SVs events per window based on multi sample file (option 5)" << std::endl;
		 std::cerr << "30: Generate pairwise distance matrix of SVs merged file (option 5)" << std::endl;
		 std::cerr << "31: NA" << std::endl;
		 std::cerr << "32: NA" << std::endl;
		 std::cerr << "33: Convert MUMmer diff to vcf" << std::endl;
		 std::cerr << "34: Estimate ability to detect SVs" << std::endl;*/
	}
	return 0;
}
