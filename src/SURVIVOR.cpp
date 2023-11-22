//============================================================================
// Name        : SURVIVOR.cpp
// Author      : Fritz Sedlazeck
// Version     :
// Copyright   : MIT license
// Description : SURVIVOR toolkit
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
#include "phasing/Phasing_vcf.h"
#include "analysis_sv/Density_VCF.h"

Parameter* Parameter::m_pInstance = NULL;

//todo: merge 1: Check all subtypes during overlap. 2: Manage split up after merge.
//todo: maybe a low mem version..?? Just storing the variant ID per input, then go back in the end..
//make file: LIBS +=-lz

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
				combine_calls_svs(std::string(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), std::string(argv[9]));
			} else {
				std::cerr << "File with VCF names and paths" << std::endl;
				std::cerr << "max distance between breakpoints (0-1 percent of length, 1- number of bp) " << std::endl;
				std::cerr << "Minimum number of supporting caller" << std::endl;
				std::cerr << "Take the type into account (1==yes, else no)" << std::endl;
				std::cerr << "Take the strands of SVs into account (1==yes, else no)" << std::endl;
				std::cerr << "Disabled." << std::endl;
				std::cerr << "Minimum size of SVs to be taken into account." << std::endl;
				std::cerr << "Output VCF filename" << std::endl;
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
				std::cerr << "Original SNP file" << std::endl;
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

		} else if (strcmp(argv[1], "genComp") == 0) {
			if (argc == 5) {
				summary_venn(std::string(argv[2]), bool(atoi(argv[3])==1), std::string(argv[4]));
			} else {
				std::cerr << "Merged Vcf file" << std::endl;
				std::cerr << "Normalize output (1==yes, else no)" << std::endl;
				std::cerr << "Output: pariwise overlap matrix" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "parent_phasing") == 0) {
			if (argc == 7) {
				parental_phasing(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]),std::string(argv[5]), std::string(argv[6]));
			} else {
				std::cerr << "Merged parental Vcf file" << std::endl;
				std::cerr << "Hapcut2 final file" << std::endl;
				std::cerr << "GATK sub file" <<std::endl;
				std::cerr << "SNP vcf propant" <<std::endl;
				std::cerr << "Output matrix" << std::endl;
			}
			exit(0);
		} else if (strcmp(argv[1], "genDensity") == 0) {
			if (argc == 5) {
				density_VCF(std::string(argv[2]), atoi(argv[3]), std::string(argv[4]));
			} else {
				std::cerr << "Population VCF file (SURVIVOR merge)" << std::endl;
				std::cerr << "Window size" << std::endl;
				std::cerr << "Output file" << std::endl;
			}
			exit(0);
		}

		/*else if (strcmp(argv[1], "scrub") == 0) {
		 if (argc == 5) {
		 //	scrup_svs(std::string(argv[2]), string(argv[3]), std::string(argv[4]));
		 } else {
		 std::cerr << "VCF file SV" << std::endl;
		 std::cerr << "VCF file SNP" << std::endl;
		 std::cerr << "Output SV vcf file" << std::endl;
		 }
		 exit(0);
		 }*/
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
	std::cerr << "\tmerge\tCompare or merge VCF files to generate a consensus or multi sample VCF files." << std::endl;
	std::cerr << "\tgenComp\tGenerates a pairwise comparison matrix based on any multi sample VCF file" << std::endl;
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
	std::cerr << "\tconvertAssemblytics\tConverts Assemblytics to a VCF file" << std::endl;

	exit(0);
}

//todo output: pb 0 -> 1

int main(int argc, char *argv[]) {

	official_interface(argc, argv);
	return 0;
}
