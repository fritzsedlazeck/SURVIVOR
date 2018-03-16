/*
 * Update_bam_pacbio.h
 *
 *  Created on: Mar 15, 2018
 *      Author: sedlazec
 */

#ifndef CONVERT_UPDATE_BAM_PACBIO_H_
#define CONVERT_UPDATE_BAM_PACBIO_H_

#include "../vcfs/Merge_VCF.h"
#include "../structs.h"
#include "../simulator/Eval_vcf.h"
#include <math.h>
#include <iosfwd>

void process_sam_forpacbio(std::string unmapped_sam, std::string mapped_sam, std::string output_sam);


#endif /* CONVERT_UPDATE_BAM_PACBIO_H_ */
