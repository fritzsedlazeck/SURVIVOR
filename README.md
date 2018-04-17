# SURVIVOR
SURVIVOR is a tool set for simulating/evaluating SVs, merging and comparing SVs within and among samples, and includes various methods to reformat or summarize SVs.

Please see our github wiki for more information (https://github.com/fritzsedlazeck/SURVIVOR/wiki ) 
**************************************
## Cite:

If you use it in your study please cite:

**Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast.**   
Jeffares, Daniel C; Jolly, Clemency; Hoti, Mimoza; Speed, Doug; Shaw, Liam; Rallis, Charalampos; Balloux, Francois; Dessimoz, Christophe; Bähler, Jürg; Sedlazeck, Fritz J.   
Nature communications, Vol. 8, 14061, 24.01.2017, p. 1-11. DOI:10.1038/NCOMMS14061

**************************************

## INSTALL:

Download SURVIVOR:
<pre>
git clone https://github.com/fritzsedlazeck/SURVIVOR.git
cd SURVIVOR/Debug
make
</pre>

**************************************

## USAGE:
```
./SURVIVOR ID
```
to see the individual parameters for each option.

choose the ID from these options:
```
Program: SURVIVOR (Tools for Structural Variations in the VCF format)
Version: 1.0.3

Usage: SURVIVOR <command> [options]

Commands:
-- Simulation/ Evaluation
	simSV	Simulates SVs and SNPs on a reference genome.
	scanreads	Obtain error profiles form mapped reads for simulation.
	simreads	Simulates long reads (Pacio or ONT).
	eval	Evaluates a VCF file after SV calling over simulated data.

-- Comparison/filtering
	merge	Compare or merge VCF files to generate a consensus or multi sample vcf files.
	filter	Filter a vcf file based on size and/or regions to ignore
	stats	Report multipe stats over a VCF file
	compMUMMer	Annotates a VCF file with the breakpoints found with MUMMer (Show-diff).

-- Conversion
	bincov	Bins coverage vector to a bed file to filter SVs in low MQ regions
	vcftobed	Converts a VCF file to a bed file
	bedtovcf	Converts a bed file to a VCF file
	smaptovcf	Converts the smap file to a VCF file (beta version)
	bedpetovcf	Converts a bedpe file ot a VCF file (beta version)
	hapcuttovcf	Converts the Hapcut2 final file to a VCF file using the original SNP file provided to Hapcut2
	convertAssemblytics	Converts Assemblytics to a VCF file```
```
**************************************
## CONTACT:

If you have questions or encounter a problem please contact:
fritz.sedlazeck@gmail.com
