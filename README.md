# SURVIVOR
Toolset for SV detection for short reads

**************************************

INSTALL:

Download SURVIVOR:
```
git clone https://github.com/fritzsedlazeck/SURVIVOR.git
```

  cd SURVIVOR/Debug
  
  make

**************************************

USAGE:
```
./SURVIVOR ID
```
to see the individual parameters for each option.

choose the ID from these options:
```
Possible options
1: Simulate SV on genome
2: Simulate PacBio reads
3: Evaluate SV calling
4: Merge SV calls (vcf) 
5: Merge + Consensus call from different callers/vcf files
6: Extract genes influenced by SVs
7: Filter and convert SV calls from Delly
8: Filter and convert SV calls from Lumpy
9: Filter and convert SV calls from Pindel
10: Convert SV calls from PBHoney (tails)
11: Convert SV calls from Assemblytics
12: Summarize MQ 0 coverage to bed file
13: Summarize SVs events in VCF file
```

**************************************
CONTACT:

If you have questions or encounter a problem please contact:
fritz.sedlazeck@gmail.com
