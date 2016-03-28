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
1: Simulate  SV on genome
2: Simulate PacBio reads
3: Evaluate SV calling
4: Merge SV calls (vcf) 
5: Consensus call from 2/3 callers
6: Extract genes influenced by SVs
7: Filter and convert SV calls from Delly
8: Filter and convert SV calls from Lumpy
9: Filter and convert SV calls from Pindel
10: Summarize MQ 0 coverage to bed file

**************************************
CONTACT:

If you have questions or encounter a problem please contact:
fritz.sedlazeck@gmail.com
