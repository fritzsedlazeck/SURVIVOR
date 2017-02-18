# SURVIVOR
SURVIVOR is a tool set for simulating/evaluating SVs, merging and comparing SVs within and among samples, and includes various methods to reformat or summarize SVs.

Please see our github wiki for more information (https://github.com/fritzsedlazeck/SURVIVOR/wiki ) 
**************************************
## Cite:

If you use it in your study please cite:

**Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast.**   
Jeffares, Daniel C; Jolly, Clemency; Hoti, Mimoza; Speed, Doug; Shaw, Liam; Rallis, Charalampos; Balloux, Francois; Dessimoz, Christophe; Bähler, Jürg; Sedlazeck, Fritz J. Sedlazeck   
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
## CONTACT:

If you have questions or encounter a problem please contact:
fritz.sedlazeck@gmail.com
