#  A repository of useful bioinformatics notes and codes
## Helpful notes and tutorials
- Introduction to bash and cheat sheets: https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/Bash_lecture.pdf
- Tutorials for analysis of low coverage genome data: https://github.com/nt246/lcwgs-guide-tutorial/tree/main
  - data processing, genotype snp calling, linkage disequilibrium, population structure, popgen summary stats

## Bits and pieces
Here is where I keep bit and pieces of codes and programs that I've tried as preliminary analysis or side projects.
#### Short read genome assembly
Side project comparing various short read assemblers (DiscovarDeNovo, megahit, GATB minia, and spades) using both whole genome resequencing reads and anchor hybrid enrichment data. I also tested if deduplication of reads improves the assemblies. Quality of the assemblies were assessed by looking at stats such as L50, N50 (stats.sh from bbmap) and completeness of orthologs recovered (BUSCO, compleasm).
