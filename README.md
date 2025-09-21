#  A repository of useful bioinformatics notes and codes
## Helpful notes and tutorials
- Introduction to bash and cheat sheets: https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/Bash_lecture.pdf
- Tutorials for analysis of low coverage genome data: https://github.com/nt246/lcwgs-guide-tutorial/tree/main
  - data processing, genotype snp calling, linkage disequilibrium, population structure, popgen summary stats

## Bits and pieces
Here is where I keep bit and pieces of codes and programs that I've tried as preliminary analysis or side projects.
#### [Short read genome assembly](https://github.com/huiqingyeooo/genomics/genome_assembly)
Side project comparing various short read assemblers (DiscovarDeNovo, megahit, GATB minia, and spades) using both whole genome resequencing reads and anchor hybrid enrichment data. I also tested if deduplication of reads improves the assemblies. Quality of the assemblies were assessed by looking at stats such as L50, N50 (stats.sh from bbmap) and completeness of orthologs recovered (BUSCO, compleasm).

#### [Assign taxonomies and filter reads](https://github.com/huiqingyeooo/genomics/BASTA)
Used BASTA to assign taxonomies to sequences or groups of sequences based on the Last Common Ancestor (LCA) of a number of best hits. Reads can then be filtered by BASTA annotations using filter_basta_fasta.py.

#### [Extract genes from genomic reads](https://github.com/huiqingyeooo/genomics/extract_genes)
Methods for assembling and extracting genes/loci from NGS data (fq.gz) using a reference database, without the need for genome assembly or alignment. Useful for looking at genes of interest (e.g., obtaining barcodes or constructing multi-loci phylogenies)
[aTRAM](https://github.com/juliema/aTRAM) and [GeneMiner2](https://github.com/sculab/GeneMiner2/blob/master/manual/EN_US/command_line.md)
