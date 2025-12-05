#  A repository of useful bioinformatics notes and codes
## Helpful notes and tutorials
- Introduction to bash and cheat sheets: https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/Bash_lecture.pdf
- Tutorials for analysis of low coverage genome data: https://github.com/nt246/lcwgs-guide-tutorial/tree/main
  - data processing, genotype snp calling, linkage disequilibrium, population structure, popgen summary stats

## Bioinformatics codes
Here is where I keep bit and pieces of codes and programs that I've tried as preliminary analysis or side projects.
#### [Short read genome assembly](https://github.com/huiqingyeooo/genomics/tree/main/genome_assembly)
Side project comparing various short read assemblers (DiscovarDeNovo, megahit, GATB minia, and spades) using both whole genome resequencing reads and anchor hybrid enrichment data. I also tested if deduplication of reads improves the assemblies. Quality of the assemblies were assessed by looking at stats such as L50, N50 (stats.sh from bbmap) and completeness of orthologs recovered (BUSCO, compleasm).

#### [Assign taxonomies and filter reads](https://github.com/huiqingyeooo/genomics/tree/main/BASTA)
Used BASTA to assign taxonomies to sequences or groups of sequences based on the Last Common Ancestor (LCA) of a number of best hits. Reads can then be filtered by BASTA annotations using filter_basta_fasta.py.

#### [Extract genes from genomic reads](https://github.com/huiqingyeooo/genomics/tree/main/extract_genes)
1. Methods for assembling and extracting genes/loci from NGS data (fq.gz) using a reference database, without the need for genome assembly or alignment. <br>
Useful for looking at genes of interest (e.g., obtaining barcodes or constructing multi-loci phylogenies). <br>
Programs: [aTRAM](https://github.com/juliema/aTRAM) and [GeneMiner2](https://github.com/sculab/GeneMiner2/blob/master/manual/EN_US/command_line.md)

2. Extract sequences from regions of interest<br>
Align reads to reference genome, call variants and create consensus fasta files (allowing for IUPAC ambiguity codes). Identify regions containing genes of interest from gff files, and use coordinates to extract the relevant regions from the fasta files.

#### [Align scaffolds to a reference genome and assess synteny](https://github.com/huiqingyeooo/genomics/tree/main/assembly_scaffolding)
Using RagTag scaffold to align scaffolds to a specified reference genome, and assessing synteny with syri. Useful when the organism that you are working on does not have a chromosomal-level genome assembly.

#### [Compare allele concordance by site between two vcf files](https://github.com/huiqingyeooo/genomics/tree/main/allele_concordance)
Compares allele concordance between two vcf files, has the flexibility to score allele matches even if the genotypes do not fully match (e.g. 0/1 from vcf1 and 0/0 from vcf2 will be scored as a partial match. The codes also extracts read depth and scores heterozygosity by site.
