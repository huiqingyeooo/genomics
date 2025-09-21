#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=72:00:00
#SBATCH --mem=100G

echo "starting run at: `date`"

####### Set environment variables ###############
source ~/software/init-conda-sl
source activate /work/soghigian_lab/apps/conda/envs/ragtag
PATH=/work/soghigian_lab/apps/conda/envs/ragtag/bin:$PATH
PATH=/work/soghigian_lab/apps/minimap2-2.23_x64-linux:$PATH
module load gcc/10.2.0

PATH=/work/soghigian_lab/apps/bioawk:$PATH
PATH=/work/soghigian_lab/apps/conda/envs/syri/bin:$PATH

####### Script ##################################
#cp /work/soghigian_lab/huiqing.yeo/flycatcher/flycatcher/GCF_003031625.1_ASM303162v1_genomic.fna.gz /work/soghigian_lab/huiqing.yeo/flycatcher/refGenome_ragtag
#cp /work/soghigian_lab/huiqing.yeo/flycatcher/finch/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz /work/soghigian_lab/huiqing.yeo/flycatcher/refGenome_ragtag

cd /work/soghigian_lab/huiqing.yeo/flycatcher/refGenome_ragtag
#gzip -d GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz
#gzip -d GCF_003031625.1_ASM303162v1_genomic.fna.gz

#define variables
ref=GCF_003957565.2_bTaeGut1.4.pri_genomic.fna #zebra finch
asm=GCF_003031625.1_ASM303162v1_genomic.fna #willow flycatcher
asmname=flycatcher

#scaffold to assembly to reference genome
# -u flag: add suffix to unplaced sequence headers

#ragtag.py scaffold ${ref} ${asm} -t 25 -u

#change working directory to the ragtag_output directory
cd ./ragtag_output

#create a txt file of chromosome names from the reference (total 41 chromosomes, https://www.ncbi.nlm.nih.gov/assembly/GCF_003957565.2/?shouldredirect=false)
#grep ">" ../$ref | awk 'sub(/^>/, "")' | awk 'NR==1, NR==41 {print$1}' > ref.chr.names.txt

#loop through all fasta and remove '_RagTag' from all scaffold names
#for scaffold in `ls ragtag.*.fasta`; do
#name=$(basename -s .scaffold.fasta $scaffold)
#sed 's/_RagTag//' ${scaffold} > ${name}_rename.scaffold.fasta
#done

#rename chromosomes in both reference genome and the newly renamed reference scaffolded assembly so that they have matching names
#if need to rerun code, remember to delete previous file (or change ">>" to ">"?)
#for i in `cat ref.chr.names.txt`; do
#cat ../${ref} | bioawk -c fastx -v chr="$i" '$name==chr{print">chr"$name;print$seq}' >> reference_chromosome.fasta
#cat ragtag_rename.scaffold.fasta | bioawk -c fastx -v chr="$i" '$name==chr{print ">chr"$name; print $seq}' >> ${asmname}_qry.fasta
#done

#use minimap to align the renamed reference and reference scaffolded assembly and output as *.sam
#-ax will need to be changed depending on how distantly related the reference species and the species of interest are. if it's close (or the same) just use asm5. if more distant, try asm20
# empidonax flycatchers (Suboscine, Tyranidae) are not that closely related to zerba finches (Passerida, Estrilidae), I try with asm20.
minimap2 -ax asm20 --eqx -t 25 reference_chromosome.fasta ${asmname}_qry.fasta > ${asmname}_syri.sam

#deactivate ragtag
conda deactivate

#activate syri
conda activate /work/soghigian_lab/apps/conda/envs/syri

#create synteny map from the *.sam file
#-c is output from minimap 
#-r is the reference chromosome 
#-q is the is your assembly
#-k keeps intermediate files default is false
#-F is input file type. T=table, S=sam, B=bam defaults to T
#--prefix is the prefix for file names
#--nc is the number of cores to use for analysis
syri -c ${asmname}_syri.sam -r reference_chromosome.fasta -q ${asmname}_qry.fasta -F S --prefix ${asmname}_ --nc 25

#do this in an interactive session with the syri environment activated
#--sr is the output from syri. has extension *.out
#you can specify --sr multiple times if you have multiple *.out files you want to compare but they have to be in a specific order when you do the comparison and they have to listed in exactly the same order in genome.txt file
#-W is the width of the plot in inches
#-H is height of the plot in inchaes
#genomes.txt is a tab-separated-file with columns file and name. file is the filename of the assembly, name is whatever you want to call that assembly file in the plot
#make sure that the first line is commented out (i.e., #file)
#optionally, genome.txt can have an additional column called tags. for use see https://github.com/schneebergerlab/plotsr
#$genomes.txt
#file	name	tags
#ref.fa	ref	lw:1.5
#B.fa	B	lw:1.5
#C.fa	C	lw:1.5

#plotsr --sr flycatcher_syri.out --genomes genomes.txt -o output_plot.png -W 11 -H 8
#plotsr --sr flycatcher_syri.out --genomes genomes.txt -o output_plot.png

#remove sam file, reference and scaffold genomes as they take up space. Do not accidentally delete reference scaffold assembly.

echo "Job finised with exit code $? at: `date`"
