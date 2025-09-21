#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=50G

echo "starting run at: `date`"

####### Set environment variables ###############
source ~/software/init-conda
#conda create -p /work/soghigian_lab/apps/conda/envs/blast
conda activate /work/soghigian_lab/apps/conda/envs/blast
PATH=/work/soghigian_lab/apps/bcftools-1.16:$PATH
PATH=/work/soghigian_lab/apps/samtools-1.15:$PATH

refGenome=/work/soghigian_lab/huiqing.yeo/pipiens/refGenome/pipiens/Cpip29_nodebris_gapfill.fasta
bamFiles=/work/soghigian_lab/huiqing.yeo/pipiens/1_alignedReads/canada/bam
outputFiles=/work/soghigian_lab/huiqing.yeo/pipiens/genes

# NOTE: Run each part separately by hashing out each section

####### Script ###################################
cd ${outputFiles}

### Part 1: Identify coordinates of target genes on Cpip reference genome
# Make blast reference database (only need to do it once)
#cp /work/soghigian_lab/john.soghigian/culex_genome/dcon/purged/canu_ph_hic/Cpip29_chrs_done/Cpip29_nodebris_mask_both.complete.augustus.cdsexons .
#makeblastdb -in Cpip29_nodebris_mask_both.complete.augustus.cdsexons -dbtype nucl -out Cpip29_nodebris_mask_both.complete.augustus.cdsexons

# blast and identify gene regions mapping to vgsc and ace-1 genes
#blastn -db Cpip29_nodebris_mask_both.complete.augustus.cdsexons \
#-query target_genes.fasta -out blast_target.out -outfmt 6 -max_target_seqs 1

# find coordinates of the target genes using the gff files
#grep "g1192" Cpip29_nodebris_mask_both.complete.augustus.gff -A 10
#g1192 scaffold_1 22486973 22526976 vgsc
#grep "g24107" Cpip29_nodebris_mask_both.complete.augustus.gff -A 10
#g24107 scaffold_2 69255445 69284836 ace-1


### Part 2: Extract genes
# env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel sbatch extract_genes.slurm :::: bam.list

# apply quality filter of MAPQ>20 to bam files
cd ${bamFiles}
samtools view -q 20 ${1}_sorted.bam -o ${1}_sorted_q20.bam
samtools index ${1}_sorted_q20.bam

# call variants in bam files and generate vcf files
cd ${outputFiles}
bcftools mpileup -f ${refGenome} --regions scaffold_1:22486973-22526976 ${bamFiles}/${1}_sorted_q20.bam | \
  bcftools call -mv -Oz -o ${outputFiles}/vgsc_${1}.vcf.gz
  
bcftools mpileup -f ${refGenome} --regions scaffold_2:69255445-69284836 ${bamFiles}/${1}_sorted_q20.bam | \
  bcftools call -mv -Oz -o ${outputFiles}/ace1_${1}.vcf.gz

# index vcf files
bcftools index vgsc_${1}.vcf.gz
bcftools index ace1_${1}.vcf.gz

# extract target regions
# NOTE: samtools faidx command needed for specifying gene regions. Otherwise bcftools consensus will suffice
samtools faidx ${refGenome} scaffold_1:22486973-22526976 | bcftools consensus -H I ${outputFiles}/vgsc_${1}.vcf.gz > ${outputFiles}/vgsc_${1}.fasta
samtools faidx ${refGenome} scaffold_2:69255445-69284836 | bcftools consensus -H I ${outputFiles}/ace1_${1}.vcf.gz > ${outputFiles}/ace1_${1}.fasta

# housekeeping, delete intermediate files
rm ${bamFiles}/${1}_sorted_q20.bam.bai
rm ${bamFiles}/${1}_sorted_q20.bam
rm ${outputFiles}/vgsc_${1}.vcf.gz.csi
rm ${outputFiles}/ace1_${1}.vcf.gz.csi
#rm ${outputFiles}/vgsc_${1}.vcf.gz
#rm ${outputFiles}/ace1_${1}.vcf.gz

### Part 3: Rename headers and concatenate fasta
for f in vgsc_*.fasta; do
#get base name without extension
base=`echo $f | awk -F"." '{print $1}'`
#add file name to start of sequence header
sed -i "s/^>/>${base%.faa}_/g" "${f}";
done

for f in ace1_*.fasta; do
#get base name without extension
base=`echo $f | awk -F"." '{print $1}'`
#add file name to start of sequence header
sed -i "s/^>/>${base%.faa}_/g" "${f}";
done

#concatenate individual fasta files
cat vgsc_*.fasta > vgsc.fasta
cat ace1_*.fasta > ace1.fasta
#rm vgsc_*.fasta ace1_*.fasta

#remove scaffold names (optional)
sed -i "s/_scaffold_1:22486973-22526976//g" vgsc.fasta
sed -i "s/_scaffold_2:69255445-69284836//g" ace1.fasta

echo "Job finised with exit code $? at: `date`"
