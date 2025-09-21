#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH --job-name=blast_TEMPLATE

####### Set environment variables ###############
source ~/software/init-conda
conda activate /work/soghigian_lab/apps/conda/envs/blast
PATH=/work/soghigian_lab/apps/bioawk/:$PATH
PATH=/work/soghigian_lab/apps/bbmap/:$PATH

#cd /work/soghigian_lab/huiqing.yeo/geneminer/blast
#list="ARG03
#ARG04
#ARG06
#NC01
#NC02
#pscil
#SP06
#SP07
#SP101
#SP102
#W66
#YF100"

#for file in $list;
#do
#sed "s/TEMPLATE/${file}/g" blast_TEMPLATE.slurm > blast_${file}.slurm
#sbatch blast_${file}.slurm
#done

######## Script ###################################
echo "starting run at: `date`"

# prep blast database
# cd /work/soghigian_lab/huiqing.yeo/geneminer/blast
# makeblastdb -in blast.fasta -dbtype nucl -out blastdb ###Make blast db

cd /work/soghigian_lab/huiqing.yeo/geneminer/blast

# Blast
blastn -query /work/soghigian_lab/data/ahe/assemblies/old_ahe/TEMPLATE*.fasta \
-db /work/soghigian_lab/huiqing.yeo/geneminer/blast/blastdb \
-out TEMPLATE_blastn -evalue 1e-5 -max_target_seqs 3 -task blastn \
-num_threads 12 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"

# Sort blast results in descending pident
sort -s -n -t$'\t' -k 3,3 -r TEMPLATE_blastn > TEMPLATE_blastn_sorted

echo "Job finised with exit code $? at: `date`"
