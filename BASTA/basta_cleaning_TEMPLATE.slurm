#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=100:00:00
#SBATCH --mem=120G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/basta

module load gcc/10.2.0
module load perl/5.32.0
PATH=$PATH:/work/soghigian_lab/apps/KrakenTools/
PATH=/work/soghigian_lab/apps/bbmap:$PATH
PATH=/work/soghigian_lab/apps/SPAdes-3.15.5-Linux/bin:$PATH
PATH=/work/soghigian_lab/apps/fastp:$PATH
PATH=/work/soghigian_lab/apps/mafft-7.490-without-extensions/bin:$PATH
PATH=/work/soghigian_lab/apps/sqlite-amalgamation-3370000:$PATH
PATH=/work/soghigian_lab/apps/hmmer-3.3.2/bin:$PATH
PATH=/work/soghigian_lab/apps/exonerate-2.2.0-x86_64/bin:$PATH
PATH=/work/soghigian_lab/apps/ncbi-blast-2.12.0+/bin:$PATH
PATH=/work/soghigian_lab/apps/Orthograph-0.7.1:$PATH

bastaInput=/work/soghigian_lab/huiqing.yeo/data/fasta_assembliesNCS
#bastaInput=/work/soghigian_lab/huiqing.yeo/data/fasta_oldAHE
#bastaInput=/work/soghigian_lab/huiqing.yeo/data/fasta_assemblies
bastaOutput=/work/soghigian_lab/huiqing.yeo/aedes_coi/28s/basta_output
jobScript=/work/soghigian_lab/huiqing.yeo/aedes_coi/28s

#cd ${bastaInput}
#for file in $(ls ./*.fasta);
#do
#taxa=$(basename $file .fasta)
#echo ${taxa}
#cd ${jobScript}
#sed "s/TEMPLATE/${taxa}/g" basta_cleaning_TEMPLATE.slurm > basta_cleaning_${taxa}.slurm
#sbatch basta_cleaning_${taxa}.slurm
#sleep 2
#done

echo "starting run at: `date`"

######## Script #######################################
cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}

# Copy taxonomy database into scratch for each sample, as program can sometimes throw errors when running parallel jobs
# https://github.com/timkahlke/BASTA/issues/24
cp -r /work/soghigian_lab/apps/BASTA/taxonomy /scratch/${SLURM_JOB_ID}/

cd ${bastaOutput}

# 1. Uses parallel to chunk fasta into blocks of 500K, and submits 5 (-j) of them each time, each using 3 (num_threads) threads on the cluster
echo "starting blast for TEMPLATE"
cat ${bastaInput}/TEMPLATE.fasta | parallel -j 5 --block 500k --recstart '>' --pipe blastn -num_threads 3 -outfmt 6 -db /work/soghigian_lab/databases/gb_nt/nt -query - > TEMPLATE.blast

# 2. Assigns taxonomies to sequences or groups of sequences based on the Last Common Ancestor (LCA) of a number of best hits
# -m: minimum number of hits a sequence must have to be assigned an LCA (default 3)
# -p: percentage of hits that are used for LCA estimation (default 100)
# -l: minimum length for good blast alignment (default 100)
# -i: minimim identity of hit to be considered good (default 80)

echo "starting basta for TEMPLATE"
basta sequence -m 1 -p 51 -l 10 -i 75 -e 0.001 TEMPLATE.blast TEMPLATE.basta gb -d /scratch/${SLURM_JOB_ID}/taxonomy

# 3. Filter a given fasta based on basta annotations
echo "extracting fasta for TEMPLATE"
/work/soghigian_lab/apps/BASTA/scripts/filter_basta_fasta.py ${bastaInput}/TEMPLATE.fasta TEMPLATE.culicidae.fasta Culicidae TEMPLATE.basta

rm -r /scratch/${SLURM_JOB_ID}/taxonomy

echo "Job finised with exit code $? at: `date`"
