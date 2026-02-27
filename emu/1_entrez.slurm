#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=100G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/entrez

cd /scratch/39702728/ref_modules
taxidFile=taxids.txt
taxID=$(cat ${taxidFile} | awk -F"\t" '{print $1}')
taxaName=$(cat ${taxidFile} | awk -F"\t" '{print $2}')

workingFolder=/scratch/39702728/ref_modules

#env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel -j 1 --link sbatch entrez.slurm ::: ${taxID} ::: ${taxaName} 

echo "starting run at: `date`"
####### script ##################################
cd ${workingFolder}

# use just the esearch command to get an idea how many hits are found
# esearch -db nucleotide -query "txid7147[Organism:exp]"

#esearch -db nucleotide -query "txid7147[Organism:exp] \
#AND (COI[Gene] OR CO1[Gene] OR COX1[Gene] OR COXI[Gene] \
#OR cytochrome c oxidase subunit I[Gene] OR cytochrome c oxidase subunit 1[Gene] \
#OR cytochrome oxidase subunit 1[Gene] OR cytochrome oxidase subunit I[Gene] \
#OR COI[Title] OR CO1[Title] OR COX1[Title] OR COXI[Title] \
#OR cytochrome c oxidase subunit I[Title] OR cytochrome c oxidase subunit 1[Title] \
#OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title] \
#NOT complete genome[Title] \
#AND 300:2000[SLEN]" \

export NCBI_API_KEY=8d70027be2ed63ac8e322bdf53042857f608

esearch -db nucleotide -query "txid${1}[Organism:exp] \
AND (COI[Gene] OR CO1[Gene] OR COX1[Gene] OR COXI[Gene] \
OR cytochrome c oxidase subunit I[Gene] OR cytochrome c oxidase subunit 1[Gene] \
OR cytochrome oxidase subunit 1[Gene] OR cytochrome oxidase subunit I[Gene] \
OR COI[Title] OR CO1[Title] OR COX1[Title] OR COXI[Title] \
OR cytochrome c oxidase subunit I[Title] OR cytochrome c oxidase subunit 1[Title] \
OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title]) \
NOT complete genome[Title] \
AND 300:2000[SLEN]" \
| efetch -format uid  > uids_${1}_${2}.txt

echo "retrieval of uids for taxid ${1} ${2} done"

# split uids text file to avoid 'too many requests' and 'connection reset' error
mkdir split_${2}
cd split_${2}

split -d -l 10000 ${workingFolder}/uids_${1}_${2}.txt uids_${1}_${2}_chunk

echo "splitting of uids for taxid ${1} ${2} done"

for i in $(ls uids_${1}_${2}_chunk*);
do
chunkNo=`echo $i | awk -F"_" '{print $4}'`
cat ${i} | efetch -db nuccore -format fasta > sequences_${1}_${2}_${chunkNo}.fasta
echo "sequences_${1}_${2}_${chunkNo}.fasta download done"
done

echo "retrieval of sequences for taxid ${1} ${2} done"

cat sequences_${1}_${2}_*.fasta > sequences_${2}_all.fasta

echo "Job finised with exit code $? at: `date`"
