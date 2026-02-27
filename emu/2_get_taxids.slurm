#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/entrez

PATH=/work/soghigian_lab/apps/seqtk/:$PATH

echo "starting run at: `date`"
####### script ##################################

cd /work/soghigian_lab/huiqing.yeo/filiaroidDiv/emu/ref_database

# concatenate coi sequences together
cat sequences_Nematoda_all.fasta sequences_Cervidae_all.fasta > ref_nem_cerv.fasta

# extract accession numbers
grep "^>" ref_nem_cerv.fasta \
| sed 's/^>//' \
| cut -d' ' -f1 \
> ref_nem_cerv.accessions

# get taxids from ncbi using list of accession numbers
split -d -l 1000 ref_nem_cerv.accessions chunk
chunk_files=$(ls chunk*)

export NCBI_API_KEY=8d70027be2ed63ac8e322bdf53042857f608

env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel -j 2 \
"epost -db nucleotide -input {} \
| esummary \
| xtract -pattern DocumentSummary -element AccessionVersion TaxId Organism > ref_nem_cerv.accessions.{}.taxid" ::: ${chunk_files}

cat ref_nem_cerv.accessions.chunk*.taxid > ref_nem_cerv.accessions.taxid

# check that number of taxids and sequences matches
num_of_taxids=$(cat ref_nem_cerv.accessions.taxid | wc -l)
num_of_sequences=$(grep ">" ref_nem_cerv.fasta | wc -l)

if [[ "$num_of_taxids" == "$num_of_sequences" ]]; then
    echo "Number of taxids = Number of sequences."
    rm *chunk*
else
	echo "Number of taxids not equals to Number of sequences. Stopping script"
	exit 1
fi

# sort taxid file to obtain the same order of accessions

TMP_FILE=$(mktemp)
while read LINE; do
    grep "$LINE" "ref_nem_cerv.accessions.taxid" >>"$TMP_FILE"
done <"ref_nem_cerv.accessions"

mv -f "$TMP_FILE" "ref_nem_cerv_sorted.accessions"

# download taxdump for nodes.dmp and names.dmp files
#conda deactivate
#conda activate /work/soghigian_lab/apps/conda/envs/taxonkit
#wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
#tar -zxvf taxdump.tar.gz -C taxdump
#conda deactivate
conda activate /work/soghigian_lab/apps/conda/envs/emu

# prepare seq2tax.map file
sed "s/ /_/g;s/\;//g;s/\,//g;s/\:/_/g;s#(##g;s#)##g;s/\./_/g;s/\-/_/g" ref_nem_cerv.fasta > ref_nem_cerv_format.fasta
grep ">" ref_nem_cerv_format.fasta > tmp1.txt
sed -i "s/>//g" tmp1.txt
cat ref_nem_cerv_sorted.accessions | awk -F"\t" '{print $2}' > tmp2.txt
paste -d'\t' tmp1.txt tmp2.txt > ref_nem_cerv_seq2tax.map
rm tmp*

# generate database
emu build-database ref_nem_cerv --sequences ref_nem_cerv_format.fasta \
--seq2tax ref_nem_cerv_seq2tax.map \
--ncbi-taxonomy /work/soghigian_lab/huiqing.yeo/filiaroidDiv/emu/ref_database/taxdump

# remove lines from seq2tax and fasta files where taxid is not found in nodes.dmp
# temporary measure
# probably mismatch between ncbi taxid and taxonomy databases for scutellonema afribrachyurus
grep -v "3686840" ref_nem_cerv_seq2tax.map > ref_nem_cerv_seq2tax_2.map
grep "3686840" ref_nem_cerv_seq2tax.map | awk -F"\t" '{print $1}' > tmp.txt
grep ">" ref_nem_cerv_format.fasta | sed "s/>//g" | grep -vxf tmp.txt > tmp_subset.txt
seqtk subseq ref_nem_cerv_format.fasta tmp_subset.txt > ref_nem_cerv_format_2.fasta


emu build-database ref_nem_cerv --sequences ref_nem_cerv_format_2.fasta \
--seq2tax ref_nem_cerv_seq2tax_2.map \
--ncbi-taxonomy /work/soghigian_lab/huiqing.yeo/filiaroidDiv/emu/ref_database/taxdump

echo "Job finised with exit code $? at: `date`"
