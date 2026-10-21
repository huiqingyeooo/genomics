# Scrapbook for code snippets
This is where I save code snippets that have been useful for me, along with notes and comments. The codes are easy enough to search or type out, but handy to have on hand to copy/paste nevertheless.

## Codes for slurm system
```bash
# Transferring files
# From local to arc
scp ${local} username@arc.ucalgary.ca:${arc}

# From arc to local
scp username@arc.ucalgary.ca:${arc} ${local}

# interactive session
salloc --mem=80G -c 1 -N 1 -n 8 -t 05:00:00

# scratch file
SCRATCH=/scratch/$SLURM_JOBID

# cancel all submitted jobs at once
squeue -u huiqing.yeo | awk '{print $1}' | xargs -n 1 scancel

# cancel specific jobs with the prefix 'index_'
squeue -u huiqing.yeo | awk '{print $1, $3}' | grep "index_" | awk '{print $1}' | xargs -n 1 scancel

# check snapshot of resource usage when job is running
arc.job-info 8436364

# to analyse resource utilization after job finishes
seff 8436364

# find out job efficiency (maximum memory utilization and wall time)
sacct -j 33677597
sacct -j 33677597 --format=JobID,MaxRSS,Elapsed
sacct -j 45920000 --format=JobName,AllocCPUS,ReqMem,MaxRSS

# to see long job names 
# change number for %.30j if want to show more than 30 characters
squeue --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R" --me

# load R
module load R/4.3.1
Rscript ${pcangsdOutput}/${angsdBatchName}/pcangsd_plot.R
```

#### Job script
Naming slurm job outputs: %j specifies job ID, %x specifies job-name

```
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --mem=80G
#SBATCH --job-name=test
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --output=slurm-angsd-b1-%j.out # this is helpful

echo "starting run at: `date`"
####### Set environment variables ###############
PATH=/work/soghigian_lab/apps/fastplong:$PATH

####### script ##################################
echo "Job finised with exit code $? at: `date`"
```

#### Summarize slurm outputs
```
for sample in $(ls slurm-*.out);
do
col1=`echo $sample`
col2=`head -n 3 ${sample} | tail -n 1`
col3=`grep "converge" ${sample}`
echo $col1 $col2 $col3 >> slurm_summary.txt
done

for i in {01..24};
do 
echo "sample${i}_q18"
echo "sample${i}_q18" >> samples_coi.txt
done
```

## Bash
```bash
# make bash file executable
chmod +x file.sh

# give permissions to folder
chmod -R 777 folder

# adds current directory to end of system's path env variable. Executables in the directory are searchable from command line
export PATH=$PATH:$(pwd)

# to list all files within a subfolder
ls -lR

# remove all files except otherwise stated (e.g., file1, file2)
rm -f !(file1.slurm|file2.slurm)

# Unzip tar.gz files
tar -xvzf 

# download files from github
wget https://raw.githubusercontent.com/xiaoming-liu/stairway-plot-v2/master/stairway_plot_v2.1.2.zip -O stairway_plot_v2.1.2.zip
```

#### Renaming files/folders
```
### Renaming files or folders ###
# Loop through each subdirectory, copy a specific file to parent directory and rename it
PARENT_DIR="/work/soghigian_lab/huiqing.yeo/geneminer"
for dir in "$PARENT_DIR"/*/; do
    # Remove trailing slash and extract subfolder name
    subfolder_name=$(basename "$dir")
    source_file="$dir/results.csv"
    dest_file="$PARENT_DIR/${subfolder_name}_results.csv"

    # Check if results.csv exists in the subdirectory
    if [[ -f "$source_file" ]]; then
        cp "$source_file" "$dest_file"
        echo "Copied $source_file to $dest_file"
    else
        echo "No results.csv in $dir"
    fi
done

# Rename files with tab-delimited text file (first col: old names, second col: new names)
awk -F'\t' 'NR > 1 { system("mv \"" $1 "\" \"" $2 "\"") }' rename.txt

# Another way to rename files with a tab-delimited text file
while read -r old_name new_name; do mv -n -- "$old_name" "$new_name"; done < rename.txt

# To rename files with string 'b1' to 'b2'
for i in *; do mv "$i" "`echo $i | sed "s/b1/b2/"`"; done
for i in *; do mv "$i" "`echo $i | sed "s#output#markDup#"`"; done

# rename files to start from ticker 1 instead of 0
ls -r | grep "OTU_" | while read file; do
    # Extract the number and base name
    num=$(echo "$file" | grep -o '[0-9]\+')
    base=$(echo "$file" | sed 's/[0-9].*//')
    extension="${file##*.}"

    # Calculate new number
    new_num=$((num + 1))
    
    # Rename
    mv "$file" "${base}${new_num}.${extension}"
done

### Renaming within a file ###
# Replace spaces within a file with underscore
tr -s '[:blank:]' '_' <file > file.formatted

# Replacing names within a file with a tab-delimited text file
cp ahetree_nt2.treefile ahetree_nt2.rename.treefile
while read -r pattern replacement; do
sed -i "s/$pattern/$replacement/" ahetree_nt2.rename.treefile
done < tipLabels.txt

# Using sed to replace filepaths within a text file
sed -i 's#/work/soghigian_lab/data/illumina/chgi/mosquitoes/plate1_lane1/##g' lane1_list.txt

```

#### Loop variations
```
list="SRR1
SRR2
SRR3"
for sample in ${list};
do
sbatch fastp_${list}.slurm
done

for file in $(cat list.txt);
do
taxa=$(basename $file _1.fastq)
taxname=`echo $taxa | awk -F"-" '{print $2}'`
sed "s/FILENAME/${taxa}/g;s/TEMPLATE/${taxname}/g" fastp_TEMPLATE.slurm > fastp_${taxa}.slurm
sbatch fastp_${taxa}.slurm
done

for sample in $(ls *_outfile1.fq.gz);
do
file=`echo $sample | awk -F"_" '{print $1}'`
sed "s/TEMPLATE/${file}/g" bowtie_pipiens_TEMPLATE.slurm > bowtie_pipiens_${file}.slurm 
sbatch bowtie_pipiens_${file}.slurm
sleep 1
done
```

#### set -e, -u, -x, -o pipefail
These settings highlight hidden errors in the script that could cause 
In particular, set -e causes bash to exit immediately if any command has a non-exit status. By default, if a line of code fails, but the last line succeeds, a successful exit code is documented. This makes it very easy to miss the error.
Detailed information here: https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425?permalink_comment_id=3799230

```bash
set -euxo pipefail
```

## GNU parallel
```
env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel --delay 1 sbatch ngsadmix_b3.slurm ::: $(seq 6) ::: $(seq 10)

seq 1 ${tot_uid_chunks} | env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel -j 2 --bar --no-notice "cat uids_chunk.{} | efetch -db nuccore -format fasta > sequences_{}.fasta"
```

## Conda

```
# activate conda environment
source ~/software/init-conda-sl

# create environment and install software
conda create --prefix /work/soghigian_lab/apps/conda/envs/emu
conda activate /work/soghigian_lab/apps/conda/envs/admixture
conda install bioconda::admixture

# delete conda env and packages within that env
conda remove -p /work/soghigian_lab/apps/miniconda3/envs/multiqc --all
```

## Tab-delimited file manipulation
```
# Split a column by ":"
# these two lines give the same output
awk -F':' 'BEGIN { OFS="\t" } {print $1, $2, $3, $4}' tmp.txt > tmp2.txt
awk -F":" '$2=$2' OFS="\t" > tmp.txt > tmp2.txt
```

## Fasta/fastq file manipulation
```
# convert multiple lines of a fasta sequence to one liner
awk '/^>/ {if (seq) print seq; print; seq=""} {if (!/^>/) seq = seq $0} END {if (seq) print seq}' b2_plinkPruned_geno10_s3.snps.tab.fasta > b2_plinkPruned_geno10_s3.snps.tab2.fasta

# remove gaps in alignment
awk '/^>/ {print; getline; gsub(/-/, "", $0); print}' concatenated2.fasta > concatenated3.fasta

# remove ? in alignment
awk '/^>/ {print; getline; gsub(/?/, "", $0); print}' concatenated3.fasta > concatenated4.fasta

# summary of fasta file
source ~/software/init-conda
seqkit stats file.fasta

# fastq to fasta
PATH=/work/soghigian_lab/apps/seqtk:$PATH
seqtk seq -a INFILE.fastq.gz > OUTFILE.fasta
```

## VCF and beagle file manipulation
```
# Subset beagle files
# beagle file with site named as "scaffold_3_xxxxx"
zcat b3_plinkPruned.beagle.gz | awk 'NR==1 {print; next} $1 ~ /^scaffold_1_/ {print}' | gzip > b3_plinkPruned_chr1.beagle.gz
zcat b3_plinkPruned.beagle.gz | awk 'NR==1 {print; next} $1 ~ /^scaffold_2_/ {print}' | gzip > b3_plinkPruned_chr2.beagle.gz
zcat b3_plinkPruned.beagle.gz | awk 'NR==1 {print; next} $1 ~ /^scaffold_3_/ {print}' | gzip > b3_plinkPruned_chr3.beagle.gz

# Subset vcf files
# subset scaffold/chromosome from vcf files
PATH=/work/soghigian_lab/apps/htslib-1.16/:$PATH
PATH=/work/soghigian_lab/apps/bcftools-1.16/:$PATH
bgzip b3_plinkPruned_geno10.vcf
tabix -p vcf b3_plinkPruned_geno10.vcf.gz
bcftools view b3_plinkPruned_geno10.vcf.gz --regions scaffold_3 -o b3_plinkPruned_geno10_s3.vcf.gz

bcftools view b2.vcf.gz --regions scaffold_1 -o b2_s1.vcf.gz
bcftools view b2.vcf.gz --regions scaffold_2 -o b2_s2.vcf.gz
bcftools view b2.vcf.gz --regions scaffold_3 -o b2_s3.vcf.gz

# grep a specific snp position in a vcf.gz file
bcftools index b7.vcf.gz
bcftools view -r scaffold_270:1834 b7.vcf.gz

# filter vcf in bcftools by minimum depth of sites across all individuals
PATH=/work/soghigian_lab/apps/bcftools-1.16/:$PATH
bcftools filter -i 'DP>=3' b8.vcf.gz -Oz -o b8_minDepth3.vcf.gz
bcftools filter -i 'DP>=5' b8.vcf.gz -Oz -o b8_minDepth5.vcf.gz
bcftools filter -i 'DP>=10' b8.vcf.gz -Oz -o b8_minDepth10.vcf.gz
```

## Databases
#### Kraken database
/work/vetmed_shared_dbs/kraken_gtdb89_genbank
