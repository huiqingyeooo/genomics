#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=60:00:00
#SBATCH --mem=150G

####### Set environment variables ###############
module load gcc/10.2.0
#module load bioconda/conda3
source ~/software/init-conda-sl
module load perl/5.32.0

PATH=/work/soghigian_lab/apps/bbmap:$PATH
PATH=/work/soghigian_lab/apps/SPAdes-3.15.3-Linux/bin:$PATH
PATH=/work/soghigian_lab/apps/fastp:$PATH
PATH=/work/soghigian_lab/apps/mafft-7.490-without-extensions/bin:$PATH
PATH=/work/soghigian_lab/apps/sqlite-amalgamation-3370000:$PATH
PATH=/work/soghigian_lab/apps/hmmer-3.3.2/bin:$PATH
PATH=/work/soghigian_lab/apps/exonerate-2.2.0-x86_64/bin:$PATH
PATH=/work/soghigian_lab/apps/ncbi-blast-2.12.0+/bin:$PATH
PATH=/work/soghigian_lab/apps/Orthograph-0.7.1:$PATH
export PERL5LIB=/work/soghigian_lab/apps/conda/envs/ahtree/lib/site_perl/5.26.2:/work/soghigian_lab/apps/conda/envs/ahtree/lib/site_perl/5.26.2/x86_64-linux-thread-multi

conda activate /work/soghigian_lab/apps/miniconda3/envs/busco_env
module load lib/boost/1.70.0-mpi
module load lib/zlib/1.2.11
export PATH=/work/soghigian_lab/apps/metaeuk/bin/:$PATH
export PATH=/work/soghigian_lab/apps/R-4.1.2/bin:$PATH
export PATH="/work/soghigian_lab/apps/Augustus/bin:$PATH"
export PATH="/work/soghigian_lab/apps/Augustus/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/work/soghigian_lab/apps/Augustus/config"
export PATH="/work/soghigian_lab/apps/ncbi-blast-2.12.0+/bin:$PATH"
export PATH="/work/soghigian_lab/apps/hmmer-3.3.2/bin:$PATH"

echo "starting run at: `date`"

########## Prepping raw reads ########## 
#raw_reads=/work/soghigian_lab/data/illumina/singapore/
#list="Ae_annandalei_E1
#Ae_gardnerii_E4
#Ae_malayensis_E2
#Ae_riversi_J1_1"

#for file in $list
#do
#cp ${raw_reads}/${file}_1.fq.gz /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/reads/
#cp ${raw_reads}/${file}_2.fq.gz /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/reads/
#done

#raw_reads=/work/soghigian_lab/data/ahe/reads
#list="W4_S194
#W10_S235
#W8_S206"

#for file in $list
#do
#base=${file%_*}
#cp ${raw_reads}/${file}_L002_R1_001.fastq.gz /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/reads/${base}_1.fq.gz
#cp ${raw_reads}/${file}_L002_R2_001.fastq.gz /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/reads/${base}_2.fq.gz
#done

#raw_reads=/work/soghigian_lab/data/ahe/rawReadsforTy
#list="NCS57_S19"
#for file in $list
#do
#base=${file%_*}
#cp ${raw_reads}/${file}_L002_R1_001.fastq.gz /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/reads/${base}_1.fq.gz
#cp ${raw_reads}/${file}_L002_R2_001.fastq.gz /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/reads/${base}_2.fq.gz
#done

#Set up the script to run

#this script will dedupe with clumpify, use fastp to clean and trim adaptors/low quality sequences, and complete an assembly.
#the output html files of fastp are stored in the current directory.
#the current directory must have a folder called reads, containing files in the format of SampleName_R1.fastq.gz and SampleName_whateverelse_R2.fastq.gz . If the files are in another format the script must be changed to accomodate this!
#the current directory must also have a folder called assemblies.  This folder will store the completed Spades assemblies.
#to run this script, he commented loop below, and paste it into the console with the # removed.  It will loop over reads and run the pipeline.

#for file in $(ls ./reads/*_R1.fq.gz);
#do
#taxa=$(basename $file _R1.fq.gz)
#sed "s/TEMPLATE/${taxa}/g" pipeline_wgrs_TEMPLATE.slurm > pipeline_wgrs_${taxa}.slurm
#sbatch pipeline_wgrs_${taxa}.slurm
#done

#Note that this is sent to run for 100 hours, and so for very large sequence runs (or actual transcriptomes) it may be better to do the assembly separately.

######################################
#NOTE THAT IF YOU CHANGE THE THREADS OR MAXIMUM MEMORY AS PART OF COMPUTING RESOURCES OR PART OF THE ACTUAL EXECUTION OF PROGRAMS, YOU MUST ADJUST THIS SCRIPT AS WELL!
######################################
cd /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup
taxname=`echo TEMPLATE_R1.fq.gz | awk -F"_" '{print $1}'`
read1=`ls reads/TEMPLATE_R1.fq.gz`
read2=`ls reads/TEMPLATE_R2.fq.gz`

mkdir ${taxname}
mkdir fastp_reports
mkdir assemblies
mkdir orthologs
mkdir orthologs/${taxname}
mkdir orthograph_done

#note that this will throw a minor error if this isn't the first Phase I script to run. This isn't a problem, but if you hate that fact, just remove this command and do it manually
#cd ${taxname}

#mkdir reads
#clumpify.sh in=../${read1} in2=../${read2} out=reads/${taxname}_R1_dedupe.fastq.gz out2=reads/${taxname}_R2_dedupe.fastq.gz dedupe

#fastp -i reads/${taxname}_R1_dedupe.fastq.gz -I reads/${taxname}_R2_dedupe.fastq.gz -m --merged_out reads/${taxname}_merged.fq.gz --out1 reads/${taxname}_outfile1.fq.gz --out2 reads/${taxname}_outfile2.fq.gz  --unpaired1 reads/${taxname}_unpa1.fq.gz --unpaired2 reads/${taxname}_unpa2.fq.gz --detect_adapter_for_pe --cut_front -g -x -w 8
#cat reads/${taxname}_unpa1.fq.gz reads/${taxname}_unpa2.fq.gz > reads/${taxname}_s.fq.gz

#spades.py -m 150 -t 12 -o assembly_${taxname} --pe1-1 reads/${taxname}_outfile1.fq.gz --pe1-2 reads/${taxname}_outfile2.fq.gz --pe1-m reads/${taxname}_merged.fq.gz --pe1-s reads/${taxname}_s.fq.gz 

#cp assembly_${taxname}/scaffolds.fasta ../assemblies/${taxname}.fasta
#cp fastp.html ../fastp_reports/${taxname}_fastp.html
#rm -r reads
#rm -r assembly_${taxname}

# Run busco to obtain retraining parameters for augustus
assembly=/work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/assemblies 

cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}

busco -m genome -c 12 --offline --augustus -f \
-i ${assembly}/${taxname}.fasta \
-o ${taxname} -l /work/soghigian_lab/apps/busco_lineages/diptera_odb10

cd /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/
mkdir busco
cp -r /scratch/${SLURM_JOB_ID}/${taxname}/run_diptera_odb10/augustus_output/retraining_parameters/BUSCO_${taxname}/ /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/busco/
cp /scratch/${SLURM_JOB_ID}/${taxname}/*short* /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/busco/${taxname}.short_summary.txt

#cd /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/busco/BUSCO_${taxname}
#rename BUSCO_ '' *
#sed -i 's/BUSCO_//g' ${taxname}_parameters.cfg

#export AUGUSTUS_CONFIG_PATH="/work/soghigian_lab/apps/Augustus/config"
#mkdir $AUGUSTUS_CONFIG_PATH/species/${taxname}
#cp *${taxname}* $AUGUSTUS_CONFIG_PATH/species/${taxname}/

#Run augustus to identify and extract gene regions first before using inputs for orthograph
#To run augustus without splitting
#for sample in $list
#do
#augustus --species=${sample} ${spades_fasta}${sample}_scaffolds.fasta > ${spades_fasta}${sample}_scaffolds.gff
#done

#Split the genome to parallelize augustus - check speed
#cd /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup
#mkdir augustus
#cd augustus
#mkdir split
#/work/soghigian_lab/apps/Augustus/scripts/splitMfasta.pl /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/assemblies/${taxname}.fasta --outputpath=split --minsize=10000000

#augDir=aug_out
#mkdir $augDir

#cd /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/augustus
#tot_seq=$(ls split/${taxname}.split.*.fa | wc -l)
#seq 1 $tot_seq | parallel -j 14 --bar --no-notice "nice augustus --softmasking=True --genemodel=complete --species=${taxname} split/${taxname}.split.{}.fa > $augDir/augustus.{}.out 2> $augDir/augustus.{}.err"
#cat $augDir/augustus.*.out > ${taxname}.tmp.out
#/work/soghigian_lab/apps/Augustus/scripts/join_aug_pred.pl < ${taxname}.tmp.out > ${taxname}.complete.augustus.gff
#rm ${taxname}.tmp.out

#perl /work/soghigian_lab/apps/Augustus/scripts/getAnnoFasta.pl ${taxname}.complete.augustus.gff --seqfile=/work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/assemblies/${taxname}.fasta

#rm -r $AUGUSTUS_CONFIG_PATH/species/${taxname}

# Generate folders and files for orthograph
#cd /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/${taxname}
#mkdir input
#mkdir configs
#cp /work/soghigian_lab/huiqing.yeo/genome_assembly/spades_dedup/augustus/${taxname}.complete.augustus.codingseq input/${taxname}.fasta

#mkdir culicologs
#mkdir done
#rm configs/orthograph_${taxname}.conf

#cp /work/soghigian_lab/databases/culicologs/orthograph.conf .
#cp /work/soghigian_lab/databases/culicologs/ref_tax.txt .

#sed "s/basetax/${taxname}/g" orthograph.conf > configs/orthograph_${taxname}.conf

#orthograph-analyzer --configfile configs/orthograph_${taxname}.conf
#orthograph-reporter --configfile configs/orthograph_${taxname}.conf
#perl /work/soghigian_lab/apps/Orthograph-0.7.1/summarize_orthograph_results.pl -i output -o done -c -t -s -u -d ref_tax.txt
#cp -r done/* ../orthologs/${taxname}/
#cd output
#tar -czvf ${taxname}.tar.gz ${taxname}/
#mv ${taxname}.tar.gz ../../orthograph_done/${taxname}.tar.gz

# don't run the below commands if you want to check that the files have been correctly generated and copied
#rm -r culicologs
#cd ${wd}
#rm -r ${taxname}/

echo "Job finised with exit code $? at: `date`"
