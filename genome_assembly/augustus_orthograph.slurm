#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=120:00:00
#SBATCH --mem=80G

####### Set environment variables ###############
#module load gcc/10.2.0
#source ~/software/init-conda-sl
#conda activate /work/soghigian_lab/apps/miniconda3/envs/busco_env
#module load lib/boost/1.70.0-mpi
#module load lib/zlib/1.2.11
#export PATH=/work/soghigian_lab/apps/metaeuk/bin/:$PATH
#export PATH=/work/soghigian_lab/apps/R-4.1.2/bin:$PATH
#export PATH="/work/soghigian_lab/apps/Augustus/bin:$PATH"
#export PATH="/work/soghigian_lab/apps/Augustus/scripts:$PATH"
#export AUGUSTUS_CONFIG_PATH="/work/soghigian_lab/apps/Augustus/config"
#export PATH="/work/soghigian_lab/apps/ncbi-blast-2.12.0+/bin:$PATH"
#export PATH="/work/soghigian_lab/apps/hmmer-3.3.2/bin:$PATH"

cd /work/soghigian_lab/huiqing.yeo/genome_assembly/spades
#mkdir retraining_parameters
busco_folder=/work/soghigian_lab/huiqing.yeo/genome_assembly/spades/retraining_parameters/
spades_fasta=/work/soghigian_lab/huiqing.yeo/genome_assembly/spades/

echo "starting run at: `date`"

list="Ae_annandalei_E1
Ae_gardnerii_E4
Ae_malayensis_E2
Ae_riversi_J1_1
W4
W8
W10"

#for sample in $list
#do
#cp -r /scratch/26327028/${sample}_busco/run_diptera_odb10/augustus_output/retraining_parameters/BUSCO_${sample}_busco/ ${busco_folder}/
#done

#cd ${busco_folder}
#for file in * ; do
#mv -v "$file" "${file#*_}"
#done

#for sample in $list
#do
#cd ${busco_folder}${sample}
#sed -i 's/BUSCO_//g' ${sample}_parameters.cfg
#sed -i 's/_busco//g' ${sample}_parameters.cfg
#mkdir $AUGUSTUS_CONFIG_PATH/species/${sample}
#cp *${sample}* $AUGUSTUS_CONFIG_PATH/species/${sample}/
#done

#To run augustus without splitting
#for sample in $list
#do
#augustus --species=${sample} ${spades_fasta}${sample}_scaffolds.fasta > ${spades_fasta}${sample}_scaffolds.gff
#done

#Split the genome to parallelize augustus - check speed
#mkdir split
#for sample in $list
#do
#/work/soghigian_lab/apps/Augustus/scripts/splitMfasta.pl ${sample}_scaffolds.fasta --outputpath=split --minsize=10000000
#done

#augDir=aug_out
#mkdir $augDir

#for sample in $list
#do
#tot_seq=$(ls split/${sample}_scaffolds.split.*.fa | wc -l)
#seq 1 $tot_seq | parallel -j 14 --bar --no-notice "nice augustus --softmasking=True --genemodel=complete --species=${sample} split/${sample}_scaffolds.split.{}.fa > $augDir/augustus.{}.out 2> $augDir/augustus.{}.err"
#cat $augDir/augustus.*.out > tmp.out
#/work/soghigian_lab/apps/Augustus/scripts/join_aug_pred.pl < tmp.out > ${sample}.complete.augustus.gff
#rm tmp.out
#perl /work/soghigian_lab/apps/Augustus/scripts/getAnnoFasta.pl ${sample}.complete.augustus.gff --seqfile=${sample}_scaffolds.fasta
#done

#mkdir input
#mkdir configs
#rm configs/orthograph_${sample}.conf

#for sample in $list
#do
#cp ${sample}.complete.augustus.codingseq ./input/${sample}.fasta
#done

#mkdir culicologs
#mkdir done

#cp /work/soghigian_lab/databases/culicologs/orthograph.conf .
#cp /work/soghigian_lab/databases/culicologs/ref_tax.txt .

#for sample in $list
#do
#sed "s/basetax/${sample}/g" orthograph.conf > configs/orthograph_${sample}.conf
#done

export PERL5LIB=/work/soghigian_lab/apps/conda/envs/ahtree/lib/site_perl/5.26.2:/work/soghigian_lab/apps/conda/envs/ahtree/lib/site_perl/5.26.2/x86_64-linux-thread-multi
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

#for sample in $list
#do
#orthograph-analyzer --configfile configs/orthograph_${sample}.conf
#orthograph-reporter --configfile configs/orthograph_${sample}.conf
#done

cd /work/soghigian_lab/huiqing.yeo/genome_assembly/spades/
perl /work/soghigian_lab/apps/Orthograph-0.7.1/summarize_orthograph_results.pl -i output -o done -c -t -s -u -d

#for sample in $list
#do
#cp -r done/* ../orthologs/${sample}/
#cd output
#tar -czvf ${sample}.tar.gz ${sample}/
#mv ${sample}.tar.gz ../../orthograph_done/${sample}.tar.gz
#rm -r ${sample}/
#done

#cd ..
#rm -r culicologs

echo "Job finised with exit code $? at: `date`"
