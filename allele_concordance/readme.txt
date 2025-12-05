PATH=/work/soghigian_lab/apps/bcftools-1.16/:$PATH

### Extract chromosome, position, GT and depth (DP) from vcf file
# First, extract the sample names
bcftools query -l b9.vcf.gz > b9_samples_names.txt
sed -i 's#/scratch/36471977/2_indelRealigned/##g;s#/scratch/36602683/2_indelRealigned/##g;s/.bam//g' b9_samples_names.txt

# Then, create the header line with "CHROM POS" followed by sample names
(echo -e "CHROM\tPOS\t$(paste -sd '\t' b9_samples_names.txt)") > b9_genotypes_depths.txt

# Finally, append the GT:DP data for each sample
bcftools query -f '%CHROM\t%POS[\t%GT:%DP]\n' b9.vcf.gz >> b9_genotypes_depths.txt

%CHROM — chromosome column
%POS — position column
[\t%GT:%DP] — for each sample, print a tab, then GT:DP
\n — new line for each variant

cp /scratch/37181873/3_angsd/b9_realigned_minQ10_minDepth1/b9_genotypes_depths.txt /scratch/37181873/8_site_concordance/b9_nanopore.txt
cp /scratch/37181873/3_angsd_illumina/b9_realigned_minQ10_minDepth1/b9_genotypes_depths.txt /scratch/37181873/8_site_concordance/b9_illumina.txt

# make sure the sample names in illumina and nanopore files match
while read -r pattern replacement; do
sed -i "s/$pattern/$replacement/" b9_nanopore.txt
done < rename.txt

#### Compare illumina and nanopore files
I have two files illumina_genotypes_depths.txt and nanopore_genotypes_depths.txt.

illumina_genotypes_depths.txt looks like this:
CHROM   POS     barcode10_q20_realigned barcode11_q20_realigned barcode12_q20_realigned barcode13_q20_realigned barcode14_q20_realigned barcode15_q20_realigned barcode16_q20_realigned barcode17_q20_realigned barcode18_q20_realigned barcode19_q20_realigned
chromosome_2    16581   0/0:4   0/0:2   ./.:0   0/0:5   0/0:33  0/1:7   0/1:3   ./.:1   0/0:4   0/0:1
chromosome_2    16590   1/1:4   0/0:2   ./.:0   0/1:4   0/0:36  0/0:7   0/0:4   0/0:1   0/0:5   0/1:1
chromosome_2    16643   0/0:4   0/1:2   ./.:0   0/1:5   0/0:36  0/0:7   0/0:4   ./.:0   1/1:4   0/0:1
chromosome_2    16675   1/1:3   0/0:2   ./.:0   0/1:4   0/0:37  0/0:7   0/0:4   0/0:2   0/0:5   0/1:1

and nanopore_genotypes_depths.txt looks like this: 
CHROM   POS     barcode10_q20_realigned barcode11_q20_realigned barcode12_q20_realigned barcode13_q20_realigned barcode14_q20_realigned barcode15_q20_realigned barcode16_q20_realigned barcode17_q20_realigned barcode18_q20_realigned barcode19_q20_realigned
chromosome_2    16581   0/1:2   0/0:2   ./.:0   0/0:5   0/0:33  0/1:7   0/1:3   ./.:1   0/0:4   0/0:1
chromosome_2    16590   ./.:1   0/0:2   ./.:0   0/1:4   0/0:36  0/0:7   0/0:4   0/0:1   0/0:5   0/1:1
chromosome_2    16642   0/0:5   0/1:2   ./.:0   0/1:5   0/0:36  0/0:7   0/0:4   ./.:0   1/1:4   0/0:1
chromosome_2    16676   0/0:3   0/0:2   ./.:0   0/1:4   0/0:37  0/0:7   0/0:4   0/0:2   0/0:5   0/1:1

I want to compare both files to see if the genotypes at each position for each sample between the two files are concordant or not.

0/0 and 0/0 --> match
1/1 and 1/1 --> match
1/0 and 1/0 --> match
0/1 and 0/1 --> match
0/1 and 1/0 --> match
1/0 and 0/1 --> match

0/1 and 1/1 --> partial_match
0/1 and 0/0 --> partial_match

0/0 and 1/1 --> no_match

if ./. is found in the illumina_genotypes_depths.txt file for that sample and position, return missing_illumina
if ./. is found in the nanopore_genotypes_depths.txt file for that sample and position, return missing_nanopore
if ./. is found in both .txt files, return missing_both

if 1/0 or 0/1 is found in the illumina_genotypes_depths.txt file for that sample and position, return yes under the column HET_illumina
if 1/0 or 0/1 is found in the nanopore_genotypes_depths.txt file for that sample and position, return yes under the column HET_nanopore

Compare the same position (POS) and sample for both files and product a text file:
CHROM	POS	SAMPLE	TYPE	DEPTH_illumina	DEPTH_nanopore	HET_illumina	HET_nanopore
chromosome_2	16581	barcode10_q20_realigned	partial_match	4	2	no	yes
chromosome_2	16590	barcode10_q20_realigned	missing_nanopore	4	1	no	no
chromosome_2	16642	barcode10_q20_realigned	match	4	5	no	no
chromosome_2	16675	barcode10_q20_realigned	no_match	3	3	no	no

# Compare alleles by running either compare_alleles.awk or compare_alleles.py
awk -f compare_alleles.awk b9_illumina.txt b9_nanopore.txt
python compare_alleles.py b9_illumina.txt b9_nanopore.txt b9_allele_comparison.csv b9_summary.csv

# compare_alleles_with_annotation.awk has annotations but doesn't work when i run it so use compare_alleles.awk instead

# Produces two files: One is a allele comparison text file specifying depth, heterozygosity and type of match for each SNP position
# Compare the same position (POS) and sample for both files and product a text file:
# CHROM	POS	SAMPLE	TYPE	DEPTH_illumina	DEPTH_nanopore	HET_illumina	HET_nanopore
# chromosome_2	16581	barcode10_q20_realigned	partial_match	4	2	no	yes
# chromosome_2	16590	barcode10_q20_realigned	missing_nanopore	4	1	no	no
# chromosome_2	16642	barcode10_q20_realigned	match	4	5	no	no
# chromosome_2	16675	barcode10_q20_realigned	no_match	3	3	no	no

# The other file is a summary file for each individual sample
# SAMPLE  match   partial_match   no_match        missing_illumina        missing_nanopore        missing_both
# 269CA_q20_realigned     4947663 1295821 177951  8440520 5210051 2445700
# 016CA_q20_realigned     5968473 1182438 211171  10222232        4014694 918698
# 036CA_q20_realigned     5794611 1231570 204793  9636372 4444386 1205974

# Subset genotype_comparison file to test out plotting in R
# !! remember to change file name to allele_comparison in future use
head -n 1 genotype_comparison.txt > header.txt
sed '1d' genotype_comparison.txt | shuf -n 5000 > random_subset_no_header.txt
cat header.txt random_subset_no_header.txt > genotype_comparison_subset.txt

local=/Users/huiqing/Documents/Cx_pipiens/nanopore
arc=/scratch/37181873/8_site_concordance/genotype_comparison_subset.txt
scp huiqing.yeo@arc.ucalgary.ca:${arc} ${local}

head -n 1 genotype_comparison_filtered.csv > header.csv
sed '1d' genotype_comparison_filtered.csv | shuf -n 10000 > random_subset_no_header.csv
cat header.csv random_subset_no_header.csv > genotype_comparison_filtered_subset.csv
rm random_subset_no_header.csv header.csv

local=/Users/huiqing/Documents/Cx_pipiens/nanopore
arc=/scratch/37181873/8_site_concordance/genotype_comparison_filtered_subset.csv
scp huiqing.yeo@arc.ucalgary.ca:${arc} ${local}
