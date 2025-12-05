#!/usr/bin/awk -f
BEGIN {
    FS = "[ \t]+"; OFS = "\t"
    # FS mean field separator. "[\t]+" means one or more spaces or tab characters. each input line will be split using whitespace
    # OFS output field separator. When awk prints multiple fields, it will separate with a tab
    
    out_detail = "genotype_comparison.txt"
    # out_detail is a variable that stores the filename for detailed output
    
    out_summary = "genotype_summary.txt"
    # defines a variable that stores filename for summary output
    
    print "CHROM", "POS", "SAMPLE", "TYPE", "DEPTH_illumina", "DEPTH_nanopore", "HET_illumina", "HET_nanopore" > out_detail
    # prints a header line to the out_summary file
}

function trim(x){sub(/^[ \t]+/,"",x);sub(/[ \t]+$/,"",x);return x}
# defines a custom AWK function named trim
# sub(/^[ \t]+/,"",x) removes leading spaces/tabs
# sub(/[ \t]+$/,"",x) removes trailing spaces/tabs
# returns a clean string

# ---- First file: illumina ----
# NR refers to the total line number across all files
# FNR refers to per-file line numner (resets to 1 for each file)
# NR == FNR --> we are still in the first file
# FNR == 1 && NR == FNR --> first line of the first file
FNR==1 && NR==FNR {
    if ($1=="CHROM"){for(i=3;i<=NF;i++){sample[$i]=i; all_samples[$i]=1} next}
    # the first file's headers are: CHROM POS sample1 sample2 sample3
	# saves which column index belongs to each sample name, sample["sample1"] = 3
	# also builts all_samples, which is a combined set of all sample names across both illumina and nanopore files
}
NR==FNR {
    # runs on every line of the first file except the header
    key=$1":"$2
    # creates a genomic coordinate key, CHROM:POS
   
    for(s in sample){illum[key,s]=$(sample[s])}
    # stores illumina genotype data into 2 dimensional array by genomic coordinate and sample name
    # illum["chr1:10500","sample1"]= genotype value
    
    seen_illum[key]=1
    # tracks which genomic locations appear in the illumina file
    next
}


# ---- Second file: nanopore ----
FNR==1 {
	# when awk begins the second file, FNR resets to 1, but NR continues increasing
	# FNR==1 runs only on the header line of the nanopore file
	
    for(i=3;i<=NF;i++){
    sample2[$i]=i; # sample --> column index for nanopore
    all_samples[$i]=1 # merge all samples between files
    }
    next
    
    # sample[] --> column positions in illumina file
    # sample2[] --> column positions in nanopore file
    # all_samples[] --> union of all sample names
}

{
    key=$1":"$2
    chrom=$1; pos=$2
    # genomic coordinates "chr1:10500", chrom=chr1, pos=10500
    for(s in all_samples){
    	# retrieve illumina and nanopore genotype values
        illum_val=illum[key,s]; nano_val=$(sample2[s])
        
        # normalize missing values to ensure consistent format
        if(illum_val==""||illum_val==".") illum_val="./.:NA"
        if(nano_val==""||nano_val==".") nano_val="./.:NA"
         
        split(illum_val,a,":"); split(nano_val,b,":")
        gtI=a[1]; dpI=a[2]; gtN=b[1]; dpN=b[2]
        # illum_val = "0/1:23" --> gtI="0/1", dpI=23
        # nano_val = "0/0: 18" --> gtN "0/0", dpN=18
        
        hetI=(gtI=="0/1"||gtI=="1/0")?"yes":"no"
        hetN=(gtN=="0/1"||gtN=="1/0")?"yes":"no"
        # classify heterozygosity

        # classify
        if(gtI=="./." && gtN=="./.") type="missing_both"
        # if both missing
        
        else if(gtI=="./.") type="missing_illumina"
        # only illumina missing
        
        else if(gtN=="./.") type="missing_nanopore"
        # only nanopore missing
        
        else if((gtI==gtN) || (gtI=="0/1"&&gtN=="1/0") || (gtI=="1/0"&&gtN=="0/1")) type="match"
        # exact genotype matches including reverse heterozygosity
        
        else if(((gtI~/[01]\/[01]/)&&(gtN~/[01]\/[01]/)) && ((gtI=="0/1"&&gtN=="0/0")||(gtI=="0/1"&&gtN=="1/1")||(gtI=="1/0"&&gtN=="0/0")||(gtI=="1/0"&&gtN=="1/1")||(gtI=="1/1"&&gtN~/0\/1|1\/0/)||(gtI=="0/0"&&gtN~/0\/1|1\/0/))) type="partial_match"
        # partial match when one genotype is heterozygous, and the other is homozygous
        
        else type="no_match"
        # otherwise there is no match

        count[s,type]++
        # count classification for summary
        
        print chrom,pos,s,type,dpI,dpN,hetI,hetN >> out_detail
    }
    seen[key]=1
    # track site as seen
}

# ---- After reading nanopore, handle illumina-only SNPs ----
END {
	# seen_illum[k] --> SNPs seen in illumina
	# seen[k] --> SNPs seen in nanopore
    for(k in seen_illum){
        if(k in seen) continue
        # loop over illumina positions
        # ensures only coordinates missing in the nanopore data are handled
        
        split(k,p,":"); chrom=p[1]; pos=p[2]
        # extract CHROM and POS
        
        for(s in all_samples){
        # process each sample, even if genotype is missing
        
            illum_val=illum[k,s]
            if(illum_val==""||illum_val==".") illum_val="./.:NA"
            split(illum_val,a,":")
            gtI=a[1]; dpI=a[2]
            # retrieve illumina genotype
            
            gtN="./."; dpN="NA"
            # nanopore is always missing here
            
            hetI=(gtI=="0/1"||gtI=="1/0")?"yes":"no"
            # heterozygosity for illumina
            
            hetN="no"
            # heterozygosity for nanopore
            
            if(gtI=="./.") type="missing_both"
            else type="missing_nanopore"
            # if illumina is missing as well, classify as "missing_both"
            # otherwise nanopore is missing, but illumina is not
            
            count[s,type]++
            print chrom,pos,s,type,dpI,dpN,hetI,hetN >> out_detail
        }
    }

    # ---- Write summary ----
    print "SAMPLE","match","partial_match","no_match","missing_illumina","missing_nanopore","missing_both" > out_summary
    for(s in all_samples){
        print s,
              (count[s,"match"]+0),
              (count[s,"partial_match"]+0),
              (count[s,"no_match"]+0),
              (count[s,"missing_illumina"]+0),
              (count[s,"missing_nanopore"]+0),
              (count[s,"missing_both"]+0) >> out_summary
    # the +0 forces missing entries to become 0
    }
}
