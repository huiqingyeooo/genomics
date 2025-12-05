#!/usr/bin/awk -f
BEGIN {
    FS = "[ \t]+"; OFS = "\t"
    out_detail = "genotype_comparison.txt"
    out_summary = "genotype_summary.txt"
    print "CHROM", "POS", "SAMPLE", "TYPE", "DEPTH_illumina", "DEPTH_nanopore", "HET_illumina", "HET_nanopore" > out_detail
}

function trim(x){sub(/^[ \t]+/,"",x);sub(/[ \t]+$/,"",x);return x}

# ---- First file: illumina ----
FNR==1 && NR==FNR {
    if ($1=="CHROM"){for(i=3;i<=NF;i++){sample[$i]=i; all_samples[$i]=1} next}
}
NR==FNR {
    key=$1":"$2
    for(s in sample){illum[key,s]=$(sample[s])}
    seen_illum[key]=1
    next
}

# ---- Second file: nanopore ----
FNR==1 {
    for(i=3;i<=NF;i++){sample2[$i]=i; all_samples[$i]=1}
    next
}

{
    key=$1":"$2
    chrom=$1; pos=$2
    for(s in all_samples){
        illum_val=illum[key,s]; nano_val=$(sample2[s])
        if(illum_val==""||illum_val==".") illum_val="./.:NA"
        if(nano_val==""||nano_val==".") nano_val="./.:NA"
        split(illum_val,a,":"); split(nano_val,b,":")
        gtI=a[1]; dpI=a[2]; gtN=b[1]; dpN=b[2]
        hetI=(gtI=="0/1"||gtI=="1/0")?"yes":"no"
        hetN=(gtN=="0/1"||gtN=="1/0")?"yes":"no"

        # classify
        if(gtI=="./." && gtN=="./.") type="missing_both"
        else if(gtI=="./.") type="missing_illumina"
        else if(gtN=="./.") type="missing_nanopore"
        else if((gtI==gtN) || (gtI=="0/1"&&gtN=="1/0") || (gtI=="1/0"&&gtN=="0/1")) type="match"
        else if(((gtI~/[01]\/[01]/)&&(gtN~/[01]\/[01]/)) && ((gtI=="0/1"&&gtN=="0/0")||(gtI=="0/1"&&gtN=="1/1")||(gtI=="1/0"&&gtN=="0/0")||(gtI=="1/0"&&gtN=="1/1")||(gtI=="1/1"&&gtN~/0\/1|1\/0/)||(gtI=="0/0"&&gtN~/0\/1|1\/0/))) type="partial_match"
        else type="no_match"

        count[s,type]++
        print chrom,pos,s,type,dpI,dpN,hetI,hetN >> out_detail
    }
    seen[key]=1
}

# ---- After reading nanopore, handle illumina-only SNPs ----
END {
    for(k in seen_illum){
        if(k in seen) continue
        split(k,p,":"); chrom=p[1]; pos=p[2]
        for(s in all_samples){
            illum_val=illum[k,s]
            if(illum_val==""||illum_val==".") illum_val="./.:NA"
            split(illum_val,a,":")
            gtI=a[1]; dpI=a[2]
            gtN="./."; dpN="NA"
            hetI=(gtI=="0/1"||gtI=="1/0")?"yes":"no"
            hetN="no"
            if(gtI=="./.") type="missing_both"
            else type="missing_nanopore"
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
    }
}
