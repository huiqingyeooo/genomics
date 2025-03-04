# This R script was designed to loop through a list of taxids and retrieve the best sequences of each gene cluster from Genbank using the R package phylotaR.
# Huiqing Yeo

# Also refer to some helpful tutorials:
#tutorial:https://docs.ropensci.org/phylotaR/articles/phylotaR.html
#https://docs.ropensci.org/phylotaR/reference/drop_by_rank.html

# Note: I forked the phylotaR repository and made the following changes to the lines of code in R/stage4.R.
# This is to get phylotaR to form clusters even when there are less than 3 sequences present and report representative barcodes for all clusters.
#info(lvl = 1, ps = ps, "Dropping all clusters of < 1 sqs ...") #changed <3 sqs to <1 sqs
#all_clstrs <- all_clstrs[nsqs >= 1] # changed >=3 to >=1

setwd()
rm(list = ls())
library(phylotaR); library(ggplot2)

# Set paths
mainDir<-"/path/to/phylotar/work/folder"
ncbi<-"/path/to/ncbi/blast/program"
cluster.plots<-"/output/directory/for/cluster/plots"
fasta<-"/output/directory/for/fasta/sequences"
dat<-read.csv("txid_subset.csv"); txid.list<-dat$id # A list of genbank taxids corresponding to the taxa of interest

# Setup
for (taxa in txid.list){
    wd=paste0(dat[dat[,2]==taxa,1]) #wd=corresponding taxon name
    dir.create(file.path(mainDir, wd)) #create folder for each taxon
    txid=taxa
    setup(wd = wd, txid = txid, ncbi_dr = ncbi, v = TRUE, ncps=2, overwrite = TRUE)
    print(paste0("***** Setup done for:",wd))
}

# Run clustering pipeline
for (taxa in txid.list){
    wd=paste0(dat[dat[,2]==taxa,1])
    run(wd = wd)
    print(paste0("***** Run done for:",wd))
}


# Summarize results and plot to visualise taxa vs gene (cluster) completeness
df_total = data.frame()
for (taxa in txid.list){
  tryCatch({
    wd=paste0(dat[dat[,2]==taxa,1])
    all_clusters <- read_phylota(wd)
    smmry <- data.frame(summary(all_clusters)) #summarize phylota results
    smmry$taxa <- taxa  #add taxa id
    smmry$wd <- wd #add taxa name
    df <- data.frame(smmry)
    df_total <- rbind(df_total,df) #append it to dataframe
    print(paste0("***** Summary done for:",wd))
  }, error=function(e){cat("ERROR :",conditionMessage(e), wd, "\n")}) #added tryCatch function to skip taxa with only one cluster or no sequences

     species_txids <- get_txids(all_clusters, txids = all_clusters@txids, rnk = 'species') # get species-level taxonomic names
     species_txids <- unique(species_txids) #obtain unique txids 
     species_txids <- species_txids[species_txids !=  ''] # dropping missing (??)
     species_nms <- get_tx_slot(all_clusters, species_txids, slt_nm = 'scnm') #get species names
     species_nms <- sort(species_nms, decreasing = TRUE) # sort alphabetically for plotting
     tryCatch({
     p <- plot_phylota_pa(phylota = all_clusters, cids = all_clusters@cids, txids = species_txids,
                         txnms = species_nms) # generate geom_object
    }, error=function(e){cat("ERROR :",conditionMessage(e),wd, "\n")}) #added tryCatch function to skip taxa without sequences
     ggsave(paste0(cluster.plots,wd,".jpg"), p, width=12, height=8, units="in")
     print(paste0("***** Plot done for:",wd))
}

# Save dataframe
write.csv(df_total,"smmry_all.csv")

# Obtain best fasta sequences for cluster 1
for (taxa in txid.list){
  wd=paste0(dat[dat[,2]==taxa,1])
  all_clusters <- read_phylota(wd)
  tryCatch({
  smmry <- data.frame(summary(all_clusters)) #summarize phylota results
  cid_keep <- smmry[1, 'ID'] #select cluster 1
  selected <- drop_clstrs(phylota = all_clusters, cid = cid_keep) #drop clusters
  reduced <- drop_by_rank(phylota = selected, rnk = 'species', n = 1, #choose best seq per species
                          choose_by = c('pambgs', 'nncltds'), #criteria: ambiguous bases, length of sequence. Excluded: age of sequence
                          greatest=c(FALSE,TRUE)) #select seq with fewest number of ambgs bases, and longest seq length
  txids <- get_txids(phylota = reduced, cid = cid_keep, rnk = 'species') # get txids at the species level for each sequence
  scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm') # look up species names for txids
  scientific_names <- gsub('\\.', '', scientific_names) #clean names
  scientific_names <- gsub('\\s+', '_', scientific_names)
  sids <- reduced@clstrs[[cid_keep]]@sids   # Add sequence IDs
  scientific_names_cluster<-paste(scientific_names,sep="_",sids)
  write_sqs(phylota = reduced, sid = sids, sq_nm = scientific_names_cluster, 
            outfile = paste0(fasta,wd,'_c1_bestSeqs.fasta')) #export fasta
  print(paste0("***** Exported best sequences for:",wd))
  }, error=function(e){cat("ERROR :",conditionMessage(e), wd, "\n")}) #added tryCatch function to skip taxa without sequences
}

# Note: errors when extracting clusters 1-3 (cluster id 0-2) can be safely ignored.
# e.g., ERROR : [NA] not in records albuginosus - refers to groups that only have one cluster
# e.g., ERROR : 'data' must be of a vector type, was 'NULL' cancraedes - refers to groups where the txid was specified, but no sequences were available in genbank
