# This R code deals with adding a list of tips to a phylogeny to create a backbone tree for constrained tree search in iqtree2.
# In this project, I am adding additional species with sequence data (from Soghigian et al 2023 and Genbank) to the Soghigian et al 2017 phylogeny.
# Huiqing Yeo

rm(list=ls())
setwd()
library(phangorn)
library(phytools)
library(treeio)
library(tidytree)
library(dplyr)
library(ape)
library(ggtree)

# Read in Soghigian et al 2017 tree
tree<-read.tree("12862_2017_1092_MOESM14_ESM.nwk")

# Collapse nodes with bootstrap <95
tree95<-as.polytomy(tree, feature='node.label', fun=function(x) as.numeric(x) < 95)

# Export tree to collapse genus/subgenus into polytomies in treeviewer
ape::write.tree(tree95, file="S2017.collapsed95.treefile")

# Read in polytomy treefile
tree<-read.tree("S2017.collapsed95.plytmy.tre")

# Rename tip labels in tree to standardized species names
names.df<-read.csv("list_S2017.csv", header=T) #csv file with two columns: a column with the tip labels, and a column to rename the tip labels
old<-names.df$name; new<-names.df$standardizedNames
tree$tip.label[match(old, tree$tip.label)] <- new
write.tree(tree) #check names have been changed

# Drop Aedes_Aedimorphus_dentatus (wrong placement in S2017 tree)
tree.1<-drop.tip(tree, c("Aedes_Aedimorphus_dentatus"), root.edge = 0)
tree<-tree.1; tips<-tree$tip.label

# Subset representative tips from each genus/subgenus/monophyletic group
sp.list<-read.csv("addTips.csv", header=T, row.names = 1)
dat<-sp.list%>% filter(position %in% c("keep"))
rownames(dat) <- dat$spp
dat.keep.ml<-geiger::treedata(tree, dat, sort=TRUE, warnings=TRUE)
tree<-dat.keep.ml$phy

### Add dropped tips from S2017, and tips from S2023 and genbank to S2017 tree
sp.list$position<-as.factor(sp.list$position)

### Add s2017 taxa by specifying tips
dat<-sp.list%>% filter(position %in% c("4"))
list.2017<-dat$spp

for (taxa in list.2017){
  print(paste0("***** Working on:",taxa))
  tip=paste0(sp.list[sp.list[,1]==taxa,2])
  node<-which(tree$tip.label==tip)
  len<-tree$edge.length[which(tree$edge[,2]==node)] 
  merged.tree<-bind.tip(tree=tree,tip.label=taxa,edge.length=len,where=node,position=0)
  tree<-merged.tree
  print(paste0("***** Added tip for:",taxa))
}

### Add taxa by specifying nodes
#https://www.mail-archive.com/search?l=r-sig-phylo@r-project.org&q=subject:%22Re%5C%3A+%5C%5BR%5C-sig%5C-phylo%5C%5D+Adding+random+branches+to+tree+recursively%22&o=newest&f=1

dat1<-sp.list%>% filter(position %in% c("0","1"))
list1<-dat1$spp

for (taxa in list1){
  print(paste0("***** Working on:",taxa))
  mrca.1=paste0(sp.list[sp.list[,1]==taxa,2])
  mrca.2=paste0(sp.list[sp.list[,1]==taxa,3])
  pos=as.integer(paste0(sp.list[sp.list[,1]==taxa,4]))
  node<-getMRCA(tree, c(mrca.1, mrca.2))
  tip<-which(tree$tip.label==mrca.1)
  len<-tree$edge.length[which(tree$edge[,2]==node)] #need to check why tree$edge[,2]
  merged.tree<-bind.tip(tree=tree,tip.label=taxa,where=node,edge.length=len,position=(pos*len*0.5))
  tree<-merged.tree
  print(paste0("***** Added tip for:",taxa))
}

### Add taxa by specifying tips
dat2<-sp.list%>% filter(position %in% c("2","3"))
list2<-dat2$spp

for (taxa in list2){
  print(paste0("***** Working on:",taxa))
  tip=paste0(sp.list[sp.list[,1]==taxa,2])
  node<-which(tree$tip.label==tip)
  len<-tree$edge.length[which(tree$edge[,2]==node)] 
  merged.tree<-bind.tip(tree=tree,tip.label=taxa,edge.length=len,where=node,position=0)
  tree<-merged.tree
  print(paste0("***** Added tip for:",taxa))
}

# Drop Opifex species and add them to root of Aedes clade A
#(Opifex,(All other clade A Aedes and other groups))
tree.2<-drop.tip(tree, c("Opifex_Nothoskusea_chathamicus","Opifex_Opifex_fuscus"), root.edge = 0)
tree<-tree.2

# Add Opifex_Nothoskusea_chathamicus
node<-getMRCA(tree,c("Aedes_Ochlerotatus_procax","Aedes_Ochlerotatus_nigrithorax"))
len<-tree$edge.length[which(tree$edge[,2]==node)]
merged.tree<-bind.tip(tree=tree,tip.label="Opifex_Nothoskusea_chathamicus",where=node,
                      edge.length=len,position=len*0.5)
tree<-merged.tree

# Add Opifex_Opifex_fuscus
node<-which(tree$tip.label=="Opifex_Nothoskusea_chathamicus")
len<-tree$edge.length[which(tree$edge[,2]==node)]
merged.tree<-bind.tip(tree=tree,tip.label="Opifex_Opifex_fuscus", edge.length=len,where=node,position=0)
tree<-merged.tree

# Add colour to plot to visualise which the added tips and check if they are in the correct place
library(ggtree)
dd<-as.data.frame(cbind(sp.list$spp,sp.list$source))
str(dd)
dd$V2<-as.factor(dd$V2)
colnames(dd)<-c("spp","source")
p<-ggtree(tree)+geom_tiplab()+geom_text2(aes(subset=!isTip, label=label))
p.col <- p %<+% dd + 
  xlim(0,0.75)+
  geom_tiplab(aes(fill = factor(source)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) # size of label border

pdf("S2017.S2023.genbank.monophyly.pdf", width=15, height=55)
p.col
dev.off()

# subset tips to taxa with sequences
sp.list<-read.csv("./list_S2017S2023genbank.csv", header=T, row.names = 1)
constrain.ml<-geiger::treedata(tree, sp.list, sort=TRUE, warnings=TRUE)
tree$tip.label
tree<-constrain.ml$phy #367 tips 116 internal nodes

# Export tree file
tree$node.label <- NULL
ape::write.tree(tree, file='S2017.S2023.genbank.monophyly.treefile')

# Move on to iqtree, then dating tree (dateTree.R)
