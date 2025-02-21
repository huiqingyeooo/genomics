# This R script is used to generate a time-calibrated ultrametric tree using a series of node ages.
# Huiqing Yeo

# Useful link: http://blog.phytools.org/2024/07/obtaining-time-calibrated-ultrametric.html
rm(list=ls())
setwd("")
library(phytools)
library(phangorn)

# read in tree
tree<-read.tree("20_S2017S2023genbank_monophyly.treefile")

# root tree
node<-which(tree$tip.label=="Culex_Culex_quinquefasciatus")
len<-tree$edge.length[which(tree$edge[,2]==node)] 
ml.tree<-bind.tip(tree=tree,tip.label="outgroup",where=node,edge.length=len,position=0.1)
ml.tree<-root(ml.tree,outgroup="Culex_Culex_quinquefasciatus")
ml.tree<-drop.tip(ml.tree,"outgroup")

# find nodes and add ages
nodes<-c(findMRCA(ml.tree,c("Culex_Culex_quinquefasciatus","Psorophora_Grabhamia_cingulata")), #1
         findMRCA(ml.tree,c("Psorophora_Grabhamia_cingulata","Psorophora_Psorophora_howardii")), #2
         findMRCA(ml.tree,c("Aedes_Stegomyia_aegypti","Aedes_Stegomyia_albopictus")), #3
         findMRCA(ml.tree,c("Verrallina_Verrallina_lineata","Aedes_Aedes_cinereus")), #4
         findMRCA(ml.tree,c("Aedes_Fredwardsius_vittatus","Opifex_Opifex_fuscus")), #5
         findMRCA(ml.tree,c("Aedes_Ochlerotatus_spencerii","Aedes_OchlerotatusChrysoconops_pallens")), #6
         findMRCA(ml.tree,c("Aedes_Fredwardsius_vittatus","Aedes_Aedimorphus_centropunctatus")), #7
         findMRCA(ml.tree,c("Aedes_OchlerotatusFinlayaSA_candidoscutellum","Aedes_Rampamyia_quinquelineatus")), #8
         findMRCA(ml.tree,c("Aedes_Kompia_purpureipes","Aedes_Lewnielsenius_muelleri")), #9
         findMRCA(ml.tree,c("Aedes_OchlerotatusCulicada_canadensis","Aedes_Ochlerotatus_fitchii")), #10
         findMRCA(ml.tree,c("Haemagogus_Conopostegus_leucocelaenus","Aedes_OchlerotatusCulicelsa_taeniorhynchus"))) #11
age.min=c(71.2449,28.9922,25.0253,25.1549,52.0234,
          26.0171,37.4407,25.0169,14.6118,22.8726,26.4)
age.max=c(108.402,44.8503,38.4289,39.3621,77.4113,
          38.3484,55.3936,37.0298,27.603,33.9044,39.1343)
age.mean=c(87.6072,35.8929,30.9185,31.5328,62.7468,
           31.4469,45.212,30.2981,20.7457,27.7271,31.9943)

pdf("S2017S2023genbank_subgen.agePoints.pdf", width=15, height=55)
plotTree(ladderize(ml.tree),fsize=0.8)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv) #get last plotted tree http://blog.phytools.org/2014/03/extracting-phylo-object-from-last.html
points(obj$xx[nodes],obj$yy[nodes],pch=21,bg=palette(rainbow(11)),cex=2)
legend("bottomleft",paste("(",age.max,", ", age.min,") mya",sep=""),pch=21,
       pt.bg=palette(rainbow(11)),pt.cex=2,bty="n")
dev.off()

# create calibration points
# with cx pip outgroup
calibration.relative<-makeChronosCalib(ml.tree,node=nodes,age.min=age.min,age.max=age.max)
calibration.relative

calibration.absolute<-makeChronosCalib(ml.tree,node=nodes,age.min=age.mean,age.max=age.mean)
calibration.absolute

# obtain time-calibrated ultrametric tree
# fit model with penalized likelihood
pl.tree.rel<-chronos(ml.tree,lambda=1,calibration=calibration.relative,model="relaxed",
                 control = chronos.control(iter.max=1e4,dual.iter.max=20))

pl.tree.abs<-chronos(ml.tree,lambda=1,calibration=calibration.absolute,model="relaxed",
                 control = chronos.control(iter.max=1e4,dual.iter.max=20))

### log-Lik (log likelihood) 
# measure of how well model fits the data
# higher log-likelihood indicates better fit to the data
### PHIIC (phylogenetic information criterion)
# model selection criteria that penalizes for model complexity to prevent overfitting
# goal is to balance fit and model simplicity, lower PHIIC values are better

# Compare with another lambda value (smoothing parameter)
pl.tree.rel.10<-chronos(ml.tree,lambda=10,calibration=calibration.relative,model="relaxed",
                        control = chronos.control(iter.max=1e4,dual.iter.max=20))

# Plot tree
pdf("S2017S2023genbank_subgen.chronos.pdf", width=15, height=55)
plotTree(pl.tree.rel,direction="leftwards",
         xlim=c(174,-40),ftype="i",mar=c(4.1,1.1,0.1,1.1),
         fsize=0.8)
axis(1)
title(xlab="millions of years before present")
abline(v=seq(0,150,by=50),lty="dotted",col="grey")
dev.off()

# Comparison of edge lengths of the two trees
plot(pl.tree.rel$edge.length,pl.tree.rel.10$edge.length,
     pch=21,bg="grey",cex=1.2,bty="n",
     xlab=expression(paste(lambda,"= 1")),
     ylab=expression(paste(lambda,"= 10")))
lines(c(0,80),c(0,80))
legend("topleft","1:1 line",lty="solid",bty="n")
grid()
title(main=expression(paste(
  "Comparison of edge lengths with two different ", lambda," values")))

# Export time calibrated tree
ape::write.tree(pl.tree.rel, file="S2017S2023genbank_subgen.chronos.treefile")

# Subset tips in time calibrated tree to RHK dataset
sp.list<-read.csv("./constraint_tree/list_RHK.csv", header=T, row.names = 1) #read in list of RHK species
rhk.ml<-geiger::treedata(pl.tree.rel, sp.list, sort=TRUE, warnings=TRUE)
tree<-rhk.ml$phy

library(ggtree)
p<-ggtree(rhk.ml$phy)+geom_tiplab()+geom_text2(aes(subset=!isTip, label=label))
pdf("S2017S2023genbank_subgen.chronos.RHK.pdf", width=45, height=50)
p
dev.off()

# Rename tip labels to match scales and setae dataset
names.df<-read.csv("rename.csv", header=T)
old<-names.df$standardizedNames; new<-names.df$names
tree$tip.label[match(old, tree$tip.label)] <- new
ape::write.tree(rhk.ml$phy, file='S2017S2023genbank_subgen.chronos.RHK.treefile')
