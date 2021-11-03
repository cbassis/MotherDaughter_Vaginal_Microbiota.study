##Code for a heatmap
##Modeled after code written by Anna Seekatz: https://github.com/aseekatz/ERIN.recurrence/blob/master/Rcode_erinsubset/erinsubset_Fig3.heatmap.R
##September 2021 correction

library(RColorBrewer)
library(gplots)
library(vegan)
library(plyr)

setwd("~/Dropbox (University of Michigan)/from_box/Mother_daughter_analysis/121317_mother_daughter_copy/September_2021/")
shared<-read.table(file = "../motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.pick.shared", header = TRUE, row.names = 2)
dim(shared)
otu<-subset(shared, select =-c(label, numOtus))
otu.matrix<-as.matrix(otu)
otu.rel<-otu.matrix/rowSums(otu.matrix)
otu.rel.max<-apply(otu.rel,2,max)
otu.rel.filtered<-otu.rel[,otu.rel.max>0.02]
myCol3 <- colorRampPalette(brewer.pal(8,"Greys"))(8)
myBreaks2 <- c(0, 0.001, 0.003, 0.01, 0.05, 0.10, 0.50, 0.80, 1) 
par(mar=c(0,0,0,0))
heatmap.2(otu.rel.filtered, col=myCol3, breaks=myBreaks2, cexRow=0.5, cexCol=0.5, trace="none", dendrogram="none") 

seqok<-read.table(file="../seqok.sample.list.txt", header = TRUE)
meta<-read.table(file="../../motherdaughter.meta.heatmap.txt", header = TRUE) #rewrites meta table, downloaded this file from Kaylie on Box, edited pair designations (removed PAIR_, e.g. now says A rather than PAIR_A)
meta.seqok<-merge(seqok, meta, by.x=c("Group"), by.y=c("Sample"))
within_pair_dist<-read.table(file="../average.within.pair.distances.plus.txt", header=TRUE) #Made this file on 2/9/18, see 121517.distance.seqok.R
meta.dist<-merge(meta, within_pair_dist, by.x=c("Pair"), by.y=c("Pair"))
meta.otu.rel.filtered<-merge(meta.dist, otu.rel.filtered, by.x=c("Sample"), by.y=c("row.names"))
write.table(meta.otu.rel.filtered, "motherdaughter.otu.meta.heatmap.txt", quote=FALSE, sep = "\t", col.names = NA)

### made motherdaughter.otu.heatmap.ordered.txt in excel: starting with motherdaughter.otu.meta.heatmap.txt, copying Sample column to Time_sample and replacing baseline B with 0,
### sorting Distance (smallest to largest), then Pair (A to Z), then Subject1 (Z to A), then Time_sample ; deleted column A
### saved as motherdaughter.otu.heatmap.ordered.wcolumns.txt, deleted all columns except Sample and OTU rel abundances, saved as motherdaughter.otu.heatmap.ordered.txt

otu.ordered.rel.filtered<-read.table("motherdaughter.otu.heatmap.ordered.txt", header=TRUE, row.names=1) ##note: this is the updated file with corrected relative abundances
meta2<-read.table(file = "motherdaughter.otu.heatmap.ordered.wcolumns.txt", header = TRUE, row.names = 1) ##note: this is the updated file with corrected relative abundances
otu.ordered.rel.filtered.matrix<-as.matrix(otu.ordered.rel.filtered)

##To check out heatmap:
heatmap.2(otu.ordered.rel.filtered.matrix, col=myCol3, breaks=myBreaks2, cexRow=0.2, cexCol=0.5, trace="none", density.info="none", dendrogram="none", Rowv=FALSE)

#Made OTUlist.txt from OTUs in motherdaughter.otu.heatmap.ordered.txt, deleted " from motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy, saved as unix file
###Made motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy.subset.txt with Alyx Schubert's code:
## in terminal: perl ~/Box/Microbiome_Explorer_Program/Aslam_102615_no_females/Cecal_male_103015_Microbiome/Alyxcode.pl ~/Box/Mother_daughter_analysis/121317_mother_daughter_copy/September_2021/OTUlist.txt ~/Box/Mother_daughter_analysis/motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy 
#9/3/21 output in "~/Box Sync/Mother_daughter_analysis/", copied to working directory, "~/Box Sync/Mother_daughter_analysis/121317_mother_daughter_copy/September_2021":motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy.subset.txt
## In Excel added headers, numbers in Phylum column to get in order that I want

##Makes column/taxonomy labels
tax<-read.table(file="motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy.subset.txt", header=TRUE) ##note: this is the updated file although may not be different than original
tax<- tax[order(tax$OTU),]
tax.col<-tax[,c("OTU", "Phylum", "Size")]
tax.col$color<-mapvalues(tax.col$Phylum, from = c("01_Firmicutes","02_Firmicutes", "03_Actinobacteria", "04_Actinobacteria", "05_Bacteroidetes", "06_Fusobacteria", "07_Tenericutes"), to = c("hotpink", "red", "deepskyblue1","dodgerblue4", "darkgreen", "darkorange1", "tan"))
tax.col<- tax.col[order(tax.col$Phylum, -tax.col$Size) , ]
tax.col<-t(tax.col)
phyla.col<-tax.col[4,]
col.order<-as.character(tax.col[1, ])  	
otu.ordered.rel.filtered.matrix<-otu.ordered.rel.filtered.matrix[,col.order]
rbind(colnames(otu.ordered.rel.filtered.matrix), tax.col)
cbind(row.names(otu.ordered.rel.filtered.matrix), meta2)  
heatmap.2(otu.ordered.rel.filtered.matrix, col=myCol3, breaks=myBreaks2, cexRow=0.5, cexCol=0.5, trace="none", density.info="none", Colv=F, Rowv=F, dendrogram="none", ColSideColors = as.character(phyla.col))

### To get colors for key:
display.brewer.pal(8, "Greys")

##Makes row labels
#Daughter's and mother's birth mode
meta2$birth<-mapvalues(meta2$X15, from = c("1", "2"), to = c("yellow", "white"))
heatmap.2(otu.ordered.rel.filtered.matrix, col=myCol3, breaks=myBreaks2, cexRow=0.5, cexCol=0.5, trace="none", density.info="none", Colv=FALSE, Rowv=FALSE, labRow=FALSE, dendrogram="none", key=FALSE, RowSideColors = as.character(meta2$birth), ColSideColors = as.character(phyla.col))
heatmap.2(otu.ordered.rel.filtered.matrix, col=myCol3, breaks=myBreaks2, cexRow=0.3, cexCol=0.5, trace="none", density.info="none", Colv=FALSE, Rowv=FALSE, dendrogram="none", key=FALSE, RowSideColors = as.character(meta2$birth), ColSideColors = as.character(phyla.col))
#Reproductive stage
#based on X6 from motherdaughter.otu.heatmap.ordered.wcolumns.txt: 1-6=Reproductive, 7=Menopausal, na=Premenarche
meta3<-read.table(file = "../Reproductive.stage.txt", header = TRUE, row.names = 1)
meta3$stage<-mapvalues(meta3$X6, from = c("Premenarche", "Reproductive", "Postmenopausal"), to = c("white", "red", "yellow"))
heatmap.2(otu.ordered.rel.filtered.matrix, col=myCol3, breaks=myBreaks2, cexRow=0.5, cexCol=0.5, trace="none", density.info="none", Colv=FALSE, Rowv=FALSE, labRow=FALSE, dendrogram="none", key=FALSE, RowSideColors = as.character(meta3$stage), ColSideColors = as.character(phyla.col))
heatmap.2(otu.ordered.rel.filtered.matrix, col=myCol3, breaks=myBreaks2, cexRow=0.3, cexCol=0.5, trace="none", density.info="none", Colv=FALSE, Rowv=FALSE, dendrogram="none", key=FALSE, RowSideColors = as.character(meta3$stage), ColSideColors = as.character(phyla.col))

#Then used Adobe Illustrator CS6 to finalize labels and colors