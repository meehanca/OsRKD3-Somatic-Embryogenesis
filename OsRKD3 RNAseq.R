
setwd("/Users/u1762824/Documents/Rice/Black_rice/original_counts")

####################################################################
#Library
####################################################################
library(DESeq2)
library(ggplot2)
library(gplots)
library(reshape2)
library(pheatmap)
library(VennDiagram)
library(ggrepel)

sampleDataFilename <- 'sampleTable.txt'
sampleTable = read.table(sampleDataFilename,header=TRUE)
head(sampleTable)
htseqDir<-getwd()

##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = htseqDir,design = ~  condition)
## design<- you say to the test to do everything in relations to condition
## if you have more than one conditions you want to differentiate (for example different genotypes) you change design = ~  condition + genotype
##  And perform the analysis (details in the manual)

##  And perform the analysis (details in the manual)
dds<-DESeq(ddsHTSeq)

####################################################################
# Do PCA
####################################################################
#principal component analysis

vst = vst(dds)

v <- plotPCA(vst, intgroup=c("condition"))
v<- v+ geom_label_repel(aes(label = name))
v
pcaData <- DESeq2::plotPCA(vst, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA_parents.pdf", height = 6, width = 6)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  #scale_colour_manual(name="",values = c("a12"="goldenrod2", "gd33"="darkslateblue", "f1"="saddlebrown"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()
#dev.off()
####################################################################
#Plotting Reps
####################################################################
plot_reps =  function(dds,x=1,y=2,cond_choice=1, cond='condition'){
  ##  Estimate the size factors for normalisation
  dds<-estimateSizeFactors(dds)
  
  ## Extract the normalised counts for the condition you want
  rep_values<- counts(dds, normalized=TRUE)[,dds[[cond]]==cond_choice]
  
  # Take logs of these values
  vals <- log2(rep_values[,c(x,y)] + 0.5)
  # And plot
  plot(vals,pch=16, cex=0.4,xlab=paste('rep',x),ylab=paste('rep',y))
  grid(col = "darkgray", lty = "solid",lwd = par("lwd"), equilogs = TRUE)
  title(paste("Comparison of",cond_choice,"replicates"))
}

par(mfrow = c(1,3))

plot_reps(dds, x=1, y=2, cond_choice="TD")

plot_reps(dds, x=1, y=2, cond_choice="TM")

plot_reps(dds, x=1, y=2, cond_choice="ND")
####################################################################
#DEGs
####################################################################
filter_degs <- function(res){
  summary(res)
  res2 = res[!(is.na(res$padj)),]
  res2 = res2[res2$padj < 0.05,]
  return(res2)
}

resultsNames(dds)
TD_TM_DEGs = results(dds, contrast= c("condition", "TD", "TM"), alpha = 0.05, pAdjustMethod = "BH", lfcThreshold = 1)
TD_TM_DEG = filter_degs(TD_TM_DEGs)

TD_ND_DEGs = results(dds, contrast= c("condition", "TD", "ND"), alpha = 0.05, pAdjustMethod = "BH", lfcThreshold = 1)
TD_ND_DEG = filter_degs(TD_ND_DEGs)

ND_TM_DEGs = results(dds, contrast= c("condition", "ND", "TM"), alpha = 0.05, pAdjustMethod = "BH", lfcThreshold = 0)
ND_TM_DEG= filter_degs(ND_TM_DEGs)

summary(TD_TM_DEG)
head(TD_TM_DEG)

summary(TD_ND_DEG)
head(TD_ND_DEG)

summary(ND_TM_DEG)
head(ND_TM_DEG)

####################################################################
#Up and Downregulation
####################################################################

TD_TM_DEG_up <- TD_TM_DEG[TD_TM_DEG[,2]>0,]
TD_TM_DEG_down <- TD_TM_DEG[TD_TM_DEG[,2]<0,]

TD_ND_DEG_up <- TD_ND_DEG[TD_ND_DEG[,2]>0,]
TD_ND_DEG_down <- TD_ND_DEG[TD_ND_DEG[,2]<0,]

TM_ND_DEG_up <- ND_TM_DEG[ND_TM_DEG[,2]>0,]
TM_ND_DEG_down <- ND_TM_DEG[ND_TM_DEG[,2]<0,]

Shared_counts_up <- Shared_counts[Shared_counts[,2]>0,]
Shared_counts_down <- Shared_counts[Shared_counts[,2]<0,]
####################################################################
#MA Plots
####################################################################
par(mfrow=c(3,1))
DESeq2::plotMA(TD_TM_DEGs, ylim=c(-10,15), main='TD_TM_DEGs')
DESeq2::plotMA(TD_ND_DEGs, ylim=c(-10,15), main='TD_ND_DEGs')
DESeq2::plotMA(ND_TM_DEGs, ylim=c(-10,15), main='NDvsTM_DEGs')
####################################################################
#Volcano
####################################################################
library(EnhancedVolcano)

EnhancedVolcano(Shared,
                lab = rownames(Shared),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-12, 15),
                ylim = c(0,60))
####################################################################
#Venn
####################################################################
vp = venn.diagram(list(TD_TM_DEG = row.names(TD_TM_DEG), TD_ND_DEG = row.names(TD_ND_DEG)),fill = c("red", "blue"),
                  alpha = c(0.5, 0.5), cex = 2,lty =2, 
                  filename = NULL)

grid.newpage()
grid.draw(vp)
###################################################################
#Preparing lists of diffrentially expressed genes for online tools
####################################################################

TD_TM_DEG_write <- rownames(TD_TM_DEG)
TD_ND_DEG_write <- rownames(TD_ND_DEG)

write.table(TD_TM_DEG, file = "TD_TM_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(TD_ND_DEG, file = "TD_ND_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(TM_ND_DEG_up, file='TM_ND_DEG_up.txt',quote=F,row.names = T,col.names = T )
write.table(TM_ND_DEG_down, file='TM_ND_DEG_down.txt',quote=F,row.names = T,col.names = T )

write.table(TM_ND_DEG_up, file='TM_ND_DEG_up.txt',quote=F,row.names = T,col.names = T )
write.table(TM_ND_DEG_down, file='TM_ND_DEG_down.txt',quote=F,row.names = T,col.names = T )

Shared_DEGs<-subset(TD_TM_DEG, rownames(TD_TM_DEG) %in% rownames(TD_ND_DEG))

#counts
write.table(TD_TM_DEG, file = "TD_TM_DEG_counts.txt", quote = F, row.names = T, col.names = T)
write.table(TD_ND_DEG, file = "TD_ND_DEG_counts.txt", quote = F, row.names = T, col.names = T)

#### FPKM ####

Shared_counts_BR5<-read.table('Shared_counts_BR5_sorted.txt')
Shared_counts_BR6<-read.table('Shared_counts_BR6_sorted.txt')
Shared_gene_length<-read.table('Shared_gene_length_sorted.txt')

Shared_counts_BR5[,2]<-Shared_counts_BR5[,2]/Shared_gene_length[,2]
Shared_counts_BR6[,2]<-Shared_counts_BR6[,2]/Shared_gene_length[,2]

write.table(Shared_counts_BR5,'Shared_counts_BR5.txt',quote=F,col.names=F,row.names=F, sep='\t')
write.table(Shared_counts_BR6,'Shared_counts_BR6.txt',quote=F,col.names=F,row.names=F,sep='\t')

BR5<-(counts(dds)[,5])/(29375748/1000000)
BR6<-(counts(dds)[,6])/(27644884/1000000)

kilobase<-read.delim('gene_length.txt',header=F)
kilobase[,2]<-kilobase[,2]/1000

write.table(kilobase,'kilobase.txt',quote=F,col.names=F,row.names=F)
write.table(BR6,'counts_BR6',quote=F,col.names=F,row.names=T)

######New 

conversion <- read.delim('conversion', header=F)
conversion<-conversion[,c(2,1)]

LOC_BR5 <- merge(conversion, Shared_counts_BR5, by=1, all=TRUE)  # merge by row names (by=0 or by="row.names")
LOC_BR6 <- merge(conversion, Shared_counts_BR6, by=1, all=TRUE)

LOC_BR5<-unique(LOC_BR5[,c(2,3)])
LOC_BR6<-unique(LOC_BR6[,c(2,3)])

write.table(LOC_BR5,'LOC_BR5',quote=F,col.names=F,row.names=F)
write.table(LOC_BR6,'LOC_BR6',quote=F,col.names=F,row.names=F)
####################################################################
#Heatmap
####################################################################

counts = counts(dds , normalized = TRUE)
counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]#it removes genes that are not express and have no variance
colnames(counts) <- c("TM1","TM2","ND1","ND2","TD1","TD2")
Shared_counts <- (counts[rownames(Shared_DEGs),])
pheatmap((log2(counts+1)), scale = "row",border_color=NA,show_rownames = F,
         main = 'Shared DEGs expression across samples')

counts = counts(dds , normalized = TRUE)
counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]
colnames(counts) <- c("TM1","TM2","ND1","ND2","TD1","TD2")
TD_ND_DEG_counts <- (counts[rownames(TD_ND_DEG),])
pheatmap((log2(TD_ND_DEG_counts+1)), scale = "row",border_color=NA,show_rownames = F, main='TD vs ND DEGs expression across samples')


counts = counts(dds , normalized = TRUE)
counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]
colnames(counts) <- c("TM1","TM2","ND1","ND2","TD1","TD2")
TD_TM_DEG_counts <- (counts[rownames(TD_TM_DEG),])
pheatmap((log2(TD_TM_DEG_counts+1)), scale = "row",border_color=NA,show_rownames = F, main='TD vs TM DEGs expression across samples')

####################################################################
#Up and Downregulation
####################################################################

TD_TM_DEG_up <- TD_TM_DEG[TD_TM_DEG[,2]>0,]
TD_TM_DEG_down <- TD_TM_DEG[TD_TM_DEG[,2]<0,]

TD_ND_DEG_up <- TD_ND_DEG[TD_ND_DEG[,2]>0,]
TD_ND_DEG_down <- TD_ND_DEG[TD_ND_DEG[,2]<0,]


Shared<- TD_TM_DEG[rownames(TD_TM_DEG) %in% rownames(TD_ND_DEG),]
Shared_up <- TD_TM_DEG_up[rownames(TD_TM_DEG_up) %in% rownames(TD_ND_DEG_up),]
Shared_down <- TD_TM_DEG_down[rownames(TD_TM_DEG_down) %in% rownames(TD_ND_DEG_down),]

write(rownames(Shared_up), "Shared_up.txt")
write(rownames(Shared_down), "Shared_down.txt")
write(rownames(Shared), 'Shared.txt')
####################################################################
#Heatmap UP/DOWN 
####################################################################

#TDvsND

counts = counts(dds , normalized = TRUE)
counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]#it removes genes that are not express and have no variance
colnames(counts) <- c("TM1","TM2","ND1","ND2","TD1","TD2")

TD_ND_DEG_up <- (counts[rownames(TD_ND_DEG_up),])
pheatmap((log2(TD_ND_DEG_up+1)), scale = "row",border_color=NA,show_rownames = F,
         main = 'TDvsND upregulated DEGs')

TD_ND_DEG_down <- (counts[rownames(TD_ND_DEG_down),])
pheatmap((log2(counts(dds)+1)), scale = "row",border_color=NA,show_rownames = F,
         main = 'TDvsND downregulated DEGs')

#TDvsTM

TD_TM_DEG_up <- (counts)
pheatmap((log2(TD_TM_DEG_up+1)), scale = "row",border_color=NA,show_rownames = F,
         main = 'Global expression across samples',cluster_rows = TRUE)
[rownames(TD_TM_DEG_up)
TD_TM_DEG_down <- (counts[rownames(TD_TM_DEG_down),])
pheatmap((log2(TD_TM_DEG_down+1)), scale = "row",border_color=NA,show_rownames = F,
         main = 'TDvsTM downregulated DEGs')

#Shared_DEGs

Shared_up <- (counts[rownames(Shared_up),])
pheatmap((log2(Shared_up+1)), scale = "row",border_color=NA,show_rownames = F,
         main = 'Shared upregulated DEGs')

Shared_down <- (counts[rownames(Shared_down),])
pheatmap((log2(Shared_down+1)), scale = "row",border_color=NA,show_rownames = F,
         main = 'Shared downregulated DEGs')

Shared <- (counts[rownames(Shared),])
pheatmap((log2(Shared+1)), scale = "row",border_color=NA,show_rownames = F,
         main = 'Shared DEGs')

#Shared_tissue_expression_up 

Shared_up_converted <- read.delim('Shared_up_finalist.txt')
Shared_down_converted <- read.delim('Shared_down_finalist.txt')

####################################################################
#Heatmap Tissue Expression
####################################################################
setwd("~/Documents/Rice/Black_rice RNAseq/original_counts/expression_matrix_construction")
rice_expression_matrix <- read.delim('final.txt',header=F)
rownames(rice_expression_matrix) <- rice_expression_matrix[,1]
columns <- c(2:17,66:72,81:82)
rice_expression_matrix[,81:82]<-rice_expression_matrix[,81:82]*1000
rice_expression_matrix <- rice_expression_matrix[,columns] 
colnames(rice_expression_matrix) <- c("Leaves 20 days","Post-emergence inflorescence","Pre-emergence inflorescence",
                                      "Anther","Pistil","Seed-5 DAP","Embryo 25 DAP","Endosperm 25 DAP","Seed 10 DAP","Endosperm 25 DAP","Leaves 20 days","Shoots",rep("Seedling",4),
                                      "Callus","Leaf","Panicle before flowering","Panicle after flowering","Roots","Seed","Shoot","TD1","TD2")
rice_expression_matrix_filtered <- rice_expression_matrix[apply(rice_expression_matrix, MARGIN = 1, FUN = function(x) sd(x) != 0),]
MSU_rice <- pheatmap((log2(rice_expression_matrix_filtered+1)), main='Expression of shared DEGs across rice tissues',scale="row",cluster_rows = T,show_rownames = F)


#we are taking log2 to reduce the difference, srink the values 
# you can remove the scale if you want to remove that scale you will see the massive difference 
# +1 is important as well, but I don't know why :-(
#cluster_rows = F
#cluster_cols = F
# these arguments are by default set as T, if you change them to false interesting things can happen 

######## Dendograms #########

plot_tree <- function(pheatmap_Data){
  plot(pheatmap_Data$tree_row)
  max <- round(range(pheatmap_Data$tree_row$height))
  #change seq() to make different intersections on the tree
  for(i in seq(0.5,max[2],1)){
    abline(a = i, b = 0, col = "blue", lty = i+1 )
  }
}
produce_clusters <- function(pheatmap_Data, value){
  #grab gene names + cluster information from heatmap, melt + transform into dataframe
  
  melted.pheatmap_Data <- melt(cutree(pheatmap_Data$tree_row, h = value))
  melted.pheatmap_Data <- as.data.frame(melted.pheatmap_Data)
  #get counts for gene names + scale counts
  
  interesting_counts <- rice_expression_matrix_filtered[rownames(melted.pheatmap_Data),]
  interesting_counts = t(scale(t(as.matrix(interesting_counts)),center = TRUE, scale = TRUE))
  interesting_counts <- as.data.frame(interesting_counts)
  #move columns of count table into same order produced by pheatmap clustering
  interesting_counts <- interesting_counts[,pheatmap_Data$tree_col$order]
  #assign cluster ids to count table
  interesting_counts$cluster <- melted.pheatmap_Data$value
  interesting_counts$ids <- rownames(interesting_counts)
  melt_counts <<- melt(interesting_counts, id.vars = c("cluster", "ids"))
  #plot that shit
  c <- ggplot(melt_counts,aes(alpha = 0.5))
  
  c <- c + geom_line(aes(x = variable, y = value,
                         group = ids)) + theme_bw() + facet_grid(cluster~., scales = "fixed") +
    scale_color_manual(values = c("lightblue","goldenrod1", "orange", "red")) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = -90,
                                     hjust = 0),
          axis.title = element_blank(),
          legend.position = "top",
          axis.ticks.y = element_blank()) 
  c
}
cluster_to_heatmap <- function( cluster, rownames = F, scale = "row"){
  specific_genes <- unique((melt_counts[melt_counts$cluster==cluster,])[,2])
  pheatmap(log2((rice_expression_matrix_filtered[specific_genes,])+1), scale = "row", main = paste("Cluster", cluster), show_rownames = T, border_color = NA)
}
setwd("~/Documents/Rice/Black_rice RNAseq/original_counts")
###################################################################
#go
####################################################################
plot_go = function(goterms,name){
  goterms$percquery = goterms$queryitem/goterms$querytotal*100
  goterms$percback = goterms$bgitem/goterms$bgtotal*100
  filtered_go = goterms[goterms$FDR < 0.05,]
  #filtered_go = filtered_go[filtered_go$term_type == "P",]
  filtered_go_perc = cbind(filtered_go$percquery, filtered_go$percback)
  colnames(filtered_go_perc) = c("query","background")
  row.names(filtered_go_perc) = paste(filtered_go$Term,filtered_go$GO_acc,sep ="-->")
  meled = melt(filtered_go_perc)
  x = ggplot(meled, aes(Var1, value, fill=Var2)) +
    geom_bar(stat="identity", position="dodge")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    xlab("sig GO Term") +
    ylab("Ratio genes with term in list") +
    ggtitle(name)+
    coord_flip()
  plot(x)
  return(x)
}
#use website agriGo v2
#genes shared by TDvsTM and TDvsND
#select Oryza sativa japonica species
#run the list of Degs we have generated, make the result. download the details information (save as .txt file in the working directory or where you have saved all the data )

Induction_GO <- read.table('./Final/Go_terms_filtered_0.05_pvalue.txt', header = T,sep="\t")
#select the file you have just saved and generate plot
plot_go(Induction_GO,"0.0001")


require(randomForest)

# there are some random features to randomForest
# we can use a 'random seed' to make the random features reproducible
# you may use any seed you like, but for this example let's use the same one
set.seed(1234)

# lets try making a model
# we will predict sample type (x)
# using gene expression (y)
# and a small number of decision trees (ntree)

sampleTable <- read.table(file='sampleTable2.txt', header=T)
##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = htseqDir,design = ~  condition)

## design<- you say to the test to do everything in relations to condition
## if you have more than one conditions you want to differentiate (for example different genotypes) you change design = ~  condition + genotype
##  And perform the analysis (details in the manual)

dds<-DESeq(ddsHTSeq)

# x
t(counts(dds))

# y
colData(dds)$condition

# let's calculate the relative importance of all genes involved in making the decision tree

set.seed(1232)

name2<- randomForest(x=t(counts[rownames(Shared_DEGs),]), y=colData(dds)$condition, ntree=500, importance=T)
print(name)

gene_list1<-varImpPlot(name2, n.var=25, main="Top 25 Influential Transcripts")
write.csv(gene_list1, "geneset2.csv")
# plot importances

varImpPlot(name, n.var=25, main="Top 25 Influential Transcripts")

# now to refine the number of genes considered in each tree
# the randomForest package comes with a built in function that will try to optimize
# by trying mtry values higher or lower than the starting one, and looking for improvement
# when it can't improve by calculating a new mtry tree, it stops
mtry <- tuneRF(x = t(counts(dds)), y = colData(dds)$condition, ntreeTry=50, mtryStart=100,
               stepFactor=5,improve=1, trace=TRUE, plot=TRUE) #increases/decreases by step factor and checks for improve %improvement
print(mtry)

#extract influence metrics from random forest generation
influence <- importance(rf.50)
#order by error and output
head(influence[order(influence[,3], decreasing=T),])
influence.order <- influence[order(influence[,3], decreasing=T),]
write.table(influence.order, 
            file="randomForest_importance.txt", sep="\t", quote=F)

###################################################################
#transcription factors
###################################################################

####Upregulated

TFs <- read.table('TFs_shared_up_rice_final.txt',sep='\t',header=F)
colnames(TFs) <- c("1","TF class") 
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 48
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(TFs, aes(factor(TFs$`TF class`))) + ggtitle("Upregulated DEGs") +
  geom_bar(stat="count", position = "dodge",size=2) + theme_minimal() + theme(
    axis.text.x=element_text(angle = -90, hjust = 0, size = 11),
    axis.ticks.x=element_blank()) + 
  scale_fill_manual(values =c(rep("grey",15))) + xlab(element_blank()) + ylab("Number of Transcription Factors")

####All TFs

TFs <- read.table('TFs_shared_rice_final.txt',sep='\t',header=F)
colnames(TFs) <- c("1","TF class") 
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 36
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(TFs, aes(factor(TFs$`TF class`))) + ggtitle("All DEGs") +
  geom_bar(stat="count", position = "dodge",size=2) + theme_minimal() + theme(
    axis.text.x=element_text(angle = -90, hjust = 0),
    axis.ticks.x=element_blank()) + 
  scale_fill_manual(values =c(rep("grey",15))) + xlab(element_blank()) + ylab("Number of Transcription Factors")

####Downregulated


TFs <- read.table('TFs_shared_down_rice_final.txt',sep='\t',header=F)
colnames(TFs) <- c("1","TF class") 
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 36
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(TFs, aes(factor(TFs$`TF class`))) + ggtitle("Downregulated DEGs") +
  geom_bar(stat="count", position = "dodge",size=2) + theme_tufte() + theme(
    axis.text.x=element_text(angle = -90, hjust = 0),
    axis.ticks.x=element_blank()) + 
  scale_fill_manual(values =c(rep("grey",15))) + xlab(element_blank()) + ylab("Number of Transcription Factors")


plotCounts(dds, gene="OsRKD3", intgroup="condition",main ="OsRKD3", transform=F) 

###################################################################
#TF terms
###################################################################
TFs<- read.delim("~/Downloads/Osj_TF_list.txt")
LOC_Os <- read.delim("~/Downloads/convert_result (18).txt")
Os_TF<-unique(merge(TFs,LOC_Os,by.x=2,by.y=1)[,c(3,4)])
###################################################################
#Gage and list construction
###################################################################
library(gage)
library(ggplot2)
#Exclude lowly expressed genes for GSEA
DESeq2_negative_gene_IDs <- is.na(as.data.frame(TD_ND_DEGs$log2FoldChange))

###################################################################
list <- list()
for(i in 1:56){
  
  TF_class <- as.character(unique(Os_TF$Family))
  TF_class_name <- TF_class[i]
  
  list[[i]] <- Os_TF[grep(paste(TF_class_name),Os_TF$Family),2]
}
names(list)<-TF_class[1:56]



#Run GAGE command for all leaky and induced expressed transgenics

Enriched <- gage(counts(dds)[!DESeq2_negative_gene_IDs,],list,ref=c(3:4),samp=c(5:6),
                 rank.test = T,me.dir = F,
                 set.size=c(1,800), compare="unpaired")

Enriched_greater <- Enriched$greater[1:18,1:5]
Enriched_lesser <- Enriched$less[1:14,1:5]

Enriched_write <- rbind(Enriched_greater,Enriched_lesser)
write.table(Enriched_write,file = '../../analysis/TF_gene_enrichment.txt', quote =F, sep= "\t")
write.table(Enriched_greater,file = '../../analysis/TF_greater_gene_enrichment.txt', quote =F, sep= "\t")

q.val <- -log10(Enriched_write[,4])
q.val[16:32] <- q.val[16:32]*-1

data<-data.frame(rownames(Enriched_write), q.val)
colnames(data) <- c("TF","q.val")

colours <- c(rep("indianred1",15), rep("royalblue",17))

library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

library(wesanderson)

ggplot(data, aes(TF, q.val, fill=q.val)) + geom_bar(stat="identity") +
  scale_fill_continuous(low="blue", high="red") +
  coord_flip()

ggplot(userData, aes(month, count, fill = count)) +
  geom_bar(stat = "identity") +
  scale_x_date() + 
  scale_fill_continuous(low="blue", high="red") +
  labs(x= "Time", y="Count")


