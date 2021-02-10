library(goseq)
library(tidyr)
library(dplyr)
###################################################################
#transcription factors
###################################################################

RAP_MSU_IDs<-read.delim('~/Downloads/RAP-MSU_2019-12-17.txt',header=F)
MSU_TF_IDs<-read.delim('~/Downloads/Osj_TF_list.txt')

TFs <- merge(RAP_MSU_IDs, MSU_TF_IDs, by.x=2, by.y=1)
TFs <- TFs[,c(2,4)]
TFs <- TFs[!grepl("None", TFs$V1),]

NRG5<-read.table('~/Downloads/Os_NRG5.txt')
NRG5[,2]<-c("NRG5 regulated")

Os_cluster3<-read.table('~/Downloads/Os_cluster3.txt')
Os_cluster3[,2]<-c("Os PBM cluster")

names(NRG5)<-c("V1","Family")
names(Os_cluster3)<-c("V1","Family")

TFs <- rbind(TFs,NRG5,Os_cluster3)

TFs_<- TFs[TFs$V1 %in% factor(rownames(Shared)),]

DESeq2_negative_gene_IDs <- is.na(as.data.frame(TD_ND_DEGs$log2FoldChange))

TFs_background<- TFs[TFs$V1 %in% DESeq2_positive_gene_IDs,]

###################################################################
#Gage and list construction
###################################################################

library(gage)

###################################################################
list <- list()
for(i in 1:58){
  
  TF_class <- as.character(unique(TFs$Family))
  TF_class_name <- TF_class[i]
  
  list[[i]] <- TFs[grep(paste(TF_class_name), TFs$Family),1]
}
names(list)<-TF_class[1:58]

###################################################################

gage(counts(dds)[!DESeq2_negative_gene_IDs,],list,ref=c(3,4),samp=c(1,2,5,6),
     rank.test = T,
     set.size=c(1,800), compare="unpaired",
     same.dir = F
     )

