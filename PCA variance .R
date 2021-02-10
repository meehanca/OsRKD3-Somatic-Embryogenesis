plotPCAEx = function(object, PCx = 1, PCy = 2, cond="condition", ntop=500, labels = TRUE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(cond %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  cond.df <- as.data.frame(colData(object)[, cond, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(cond) > 1) {
    factor(apply( cond.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[cond]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PCa=pca$x[,PCx], PCb=pca$x[,PCy], group=group, cond.df, name=colnames(object))
  
  pc1 <- ggplot(data=d, aes_string(x="PCa", y="PCb", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC",PCx,": ",round(percentVar[PCx] * 100),"% variance")) +
    ylab(paste0("PC",PCy,": ",round(percentVar[PCy] * 100),"% variance")) +
    coord_fixed()
  
  pc2 <- pc1 + geom_point()
  pc3 <- pc2  +     theme(legend.position="none")
  
  #  Finally add the labels, using ggrepel so that they dont write over each other or the points  
  if (labels)
  {
    library("ggrepel")
    pc3 + geom_label_repel(aes(label = name),
                           color = "gray20",
                           data = d,
                           force = 10)
  }
}

vst <- vst(dds)

pdf("./Figures/PCA.pdf", height = 5, width = 9,onefile=FALSE)
o <- plotPCAEx(vst, cond ="batch",1,2)
o +      theme(legend.position="right")

#######################################################################
#Top 100 most variable genes IDs
#######################################################################

rv <- rowVars(assay(dds))
rownames(counts(dds)[order(rv, decreasing=TRUE)[seq_len(min(100, length(rv)))],])
