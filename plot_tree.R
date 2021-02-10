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
  
  interesting_counts <- counts[rownames(melted.pheatmap_Data),]
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
  specific_genes <<- unique((melt_counts[melt_counts$cluster==cluster,])[,2])
  pheatmap(log2((counts[specific_genes,])+1), scale = "row", main = paste("Cluster", cluster), show_rownames = rownames, border_color = NA)
}