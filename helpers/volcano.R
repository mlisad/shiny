## Create the content for the dataframe
createDEGdf <- function(results, logFC, FDR, con1, con2) {
  ## Define Conditions as an empty column
  results$Condition <- ""
  ## Fill in the corresponding condition
  results$Condition[which(results$log2FoldChange > logFC & results$padj < FDR)] <- con1
  results$Condition[which(results$log2FoldChange < -logFC & results$padj < FDR)] <- con2
  ## Return a dataframe with the logFC, pval and FDR values
  as.data.frame(results)[which(results$Condition != ""),c(2,5,6)]
}

## Create the content for the volcano image
createCustomVolcano <- function(results, logFC, FDR, con1, con2) {
  ## Define Conditions as an empty column
  results$Condition <- ""
  ## Fill in the corresponding condition
  results$Condition[which(results$log2FoldChange > logFC & results$padj < FDR)] <- con1
  results$Condition[which(results$log2FoldChange < -logFC & results$padj < FDR)] <- con2
  ## Create the appropriate factor ordering
  results$Condition <- factor(results$Condition, levels=c(con1, "", con2))

  ## Draw the volcano plot
  ggplot(as.data.frame(results), aes(x = log2FoldChange, y=-log10(padj), color=Condition)) +
    ## If the condition is empty, then the alpha value is 0.2
    geom_point(aes(alpha=ifelse(Condition=="", 0.2,1)), show.legend=FALSE) +
    ggplot2::theme_bw() +
    scale_color_brewer(name= "Condition", palette="BrBG") +
    labs(title = paste0(con2, " | ", con1)) +
    # Make sure that the background is clear.
    theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.position = "none",
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent")) #+ # get rid of legend panel bg
}

