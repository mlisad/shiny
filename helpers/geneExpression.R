## Create the gene expression plots
createCustomGeneExpressionPlot <- function(df, gene, con1, con2) {
  ## Determine the position of the filled in gene
  geneInt <- which(rownames(df) == gene)
  ## Obtain the values and transcform it to a dataframe that can be used by ggplot
  ggdf <- data.frame(sample=colnames(df), geneCount=assay(df)[geneInt,], condition=colData(df)$Description)
  
  ## Draw the ggplot
  ggplot(ggdf, aes(x = condition, y=geneCount, color=condition, fill=condition)) +
    ## Include the datapoints
    geom_point() + 
    ## Include a boxplot to determine the sample distribution
    geom_boxplot(alpha=0.8, width=0.2, position = position_dodge(0.9)) +
    ggplot2::theme_bw() +
    labs(title = gene, x= "", y="Expression(cpm)") +
    # Make sure that the background is clear.
    theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text = element_text(size=12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent"))
}