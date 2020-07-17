
# Generation of the custom PCA, based on the number of dimentions that are found.
createCustomPCA <- function(B.resDat, targetFile, intgroup, targetIdentifier) {
  ggFig <- ""
  # Determine the number of interesting dim (based on a scree plot)
  dimNum <- checkNumberOfDim(B.resDat)
  if (dimNum == 1) {
    dimNum <- 2
  }
  
  # The information of the PC is altered so it can be used in the PCA plot. 
  pc.s = summary(B.resDat)$importance[2,1:dimNum]
  for (number in 1:dimNum) {
    pc.s = summary(B.resDat)$importance[2,1:dimNum]
    assign(paste0("pc", number, ".var"), round(pc.s[[paste0("PC", as.character(number))]],2))
  }
  # The rotation values are obtained.
  # rotatedData <- as.data.frame(B.resDat$rotation)
  
  # The counter is used to obtain the targetFile information of interest and to proceed when necessary.
  counter <- 1
  # The subtargets are the previously defined colNumbers, highlighting the interesting information for the PCA plot.
  newTarget <- as.data.frame(targetFile[,intgroup, drop=F])
  for (item in 1:ncol(targetFile[intgroup])) {
    usedTarget <- newTarget[,item]
    # The plot is drawn.
    rotatedData <- data.frame(PC1 = B.resDat$x[,1], PC2 = B.resDat$x[,2], name = rownames(targetFile))
    ggFig <- ggplot(rotatedData, aes(x=PC1, y=PC2, color=factor(usedTarget))) +
      # The color is dependant on the provided 'item'.
      geom_point(size=5, alpha = 0.5) +
      # The percentage of variance are added to the x- and y-axis.
      xlab(paste0("PC1: ",pc1.var*100,"% variance")) +
      ylab(paste0("PC2: ",pc2.var*100,"% variance")) +
      # The ID's are added to the different dots, geom_text_repel makes sure that the names are not overlapping.
      geom_text_repel(aes(label=as.character(targetIdentifier)), color="black", fontface="bold",
                      box.padding = unit(0.3, "lines"),
                      point.padding = unit(0.5, "lines"),
                      segment.colour = "grey50", alpha = 0.5) +
      labs(colour=colnames(newTarget)[counter]) + theme_bw() + 
      # Make sure that the background is clear.
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    counter <- counter + 1
  }
  return(ggFig)
}

###
# Checks the number of interesting dimensions
checkNumberOfDim <- function(B.resDat) {
  # Code for scree plot to identify the right number of PCA's
  par(mfrow=c(1,1))
  # biplot(B.resDat, scale = 0)
  std_dev <- B.resDat$sdev
  pr_var <- std_dev^2
  prop_varex <- pr_var/sum(pr_var)
  # Plot the scree plot.
  # plot(prop_varex, xlab = "Principal Component",
  #      ylab = "Proportion of Variance Explained",
  #      type = "b")
  dimNum <- countVals(prop_varex)
  return(dimNum)
}


###
# Count the number of PC's that is needed to explain at least 60% of the variance.
countVals <- function(vec, perc) {
  # The value of the previous counted value
  lastVal <- 0
  numberOfDim <- 0
  
  # The vec object consists of a df that contains the different PC's and the percentage that they explain.  
  for (i in 1:length(vec)) {
    # The lastVal is the current value + the value that was saved in lastVal (counting the different PC's based on a stepwise approach).
    lastVal <- vec[i] + lastVal
    # If the lastVal exceeds the 60% the numberOfDim is set to the number of PC's that are found.
    if (lastVal*100 >= 60 & numberOfDim == 0) {
      numberOfDim <- i
    }
  }
  # Print and return the number of dimensions found.
  # print(numberOfDim)
  return(numberOfDim)
}