## process the information that is filled in, creating RDS files in return
createFiles <- function(dfDir, metaDir, con1, con2, uniqueID) {
  ## Load data
  df <- read.csv(dfDir$datapath, header = TRUE, sep = ",", row.names = 1)
  meta <- read.csv(metaDir$datapath, header = TRUE, sep = ",")
  ## Match the meta data with the raw data
  meta <- meta[match(colnames(df), meta[,1]),]
  meta <- meta[complete.cases(meta),]
  ## Build the DESeq data set
  dds <- suppressWarnings(suppressMessages(DESeqDataSetFromMatrix(countData = df,
                                                                  colData = meta,
                                                                  design = ~ Description)))
  
  ## Filter lowly expressed genes (genes with a rowsum > 1 are kept)
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  ## Perform the differential expression analysis
  dds <- suppressMessages(DESeq(dds))
  
  ## Generate a genesymbol vector
  if (startsWith(rownames(dds)[1], "ENSG")) {
    geneSymbol <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=rownames(dds),column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first"))
  } else if (startsWith(rownames(dds)[1],"ENSMUSG")) {
    geneSymbol <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys=rownames(dds),column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first"))
  } else {
    geneSymbol <- rownames(dds)
  }
  
  ## Add geneSymbol
  geneSymbol <- geneSymbol[match(rownames(assay(dds)), geneSymbol$ENSEMBL),]
  ## Generate a new deseq2 dataframe with only the transcripts that have a gene symbol
  ddsfil <- (dds)[which(!is.na(geneSymbol$SYMBOL)),]
  ## Add the CPM normalized values
  assay(ddsfil, "cpm") <- cpm(assay(ddsfil))
  ## Add the gene symbols as a rownames in the ddsfil 
  rownames(ddsfil) <- geneSymbol$SYMBOL[which(!is.na(geneSymbol$SYMBOL))]
  ## Obtain the results
  res <- results(ddsfil, contrast=c("Description", con1, con2))
  ## Save the outcome to DEA.rds
  saveRDS(as.data.frame(res), paste0(uniqueID, "DEA.rds"))
  
  ## Normalize the data
  rld <- suppressMessages(vst(dds, blind = FALSE))
  ## Calculate the principal components
  b.res <- prcomp(t(assay(rld)))
  ## Save the prcomp outcome as PCA input
  saveRDS(b.res, paste0(uniqueID, "PCA-input.rds"))
  ## Save the quantitative information
  saveRDS(ddsfil, paste0(uniqueID, "QEA.rds"))
}
