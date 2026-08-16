merge_low_abundance_ra <- function(data, threshold=1){ ##data has to be a phyloseq object with relative abundances
  otu.table <- as.data.frame(otu_table(data))
  otu.list <- row.names(otu.table[rowMeans(otu.table) < threshold,])
  merged <- merge_taxa(data, otu.list, 1)
  for (i in 1:dim(phyloseq::tax_table(merged))[1]){
    if (is.na(phyloseq::tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- paste("Others", "<" , threshold,"% RA")
      phyloseq::tax_table(merged)[i, 1:ncol(phyloseq::tax_table(merged))] <- paste("Others", "<" , threshold,"% RA")}
  }
  return(merged)
}
