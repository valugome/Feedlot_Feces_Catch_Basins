top_taxa_legend <- function(data, taxlevel = "Family",  n = 10) { #Taxlevel takes a string ("Family", "Genus", etc), and data[[level]] gets the corresponding column from data (data should be a melted df)
  # If you've "factored" the taxlevel column, this takes that column and turns it into a character vector (I did this so later on I'd get actual taxa names, not the factor numbers)
  taxlevel_column <- as.character(data[[taxlevel]])
  #Then, aggregate abundances by tax level (getting the mean RA across samples)
  taxlevel_abundance <- aggregate(Abundance ~ taxlevel_column, data = data, mean) #Here, the result is a df with two columns, first one are taxa names, second one is the average abundance
  # I'll rename the first column of taxlevel_abundance. Instead of it being "data[[taxlevel]]", this changes it to the "taxlevel" provided.
  colnames(taxlevel_abundance)[1] <- taxlevel
  taxlevel_abundance <- taxlevel_abundance[order(-taxlevel_abundance$Abundance), ]# This orders taxa by abundance
  top_taxa <- head(taxlevel_abundance[[taxlevel]], n)# This selects the top (n) taxa and returns only the names
  
  # Now, to handle Others (just making sure it is at the bottom of the list, but may not be necessary)
  
  # Use a regular expression to identify top taxa that contain "Others" in their name
  others_taxa <- grep("Others", top_taxa, value = TRUE)
  
  # If there are any "Others" in the top taxa, remove them from the top_taxa (they'll be added later on at the bottom of the list)
  if (length(others_taxa) > 0) {
    top_taxa <- top_taxa[!top_taxa %in% others_taxa]  # Remove "Others" from the list, so that on the next step I add them at the bottom of the list 
  }
  
  ##Look for "Others" again, but this time in the whole taxlevel list 
  others_taxa <- grep("Others", taxlevel_abundance[[taxlevel]], value = TRUE)
  
  # Ensure that "Others" is always included at the bottom of the list of top taxa
  top_taxa <- c(top_taxa, others_taxa)
  #The check for others_taxa ensures that if any taxa containing "Others" are present in the top_taxa list, they are moved to the bottom. However, if no "Others" taxa are found in the top taxa, the code still adds "Others" to the end of the top_taxa list.
  
  return(top_taxa)
}