merge_low_abundance_ARG_ra <- function(data, variable, level, threshold = 1) {
  # Convert to relative abundance (%)
  #data <- transform_sample_counts(data, function(x) {x / sum(x)} * 100)
  
  #TOOK OUT THE CONVERT TO RELATIVE ABUNDANCE. MUST GIVE A OBJECT ALREADY IN RELATIVE ABUNDANCE
  #Melt to long format
  melted <- psmelt(data)
  #Get taxonomic groups (given by level) with mean abundance < threshold in every variable group
  low_abund <- melted %>%
    group_by(!!sym(variable), !!sym(level)) %>% #since "variable" is a string, sym turns it into a symbol and then !! unquotes it, so that group_by can identify it as a column name
    summarise(mean_abun = mean(Abundance), .groups = "drop") %>% ##For each group (sample_type × taxonomic group), calculate mean abundance 
    group_by(!!sym(level)) %>% #Going to check if taxonomic group (level) are below threshold 
    summarise(all_below = all(mean_abun < threshold)) %>% # Check if mean abundance for the taxonomic group is less than threshold in all of the variable groups
    #If a taxonomic group has mean_abun < threshold for ALL of the variable groups, it will get a value of TRUE in all_below. Otherwise it will get FALSE
    filter(all_below) %>% # keeps only those tax groups where any_below is TRUE
    pull(!!sym(level)) #extract tax group names

  #Merge all OTUs that belong to low-abundance taxonomic groups
  tax_table_df <- as.data.frame(phyloseq::tax_table(data))
  taxa_to_merge <- tax_table_df %>%
    filter(!!sym(level) %in% low_abund)%>%
    rownames() ##Pull OTUs that beling to those low-abudance tax-groups
  merged <- merge_taxa(data, taxa_to_merge, 1)
  
  # Rename the merged "Other" taxon
  for (i in 1:nrow(phyloseq::tax_table(merged))) {
    if (is.na(phyloseq::tax_table(merged)[i, 2])) {
      taxa_names(merged)[i] <- paste0("Others <", threshold, "% RA")
      phyloseq::tax_table(merged)[i, 1:ncol(phyloseq::tax_table(merged))] <- paste0("Others <", threshold, "% RA")
    }
  }

  return(merged)
}
