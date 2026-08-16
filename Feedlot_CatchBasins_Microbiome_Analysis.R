#LOAD R PACKAGES ######
setwd('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Analysis_repository/Analyses/R Analyses')
# Installs (if needed) and loads packages from CRAN, Bioconductor, GitHub
load_packages <- function(cran_pkgs = character(0),
                          bioc_pkgs = character(0),
                          github_pkgs = character(0)) {
  
  # Make sure BiocManager and remotes are available for installs
  if (length(bioc_pkgs) > 0 && !requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (length(github_pkgs) > 0 && !requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  #  CRAN packages 
  for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing CRAN package: ", pkg)
      install.packages(pkg)
    }
  }
  
  # Bioconductor packages 
  for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing Bioconductor package: ", pkg)
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  }
  
  # GitHub packages (named vector: pkg_name = "owner/repo") 
  for (pkg in names(github_pkgs)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing GitHub package: ", github_pkgs[[pkg]])
      remotes::install_github(github_pkgs[[pkg]])
    }
  }
  
  # Finally, load all
  all_pkgs <- c(cran_pkgs, bioc_pkgs, names(github_pkgs))
  invisible(lapply(all_pkgs, function(pkg) {
    library(pkg, character.only = TRUE)
  }))
  
  message("All packages loaded")
}

##Packages lists
#CRAN
cran_pkgs = c(
  "tidyverse", "ggplot2", "stringr", "dplyr", "vegan", "cowplot",
  "ggdendro", "randomcoloR", "ggpubr", "ggsignif", "UpSetR",
  "ggtext", "ggnewscale", "rstatix", "ggrepel", "ggh4x", "svglite",
  "writexl", "paletteer", "lme4", "lmerTest", "car", "emmeans",
  "Polychrome", "colorspace", "devtools", "remotes"
)
#BioConductor
bioc_pkgs = c(
  "phyloseq", "metagenomeSeq", "ANCOMBC", "maaslin3",
  "MicrobiotaProcess", "microbiome"
)
#GitHub
github_pkgs = c(
  metagMisc      = "vmikk/metagMisc",
  pairwiseAdonis = "pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",
  maaslin3       = "biobakery/maaslin3"  
)

#Load them with the function 
load_packages(cran_pkgs, bioc_pkgs, github_pkgs)

#SOURCE FUNCTIONS#####
source('MergeLowAbundanceOthersPercentage.R')
source("top_taxa_legend_updated.R")
source('fill_taxonomy_updated.R')
source('MergeLowAbun_group_microbiome.R')


#DATA FROM KRAKEN OUTPUT#### 
##Classified counts matrix ###
counts <- read_csv('kraken_analytic_matrix.conf_0.0_merging.csv')

#Unclassified counts data frame###
unclassified_counts <- read_csv('unclassifieds_kraken_analytic_matrix.conf_0.0_merging.csv')

##Changing col names (sample IDs), to match them as they are in the metadata file
counts_newcolnames <- sapply(strsplit(colnames(counts), "_"), `[`, 1) #Splitting col names by "_", then extracting the first part of each split column name 
colnames(counts) <- counts_newcolnames #replacing col names for new ones
colnames(counts) ##good, now sample names are OK!

##Unclassified counts, change SampleID
unclassified_counts <- unclassified_counts %>%
  mutate(SampleID = sapply(strsplit(SampleID, "_"), `[`, 1) )
unclassified_counts$SampleID ##Good now

#Change zymo- to zymo. as they are in metadata
unclassified_counts[which(unclassified_counts$SampleID == "Zymo-1a"), "SampleID"] <- "Zymo.1a"
unclassified_counts[which(unclassified_counts$SampleID == "Zymo-1b"), "SampleID"] <- "Zymo.1a"
unclassified_counts$SampleID ##Good now

#Clean up
counts <- counts%>%
  rename(taxonomy=taxa)%>%##changing 'taxa' column name to taxonomy
  rename("Zymo.1a" = `Zymo-1a`,
         "Zymo.1b" = `Zymo-1b`)#Change zymo- to zymo. as they are in metadata
colnames(counts) ##okay, now the first column is "taxonomy" and the following ones are the sample IDs 

##Separating into taxonomy levels
counts_separated_tax <- counts %>%
  separate(taxonomy, 
           into = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "\\|",#splits the strings by the "|" symbol.
           fill = "right") #fill = "right", missing components are added as "NA" to the right (last columns) instead of to the left 

##Extracting just taxonomy  (columns 1:8 are taxonomy, the rest are counts)
tax.table<- counts_separated_tax  %>%
  select(1:8)
tax.table

##Filling up actual NAs and string "NA"s in the taxonomy table
filled_taxonomy <- fill_taxonomy(tax.table) ##apply the function to the taxonomy table
anyNA(filled_taxonomy) ##OK, no NAs now 
grep("^NA$", filled_taxonomy, value = T) ##OK, no "NA" strings now

##Now, to add the row names as "OTU1, OTU2, etc..." for phyloseq later on
filled_taxonomy_2<- filled_taxonomy %>%
  mutate(OTU = paste0("OTU", 1:nrow(filled_taxonomy))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") %>% ##Make OTU column into row names
  as.matrix() ##convert into matrix for phyloseq
filled_taxonomy_2

##OTU table 
otu_table <- counts[, -1]%>% #Excludes the first column (taxonomy)
  mutate(OTU = paste0("OTU", 1:nrow(counts))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") %>% ##Make OTU column into row names
  as.matrix() ##make into matrix so it is compatible with otu_table function from phyloseq
otu_table

# PHYLOSEQ OBJECT #####
OTU <-phyloseq::otu_table(otu_table, taxa_are_rows = TRUE)
TAX <-phyloseq::tax_table(filled_taxonomy_2)
phyloseq <- phyloseq(OTU, TAX)

#METADATA####
metadata <- read.csv('Metadata_Feedlot_CatchBasins.csv', 
                     check.names = F,
                     row.names = "sampleID")

metadata$SampleID<- rownames(metadata) #Making a SampleID column 
metadata$original_sample <- sub("c$", "", metadata$SampleID) #Making a column for the sample that both the metagenome and metatranscriptomes come from (match 2)


#HOST-FREE READS####
hostrem <- read.csv('HostRem_Reads_Feedlot_CatchBasins.csv')

#Change zymo- to zymo. as they are in metadata
hostrem[which(hostrem$SampleID == "Zymo-1a"), "SampleID"] <- "Zymo.1a"
hostrem[which(hostrem$SampleID == "Zymo-1b"), "SampleID"] <- "Zymo.1b"

#Merge with metadata
hostrem <- hostrem %>%
  select(SampleID, Hostrem_output_total_num_seqs)%>%
  left_join(metadata, by = "SampleID")%>%
  rename(HostFree_Reads=Hostrem_output_total_num_seqs)#Hostrem_output_total_num_seqs is the total number of host free reads
rownames(hostrem) <- hostrem$SampleID

##Making into phyloseq-compatible object
sampledata_phyloseq <- sample_data(metadata) ##use phyloseq function sample_data() to make metadata into phyloseq sample data object

## Add metadata to the phyloseq object####
data <- merge_phyloseq(phyloseq, sampledata_phyloseq)
sample_names(data) #OK

#COLOR PALETTES #####
# Feedlot
feedlot_palette <- c("1" = "#fcca46", 
                     "2" = "#fe7f2d", 
                     "3" = "#233d4d", 
                     "4"= "#3b9ab2", 
                     "5"= "#e1b6ff")
# Sample type
sample.type.palette <- c("Water" = "#4C72B0",
                         "Feces" = "brown") 
# Library type
gen.material.palette <- c("cDNA" = "#009E73",  
                          "DNA"  = "#CC79A7" )  

# Salmonella 
salmonella.palette <- c("positive"= "#fc8d62", "negative" = "#8da0cb")

# Sequencing batches
batch_palette <- c("no" = "#d19bac", 
                   "yes" = "#6a9c55") #sequencing batches



# PREPROCESSING ####
data #303961 taxa and 126 samples
sample_sums(data)

##ALL SAMPLES (samples + controls) ########
###Selecting only Bacteria/Archaea#######
data.bacteria <- subset_taxa(data, Domain=="Archaea" | Domain=="Bacteria")
data.bacteria #42541 taxa (dropped 191910 taxa), 126 samples = 120 samples, plus EBc2, EBc, EB, NTC2, Zymo.1a and Zymo.1b
###Selecting only viruses#######
data.viruses <- subset_taxa(data, Domain=="Viruses")
##Viruses have a different classification system, updating it here
taxanames_viruses <- c("Kingdom", "Realm", "Phylum", "Class", "Order", "Family", "Genus", "Species") 
colnames(data.viruses@tax_table) <- taxanames_viruses #replacing col names of the tax_table for new ones
colnames(data.viruses@tax_table) #OK taxonomy ranks
data.viruses #84188 taxa and 126 samples = 120 samples, plus EBc2, EBc, EB, NTC2, Zymo.1a and Zymo.1b

###Selecting only eukaryota ######
data.eukaryota <- subset_taxa(data, Domain=="Eukaryota")
colnames(data.eukaryota@tax_table) ##These are OK taxonomy ranks
data.eukaryota #177232 taxa and 126 samples 

##ONLY SAMPLES###
data.samples <- subset_samples(data, sample_type=="Water" | sample_type=="Feces")
data.samples <- prune_taxa(taxa_sums(data.samples) > 0, data.samples) 
data.samples #302969 taxa and 120 samples

#Bacteria/arcahaea
data.bacteria.samples <- subset_samples(data.bacteria, sample_type=="Water" | sample_type=="Feces")
data.bacteria.samples <- prune_taxa(taxa_sums(data.bacteria.samples) > 0, data.bacteria.samples) 
data.bacteria.samples #42002 taxa and 120 samples

#Viruses 
data.viruses.samples <- subset_samples(data.viruses, sample_type=="Water" | sample_type=="Feces")
data.viruses.samples <- prune_taxa(taxa_sums(data.viruses.samples) > 0, data.viruses.samples) 
data.viruses.samples #84040 taxa and 120 samples

#VEukarya
data.eukaryota.samples <- subset_samples(data.eukaryota, sample_type=="Water" | sample_type=="Feces")
data.eukaryota.samples <- prune_taxa(taxa_sums(data.eukaryota.samples) > 0, data.eukaryota.samples) 
data.eukaryota.samples #176927

##ZYMO MOCK COMMUNITIES #####
data.zymos <- subset_samples(data, grepl("Zymo", SampleID))
data.zymos <- prune_taxa(taxa_sums(data.zymos) > 0, data.zymos) 
data.zymos #19109 taxa and 2 samples

data.zymos.bacteria <- subset_taxa(data.zymos, Domain=="Archaea" | Domain=="Bacteria")
data.zymos.bacteria <- prune_taxa(taxa_sums(data.zymos.bacteria) > 0, data.zymos.bacteria) 
data.zymos.bacteria #13947 taxa and 2 samples 

data.zymos.fungi <- subset_taxa(data.zymos, Kingdom=="Fungi")
data.zymos.fungi <- prune_taxa(taxa_sums(data.zymos.fungi) > 0, data.zymos.fungi) 
data.zymos.fungi #1227 taxa and 2 samples 

###COMPOSITION (RA)#####
data.zymos #19109 taxa and 2 samples
##Only keeping bacteria and fungi 
data.zymos.bacteria.fungi <- merge_phyloseq(data.zymos.fungi, data.zymos.bacteria)
data.zymos.bacteria.fungi #15174 taxa and 2 samples

#Relative abundance
data.zymos.bacteria.fungi.ra <- transform_sample_counts(data.zymos.bacteria.fungi, 
                                                        function(x) x/sum(x)*100) ##Relative abundance 
#### GENUS LEVEL ####
data.zymos.bacteria.fungi.ra.genus <- tax_glom(data.zymos.bacteria.fungi.ra, taxrank = "Genus", NArm = F)

#Melt to long format at genus level
data.zymos.bacteria.fungi.ra.genus.melt <- psmelt(data.zymos.bacteria.fungi.ra.genus)

#Merge low abundance taxa
data.zymos.bacteria.fungi.ra.genus.filt <- merge_low_abundance_grouped_ra(data.zymos.bacteria.fungi.ra.genus, 
                                                                          variable = "sample_type",
                                                                          level = "Genus", 
                                                                          threshold = 0.2)
data.zymos.bacteria.fungi.ra.genus.filt #15 genera over 0.2% mean RA

#Factor order (low abundance last)
data.zymos.bacteria.fungi.ra.genus.filt.melt <- psmelt(data.zymos.bacteria.fungi.ra.genus.filt)%>%
  mutate(Genus = factor(Genus, 
                        levels = c(setdiff(Genus, 
                                           unique(grep("Others", Genus, value = TRUE))), 
                                   unique(grep("Others", Genus, value = TRUE)))))##Factoring the Phylum column so that "Others.." is the last category
levels(data.zymos.bacteria.fungi.ra.genus.filt.melt$Genus) ##ok

#Color palette at the genus level
zymo.genus.palette <- as.character(palette36.colors(15))
names(zymo.genus.palette) <- unique(data.zymos.bacteria.fungi.ra.genus.filt.melt$Genus)
zymo.genus.palette$'Staphylococcus' <- "darkblue" 
zymo.genus.palette$'Others <0.2% RA' <- 'grey90'

#Plot at the genus level
data.zymos.bacteria.fungi.ra.genus.plot <- ggplot(data.zymos.bacteria.fungi.ra.genus.filt.melt, aes(x=Sample, y= Abundance, fill = Genus)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", title = "Zymo Mock Communities") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =zymo.genus.palette) +
  guides(fill=guide_legend(title.position="top", nrow = 7))+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 30, colour = "black", face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
data.zymos.bacteria.fungi.ra.genus.plot

#What are the top genera?
top_genus_zymo <- data.zymos.bacteria.fungi.ra.genus.filt.melt %>%
  mutate(Abundance = round(Abundance, 2))%>%
  pivot_wider(
    id_cols = Genus,
    names_from = SampleID,
    values_from = Abundance)
top_genus_zymo


#### SPECIES LEVEL ####
data.zymos.bacteria.fungi.ra.species <- data.zymos.bacteria.fungi.ra
#Merge low abundance together
data.zymos.bacteria.fungi.ra.species.filt <- merge_low_abundance_grouped_ra(data.zymos.bacteria.fungi.ra.species, 
                                                                            variable = "sample_type",
                                                                            level = "Species", 
                                                                            threshold = 0.2)
data.zymos.bacteria.fungi.ra.species.filt #21 Species over 0.2% mean RA
#Melt to long format at the species level
data.zymos.bacteria.fungi.ra.species.filt.melt <- psmelt(data.zymos.bacteria.fungi.ra.species.filt)

#Factor ordering
data.zymos.bacteria.fungi.ra.species.filt.melt <- data.zymos.bacteria.fungi.ra.species.filt.melt%>%
  mutate(Species = factor(Species, levels = c(
    "Limosilactobacillus fermentum",
    "Listeria monocytogenes", 
    "unclassified Listeria", 
    "Pseudomonas aeruginosa", 
    "unclassified Pseudomonas", 
    "Bacillus spizizenii", 
    "unclassified Bacillus",
    "unclassified Bacilli",
    "unclassified Bacillales",
    "Escherichia coli", 
    "unclassified Escherichia", 
    "unclassified Enterobacteriaceae",
    "Salmonella enterica",
    "unclassified Salmonella", 
    "Enterococcus faecalis", 
    "unclassified Enterococcus",
    "Staphylococcus aureus",
    "unclassified Staphylococcus", 
    "Cryptococcus neoformans", 
    "Saccharomyces cerevisiae", 
    "unclassified Bacteria", 
    "Others <0.2% RA"
  )))
levels(data.zymos.bacteria.fungi.ra.species.filt.melt$Species) ##ok


#Color palette (based on species from the same genera)
zymo_genus_base_colors <- zymo.genus.palette
#Make hues based on families within each genus
palette_zymo_genus_df <- data.zymos.bacteria.fungi.ra.species.filt.melt %>% 
  distinct(Genus, Species) %>%
  group_by(Genus) %>%  
  arrange(Species) %>%   
  mutate(
    base_color = zymo_genus_base_colors[Genus],
    shade = seq(-0.2, 0.2, length.out = n()),  # slightly wider range helps
    color = darken(base_color, amount = shade)
  ) %>%
  ungroup()
palette_zymo_genus_df

#Set up final palette
palette_zymo_genus <- setNames(
  palette_zymo_genus_df$color,
  palette_zymo_genus_df$Species)
palette_zymo_genus
palette_zymo_genus$'Others <0.2% RA' <- 'grey90'

#Plot at the species level
data.zymos.bacteria.fungi.ra.species.plot <- ggplot(data.zymos.bacteria.fungi.ra.species.filt.melt, aes(x=Sample, y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", title = "Zymo Mock Communities") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = palette_zymo_genus) +
  guides(fill=guide_legend(title.position="top", nrow = 7))+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 30, colour = "black", face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
data.zymos.bacteria.fungi.ra.species.plot

#What are the top species?
top_species_zymo <- data.zymos.bacteria.fungi.ra.species.filt.melt %>%
  group_by(Family, Genus, Species) %>%
  summarise(
    `Mean Relative Abundance (%)` = mean(Abundance, na.rm = TRUE),
    `Standard Deviation` = sd(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Family, Genus) %>%
  arrange(Family, Genus, desc(`Mean Relative Abundance (%)`), .by_group = FALSE) %>%
  mutate(
    `Mean Relative Abundance (%) ± SD` = paste0(
      round(`Mean Relative Abundance (%)`, 2),
      " ± ",
      round(`Standard Deviation`, 3)
    )
  ) %>%
  ungroup() %>%
  select(Genus, Species, `Mean Relative Abundance (%) ± SD`)
top_species_zymo

#####SUPPLEMENTARY TABLE 4####
stable4 <- top_genus_zymo
write_xlsx(stable4, 
           "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryTable4.xlsx")

#####SUPPLEMENTARY FIGURE 4######
sfigure4 <- data.zymos.bacteria.fungi.ra.genus.plot
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryFigure4.png", 
       sfigure4, 
       device = "png", 
       width = 10, 
       height =10, 
       dpi = 500,
       bg = "white"
       )

##EBs and NTC#####
data.negatives <- subset_samples(data, grepl("EB|NTC", SampleID))
data.negatives <- prune_taxa(taxa_sums(data.negatives) > 0, data.negatives)
data.negatives #87 taxa and 4 samples

data.negatives.bacteria <- subset_taxa(data.negatives,  Domain=="Archaea" | Domain=="Bacteria")
data.negatives.bacteria <- prune_taxa(taxa_sums(data.negatives.bacteria) > 0, data.negatives.bacteria) 
data.negatives.bacteria #78 taxa and 4 samples  

###NEGATIVE CONTROLS COMPOSITION#####
data.negatives.bacteria #78 taxa and 4 samples 
data.negatives.bacteria.ra <- transform_sample_counts(data.negatives.bacteria, function(x) x/sum(x)*100) ##Relative abundance 

#Glom to genus level 
data.negatives.bacteria.ra.genus <- tax_glom(data.negatives.bacteria.ra, taxrank = "Genus", NArm = F)

#Melt to long format at genus level
data.negatives.bacteria.ra.genus.melt <- psmelt(data.negatives.bacteria.ra.genus)

#What are the top species?
top_genera_negatives <- data.negatives.bacteria.ra.genus.melt %>%
  group_by(Genus) %>%
  summarise(mean_relative_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_relative_abundance))%>%
  head(n = 30)
top_genera_negatives

#Plot
data.negatives.bacteria.ra.genus.filt <- merge_low_abundance_grouped_ra(data.negatives.bacteria.ra.genus, variable = "sample_type",
                                                                        level = "Genus", 
                                                                        threshold = 1.5)
data.negatives.bacteria.ra.genus.filt #14 genera over 1.5% mean RA
data.negatives.bacteria.ra.genus.filt.melt <- psmelt(data.negatives.bacteria.ra.genus.filt)%>%
  mutate(Genus = factor(Genus, 
                        levels = c(setdiff(Genus, 
                                           unique(grep("Others", Genus, value = TRUE))), 
                                   unique(grep("Others", Genus, value = TRUE)))))##Factoring the Phylum column so that "Others.." is the last category
levels(data.negatives.bacteria.ra.genus.filt.melt$Genus) ##ok

#Palette
negatives.genus.palette <- as.character(paletteer_d("ggsci::category20b_d3")) 
negatives.genus.palette <- setNames(negatives.genus.palette, 
                                    levels(data.negatives.bacteria.ra.genus.filt.melt$Genus))
negatives.genus.palette$"Others <1.5% RA" <- "grey80"
negatives.genus.palette$"unclassified Bacteria" <- "#FF99BFFF"

#Plot
data.negatives.bacteria.ra.genus.plot <- ggplot(data.negatives.bacteria.ra.genus.filt.melt, aes(x=Sample, 
                                                                                                y= Abundance, fill = Genus)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", title = "Negative Controls") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = negatives.genus.palette) +
  guides(fill=guide_legend(title.position="top", nrow =7))+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 30, colour = "black", face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
data.negatives.bacteria.ra.genus.plot


##BACTERIA/ARCHAEA ONLY#######
# some QC checks of the "classified" reads per samples
min(sample_sums(data.bacteria)) # 4 (NTC2
max(sample_sums(data.bacteria)) # 42,044,077 (Zymo.1b)
mean(sample_sums(data.bacteria)) #10,274,560
median(sample_sums(data.bacteria)) # 9,540,090
sort(sample_sums(data.bacteria)) # 4, 27, 37, 48 (these lowest ones are NTC2, EB, EBc, EBc2)

###ONLY SAMPLES#####
data.bacteria.samples <- subset_samples(data.bacteria, sample_type=="Water" | sample_type=="Feces")
data.bacteria.samples <- prune_taxa(taxa_sums(data.bacteria.samples) > 0, data.bacteria.samples) 
data.bacteria.samples #42002 taxa, 120 samples (dropped NTC2, EB, EBc, EBc2, Zymo.1a, Zymo.1b)

##QC checks again
min(sample_sums(data.bacteria.samples)) #5,561,216(F2W03c)
max(sample_sums(data.bacteria.samples)) #23,115,920 (F2W01c) 
mean(sample_sums(data.bacteria.samples)) #10,095,556
median(sample_sums(data.bacteria.samples)) #9,569,330 
sort(sample_sums(data.bacteria.samples)) 


#KRAKEN2 CLASSIFICATION PERCENTAGES####
#Join unclassified counts with metadata
unclassified_counts_metadata <- unclassified_counts%>%
  #Join with metadata
  left_join(metadata, by = "SampleID")%>% 
  mutate(PercentClassified = 100 - PercentUnclassified)%>%
  filter(!grepl("Zymo|EB|NTC", SampleID)) #Pick only samples
nrow(unclassified_counts_metadata) #120 samples

#Factor order 
unclassified_counts_metadata$sample_type <- factor(unclassified_counts_metadata$sample_type, levels = c("Feces", "Water"))
unclassified_counts_metadata$gen_material <- factor(unclassified_counts_metadata$gen_material, levels = c("DNA", "cDNA"))

##DNA vs cDNA(RNA) faceted by Sample type - Not significant#####
kraken2_classification_percentage.DNAvscDNA <- ggplot(unclassified_counts_metadata%>%arrange(SampleID),
                                           aes(x = gen_material,
                                               y= PercentClassified, 
                                               color = gen_material, 
                                               fill = gen_material)) +
  theme_bw() +
  labs(y= "Percentage (%)\n of reads", 
       color = "Library Type", fill = "Library Type",
       title = "KRAKEN2 CLASSIFICATION") +
  facet_grid(~sample_type, scales = "free",  labeller = as_labeller(c("Feces" = "FECES", "Water" = "CATCH BASIN"))) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  scale_fill_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.title = element_text(size = 32, face = "bold"),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc(method = "wilcox_test",
           label = "Wilcoxon, p = {p}",
           hide.ns = TRUE,
           method.args = list(paired = T),
           step.increase = 0.08,
           label.size = 6,
           tip.length = 0.02)
kraken2_classification_percentage.DNAvscDNA

#Checking Stats
unclassified_counts_metadata.water <- subset(unclassified_counts_metadata, sample_type == "Water")
#Percentages - CB
unclassified_counts_metadata.water %>%
  group_by(gen_material) %>%
  summarise(mean_percent_classified = mean(PercentClassified, na.rm = TRUE),
            min_percent_classified = min(PercentClassified, na.rm = TRUE),
            max_percent_classified = max(PercentClassified, na.rm = TRUE))
# gen_material mean_percent_classified min_percent_classified max_percent_classified
# <chr>                          <dbl>                  <dbl>                  <dbl>
#  DNA                             39.5                   33.2                   48.9
#  cDNA                            43.3                   35.4                   56.5

#Wilcox test
wilcox_test(unclassified_counts_metadata.water%>%arrange(SampleID), 
            PercentClassified ~ gen_material, 
            paired = T) #s, p = 0.000488

#Percentages - Feces
unclassified_counts_metadata.feces <- subset(unclassified_counts_metadata, sample_type == "Feces")
unclassified_counts_metadata.feces %>%
  group_by(gen_material) %>%
  summarise(mean_percent_classified = mean(PercentClassified, na.rm = TRUE),
            min_percent_classified = min(PercentClassified, na.rm = TRUE),
            max_percent_classified = max(PercentClassified, na.rm = TRUE))
# gen_material mean_percent_classified min_percent_classified max_percent_classified
# <chr>                          <dbl>                  <dbl>                  <dbl>
# DNA                             32.5                   28.4                   42.3
# cDNA                            32.3                   28.2                   40.9

#Wilcox test 
wilcox_test(unclassified_counts_metadata.feces%>%arrange(SampleID), 
            PercentClassified ~ gen_material, 
            paired = T) #ns, p = 0.131


##Water vs Feces faceted by cDNA and DNA- Significant#####
kraken2_classification_percentage.WvsF<- ggplot(unclassified_counts_metadata, 
                                                     aes(x = sample_type, 
                                                         y= PercentClassified, 
                                                         color = sample_type, 
                                                         fill = sample_type)) +
  theme_bw() +
  labs(y= "Percentage (%)\n of reads", 
       color = "Sample Type", fill = "Sample Type",
       title = "KRAKEN2 CLASSIFICATION") +
  facet_grid(~gen_material, scales = "free",  
             labeller = as_labeller(c("cDNA" = "RNA (cDNA)", "DNA" = "DNA"))) +
  # facet_wrap(~gen_material) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  scale_fill_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch Basin"))+
  scale_color_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch Basin"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        plot.title = element_text(size = 32, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))+
  geom_pwc(method = "wilcox_test",
           label = "Wilcoxon, p = {p}",
           hide.ns = TRUE,
           step.increase = 0.08,
           label.size = 6,
           tip.length = 0.02)
kraken2_classification_percentage.WvsF

#Checking Stats
unclassified_counts_metadata.DNA <- subset(unclassified_counts_metadata, gen_material == "DNA")
#Percentages - DNA
unclassified_counts_metadata.DNA %>%
  group_by(sample_type) %>%
  summarise(mean_percent_classified = mean(PercentClassified, na.rm = TRUE),
            min_percent_classified = min(PercentClassified, na.rm = TRUE),
            max_percent_classified = max(PercentClassified, na.rm = TRUE))
# sample_type mean_percent_classified min_percent_classified max_percent_classified
# Feces                          32.5                   28.4                   42.3
# Water                          39.5                   33.2                   48.9

#Wilcox test
wilcox_test(unclassified_counts_metadata.DNA, PercentClassified ~ sample_type) #S, p = 0.00000136

#Percentages - cDNA
unclassified_counts_metadata.cDNA <- subset(unclassified_counts_metadata, gen_material == "cDNA")
unclassified_counts_metadata.cDNA %>%
  group_by(sample_type) %>%
  summarise(mean_percent_classified = mean(PercentClassified, na.rm = TRUE),
            min_percent_classified = min(PercentClassified, na.rm = TRUE),
            max_percent_classified = max(PercentClassified, na.rm = TRUE))
# sample_type mean_percent_classified min_percent_classified max_percent_classified
# Feces                          32.3                   28.2                   40.9
# Water                          43.3                   35.4                   56.5

#Wilcox test 
wilcox_test(unclassified_counts_metadata.cDNA, PercentClassified ~ sample_type) #S, p = 7.22e-10

###SUPPLEMENTARY FIGURE 3#####
kraken2_classification_percentage.WvsF
kraken2_classification_percentage.DNAvscDNA

sfigure3 <- cowplot::plot_grid(
  kraken2_classification_percentage.WvsF + theme(axis.title.y = element_text(size = 20),
                                                 plot.title = element_blank()),
  kraken2_classification_percentage.DNAvscDNA + theme(axis.title.y = element_text(size = 20),
                                                      plot.title = element_blank()),
  ncol = 1, labels = "AUTO", label_size = 22) +
  #labs(title = "KRAKEN2 READ CLASSIFICATION") +
  theme(
      #plot.title = element_text(size = 32, face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
sfigure3 
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryFigure3.png", 
       plot = sfigure3, 
       device = "png", 
       dpi = 600,
       width = 10, height =10)


#COMPARING OTU COUNTS#######
sample.sums <- sample_sums(data.bacteria.samples) #making a sample sums object
data.bacteria.samples.df <- cbind(data.bacteria.samples@sam_data, sample.sums) #combining sample sums with metadata
data.bacteria.samples.df 
data.bacteria.samples.df$sampleID <- rownames(data.bacteria.samples.df) ##making a sampleID column

#Factor order 
data.bacteria.samples.df$sample_type <- factor(data.bacteria.samples.df$sample_type, levels = c("Feces", "Water"))
data.bacteria.samples.df$gen_material <- factor(data.bacteria.samples.df$gen_material, levels = c("DNA", "cDNA"))

##cDNA vs DNA faceted by water and feces- N.S.####
otu_counts_cDNAvDNA.WandF<- ggplot(data.bacteria.samples.df, aes(x = gen_material, y= sample.sums, color = gen_material, fill = gen_material)) +
  theme_bw() +
  labs(y= "OTUs per Sample", color = "Library Type", fill = "Library Type") +
  facet_grid(~sample_type, scales = "free",  labeller = as_labeller(c("Feces" = "FECES", "Water" = "CATCH BASINS"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
  scale_fill_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) 
otu_counts_cDNAvDNA.WandF

#Checking Stats
water.data.bacteria.samples <- subset(data.bacteria.samples.df, sample_type == "Water")
wilcox_test(water.data.bacteria.samples, sample.sums~gen_material) #n.s. p = 0.291
feces.data.bacteria.samples <- subset(data.bacteria.samples.df, sample_type == "Feces")
wilcox_test(feces.data.bacteria.samples, sample.sums~gen_material) #n.s. p = 0.724

##Water vs Feces faceted by cDNA and DNA- Significant#####
sequencing_depth_WvF.cDNAandDNA<- ggplot(data.bacteria.samples.df, aes(x = sample_type, y= sample.sums, color = sample_type, fill = sample_type)) +
  theme_bw() +
  labs(y= "OTUs per Sample", color = "Sample Type", fill = "Sample Type") +
  facet_grid(~gen_material, scales = "free",  
             labeller = as_labeller(c("cDNA" = "RNA (cDNA)", "DNA" = "DNA"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
  scale_fill_manual(values = sample.type.palette)+
  scale_color_manual(values = sample.type.palette)+
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))+
  geom_pwc(method = "wilcox_test",
           label = "Wilcoxon, p = {p}",
           hide.ns = TRUE,
           step.increase = 0.08,
           label.size = 6,
           tip.length = 0.02)
sequencing_depth_WvF.cDNAandDNA

#Checking Stats
cDNA.data.bacteria.samples <- subset(data.bacteria.samples.df, gen_material == "cDNA")
wilcox_test(cDNA.data.bacteria.samples, sample.sums~ sample_type) #s, p = 6.62e-06
DNA.data.bacteria.samples <- subset(data.bacteria.samples.df, gen_material == "DNA")
wilcox_test(DNA.data.bacteria.samples, sample.sums~ sample_type) #S p = 6.75e-05

##Feedlot ####
sequencing_depth_feedlot <- ggplot(data.bacteria.samples.df, aes(x = feedlot, y= sample.sums, color = factor(feedlot), fill = factor(feedlot))) +
  theme_bw() +
  labs(title = "Bacterial-Archaeal OTUs" , y= "OTUs per Sample", color = "Feedlot", fill = "Feedlot" ) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3,  shape = 18) +
  scale_fill_manual(values = feedlot_palette) +
  scale_color_manual(values = feedlot_palette) +
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.title = element_text(size = 32, face = "bold"),
        plot.title.position = "plot",
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))+
geom_pwc(method = "wilcox_test", p.adjust.method = "BH",label = "p = {p.adj}",
         hide.ns = TRUE,
         step.increase = 0.08,
         label.size = 5,
         tip.length = 0.02)
sequencing_depth_feedlot

#Checking Stats
wilcox_test(data.bacteria.samples.df, sample.sums~feedlot, p.adjust.method = "BH") #Significant between 1 and 4 (0.042), and between 3 and 4 (p =0.042)

##Checking for batch effects - reseq vs first run ####
sequencing_depth_batches <- ggplot(data.bacteria.samples.df, aes(x = re_sequenced, y= sample.sums, color = re_sequenced, fill = re_sequenced)) +
  theme_bw() +
  labs(title = "Bacterial-Archaeal OTUs", y= "OTUs per Sample",color = "Re-Sequenced", fill = "Re-Sequenced") +
  geom_boxplot(alpha = 0.1) +
  #geom_boxplot(aes(group = factor(feedlot)), alpha = 0.1) +
  geom_point(size = 3,  shape = 18) +
  scale_fill_manual(values = batch_palette) +
  scale_color_manual(values = batch_palette) +
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.title = element_text(size = 32, face = "bold"),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))+
  geom_pwc(method = "wilcox_test",
           label = "Wilcoxon, p = {p}",
           hide.ns = TRUE,
           step.increase = 0.08,
           label.size = 6,
           tip.length = 0.02)
sequencing_depth_batches

#Checking Stats
wilcox_test(data.bacteria.samples.df, sample.sums~re_sequenced) #NS (p = 0.669)


#TSS NORMALIZATION (RA)####
any(sample_sums(data.bacteria.samples)== 0) ## no samples with 0 OTUs
data.bacteria.samples.tss <- transform_sample_counts(data.bacteria.samples, function(x) x/sum(x)*100) ##Relative abundance 

#Metadata DF
data.bacteria.samples.tss.df <- as(data.bacteria.samples.tss@sam_data,"data.frame") # make DF from metadata
data.bacteria.samples.tss.df

##TAX GLOMMING #####
##PHYLUM
data.bacteria.samples_phylum.ra <- tax_glom(data.bacteria.samples.tss, taxrank = "Phylum", NArm = F) 
data.bacteria.samples_phylum.ra #6775 taxa and 120 samples
#Are there duplicates? 
length(unique(phyloseq::tax_table(data.bacteria.samples_phylum.ra)[, "Phylum"])) #5875 taxa (so No duplicates)

Unknown_phylum_abundance <- data.bacteria.samples_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Phylum, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_phylum_abundance ##1.32% abundance by Unknown Phyla

Unclassified_phylum_abundance <- data.bacteria.samples_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_phylum_abundance ##8.74% abundance by Unclassified Phyla

Classified_phylum_abundance <- data.bacteria.samples_phylum.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Phylum, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_phylum_abundance ##89.9% abundance by Classified Phyla

##Checking on excel
write.csv(data.bacteria.samples_phylum.ra@otu_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/phylum_otus.csv")
write.csv(data.bacteria.samples_phylum.ra@tax_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/phylum_taxa.csv")  

#How many unclassified?
data.bacteria.samples_phylum.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(data.bacteria.samples_phylum.ra)[, "Phylum"]),
  data.bacteria.samples_phylum.ra)
data.bacteria.samples_phylum.unclassified.ra #9 unclassified Phyla

#How many unknown?
data.bacteria.samples_phylum.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(data.bacteria.samples_phylum.ra)[, "Phylum"]),
  data.bacteria.samples_phylum.ra)
data.bacteria.samples_phylum.unknown.ra #5747 "unknown" Phyla

#Keep just classified Phyla
data.bacteria.samples_phylum.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(data.bacteria.samples_phylum.ra)[, "Phylum"]),
  data.bacteria.samples_phylum.ra)
data.bacteria.samples_phylum.classified.ra ##125 classified (not unknown or unclassified) Phyla

##CLASS
data.bacteria.samples_class.ra <- tax_glom(data.bacteria.samples.tss, taxrank = "Class", NArm = F) 
data.bacteria.samples_class.ra #8,452 classes (120 samples)
#Are there duplicates? 
length(unique(phyloseq::tax_table(data.bacteria.samples_class.ra)[, "Class"])) #8451 taxa (so one duplicate)

Unknown_class_abundance <- data.bacteria.samples_class.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Class, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_class_abundance #2.52% Abundance by Unknown classes

Unclassified_class_abundance <- data.bacteria.samples_class.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Class, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_class_abundance ##11.4% Abundance by Unclassified Classes

Classified_class_abundance <- data.bacteria.samples_class.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Class, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_class_abundance ##86.1% Abundance by Classified classes

##Checking on excel
write.csv(data.bacteria.samples_class.ra@otu_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/class_otus.csv")
write.csv(data.bacteria.samples_class.ra@tax_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/class_taxa.csv") 


#How many unclassified?
data.bacteria.samples_class.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(data.bacteria.samples_class.ra)[, "Class"]),
  data.bacteria.samples_class.ra)
data.bacteria.samples_class.unclassified.ra #81 unclassified classes

#How many unknown?
data.bacteria.samples_class.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(data.bacteria.samples_class.ra)[, "Class"]),
  data.bacteria.samples_class.ra)
data.bacteria.samples_class.unknown.ra #8215 "unknown" classes

#Keep just classified Classes
data.bacteria.samples_class.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(data.bacteria.samples_class.ra)[, "Class"]),
  data.bacteria.samples_class.ra)
data.bacteria.samples_class.classified.ra #156 classified classes

##ORDER
data.bacteria.samples_order.ra <- tax_glom(data.bacteria.samples.tss, taxrank = "Order", NArm = F) 
data.bacteria.samples_order.ra #9297 orders
#Are there duplicates? 
length(unique(phyloseq::tax_table(data.bacteria.samples_order.ra)[, "Order"])) #9293 taxa (4 duplicates)
order_taxa_vec <- as.character(phyloseq::tax_table(data.bacteria.samples_order.ra)[, "Order"])
unique(order_taxa_vec[duplicated(order_taxa_vec)]) 
#"Candidatus Fermentimicrarchaeales", "Candidatus Cenarchaeales", "Mycoplasmoidales",  "Candidatus Moduliflexales" 

Unknown_order_abundance <- data.bacteria.samples_order.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Order, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_order_abundance ##3.94% abundance by Unknown Orders

Unclassified_order_abundance <- data.bacteria.samples_order.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Order, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_order_abundance ##14% abundance by Unclassified Orders

Classified_order_abundance <- data.bacteria.samples_order.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Order, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_order_abundance ##82.1% abundance by Classified orders

#Checking on excel
write.csv(data.bacteria.samples_order.ra@otu_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/order_otus.csv")
write.csv(data.bacteria.samples_order.ra@tax_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/order_taxa.csv") 

#How many unclassified?
data.bacteria.samples_order.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(data.bacteria.samples_order.ra)[, "Order"]),
  data.bacteria.samples_order.ra)
data.bacteria.samples_order.unclassified.ra #149 unclassified orders

#How many unknown?
data.bacteria.samples_order.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(data.bacteria.samples_order.ra)[, "Order"]),
  data.bacteria.samples_order.ra)
data.bacteria.samples_order.unknown.ra #8804 "unknown" orders

#Keep just classified Orders
data.bacteria.samples_order.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(data.bacteria.samples_order.ra)[, "Order"]),
  data.bacteria.samples_order.ra)
data.bacteria.samples_order.classified.ra #344 classified orders
length(unique(phyloseq::tax_table(data.bacteria.samples_order.classified.ra)[, "Order"])) ##338 classified orders (unique - without duplicates)

##FAMILY
data.bacteria.samples_family.ra <- tax_glom(data.bacteria.samples.tss, taxrank = "Family", NArm = F) 
data.bacteria.samples_family.ra #10,529 families
#Are there duplicates? 
length(unique(phyloseq::tax_table(data.bacteria.samples_family.ra)[, "Family"])) #10523 taxa (6 duplicates)
family_taxa_vec <- as.character(phyloseq::tax_table(data.bacteria.samples_family.ra)[, "Family"])
unique(family_taxa_vec[duplicated(family_taxa_vec)]) 
#"Candidatus Fermentimicrarchaeales", 
# "Candidatus Cenarchaeales", "unclassified Mycoplasmoidales", "Mycoplasmoidales"
# Metamycoplasmataceae", "Candidatus Moduliflexaceae"

Unknown_family_abundance <- data.bacteria.samples_family.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_family_abundance #7.72% abundance by Unknown Families

Unclassified_family_abundance <- data.bacteria.samples_family.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_family_abundance ##16.9% abundance by Unclassified Families

Classified_family_abundance <- data.bacteria.samples_family.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Family, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_family_abundance ##75.4% abundance by Classified Families

#Checking on excel
write.csv(data.bacteria.samples_family.ra@otu_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/family_otus.csv")
write.csv(data.bacteria.samples_family.ra@tax_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/family_taxa.csv") 

#How many unclassified?
data.bacteria.samples_family.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(data.bacteria.samples_family.ra)[, "Family"]),
  data.bacteria.samples_family.ra)
data.bacteria.samples_family.unclassified.ra #296 unclassified families
length(unique(phyloseq::tax_table(data.bacteria.samples_family.unclassified.ra)[, "Family"])) ##295 classified families (unique - without duplicates)

#How many unknown?
data.bacteria.samples_family.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(data.bacteria.samples_family.ra)[, "Family"]),
  data.bacteria.samples_family.ra)
data.bacteria.samples_family.unknown.ra #9420 "unknown" families
length(unique(phyloseq::tax_table(data.bacteria.samples_family.unknown.ra)[, "Family"]))#9420 "unknown" taxa (unique - without duplicates)

#Keep just classified Families
data.bacteria.samples_family.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(data.bacteria.samples_family.ra)[, "Family"]),
  data.bacteria.samples_family.ra)
data.bacteria.samples_family.classified.ra #813 classified families
length(unique(phyloseq::tax_table(data.bacteria.samples_family.classified.ra)[, "Family"]))#808 classified families (unique - without duplicates)

##GENUS 
data.bacteria.samples_genus.ra <- tax_glom(data.bacteria.samples.tss, taxrank = "Genus", NArm = F) 
data.bacteria.samples_genus.ra #14574 genera
#Are there duplicates? 
length(unique(phyloseq::tax_table(data.bacteria.samples_genus.ra)[, "Genus"])) #14498 taxa (63 duplicates)
genus_taxa_vec <- as.character(phyloseq::tax_table(data.bacteria.samples_genus.ra)[, "Genus"])
unique(genus_taxa_vec[duplicated(genus_taxa_vec)]) #47 duplicated unique ones

Unknown_genus_abundance <- data.bacteria.samples_genus.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_genus_abundance ##13.4%  abundance by unknown genera

Unclassified_genus_abundance <- data.bacteria.samples_genus.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_genus_abundance ##20.6% abundance by unclassified genera

Classified_genus_abundance <- data.bacteria.samples_genus.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Genus, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_genus_abundance ##66.0% abundance by Classified Genera

#Checking on excel
write.csv(data.bacteria.samples_genus.ra@otu_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/genus_otus.csv")
write.csv(data.bacteria.samples_genus.ra@tax_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/genus_taxa.csv") 

#How many unclassified?
data.bacteria.samples_genus.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(data.bacteria.samples_genus.ra)[, "Genus"]),
  data.bacteria.samples_genus.ra)
data.bacteria.samples_genus.unclassified.ra #687 unclassified genera
length(unique(phyloseq::tax_table(data.bacteria.samples_genus.unclassified.ra)[, "Genus"])) ##684 unclassified genera (unique - without duplicates)

#How many unknown?
data.bacteria.samples_genus.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(data.bacteria.samples_genus.ra)[, "Genus"]),
  data.bacteria.samples_genus.ra)
data.bacteria.samples_genus.unknown.ra #10,323 "unknown" genera
length(unique(phyloseq::tax_table(data.bacteria.samples_genus.unknown.ra)[, "Genus"])) ##10323 unknown genera (unique - without duplicates)


#Keep just classified Genera
data.bacteria.samples_genus.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(data.bacteria.samples_genus.ra)[, "Genus"]),
  data.bacteria.samples_genus.ra)
data.bacteria.samples_genus.classified.ra #3564 classified genera
length(unique(phyloseq::tax_table(data.bacteria.samples_genus.classified.ra)[, "Genus"])) ##3491 classified genera (unique - without duplicates)


##Species - The most granular level is species
data.bacteria.samples_species.ra <- tax_glom(data.bacteria.samples.tss, taxrank = "Species", NArm = F) 
data.bacteria.samples_species.ra #42,002 species (same as data.bacteria.samples.tss)
#Are there duplicates? 
length(unique(phyloseq::tax_table(data.bacteria.samples_species.ra)[, "Species"])) #41714 taxa (there are duplicates)
species_taxa_vec <- as.character(phyloseq::tax_table(data.bacteria.samples_species.ra)[, "Species"])
unique(species_taxa_vec[duplicated(species_taxa_vec)]) #209 unique duplicated ones


Unclassified_species_abundance <- data.bacteria.samples_species.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Species, ignore.case = TRUE))%>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU across samples
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_species_abundance #28.9% Abundance by Unclassified Species 

Classified_species_abundance <- data.bacteria.samples_species.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Species, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_species_abundance ##71.1% Abundance by Classified Species

#Checking on excel
write.csv(data.bacteria.samples_species.ra@otu_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/species_otus.csv")
write.csv(data.bacteria.samples_species.ra@tax_table, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/species_taxa.csv") 

#How many unclassified?
data.bacteria.samples_species.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(data.bacteria.samples_species.ra)[, "Species"]),
  data.bacteria.samples_species.ra)
data.bacteria.samples_species.unclassified.ra #2234 unclassified Species 
length(unique(phyloseq::tax_table(data.bacteria.samples_species.unclassified.ra)[, "Species"])) ## 2217 unclassified species (unique - without duplicates)

#Keep just classified Species
data.bacteria.samples_species.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(data.bacteria.samples_species.ra)[, "Species"]),
  data.bacteria.samples_species.ra)
data.bacteria.samples_species.classified.ra #39766 classified species
length(unique(phyloseq::tax_table(data.bacteria.samples_species.classified.ra)[, "Species"])) ## 39495 classified species (unique - without duplicates)


##Putting all together 
classified_unknown_table <- data.frame(
  "Taxonomic level" = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  # Total taxa at each level
  "Number of Taxa" = c(
    length(taxa_names(data.bacteria.samples_phylum.ra)),
    length(taxa_names(data.bacteria.samples_class.ra)), 
    length(taxa_names(data.bacteria.samples_order.ra)), 
    length(taxa_names(data.bacteria.samples_family.ra)), 
    length(taxa_names(data.bacteria.samples_genus.ra)),
    length(taxa_names(data.bacteria.samples_species.ra))),
  
  # Total unique taxa at each level
  "Number of Unique Taxa" = c(
    length(unique(phyloseq::tax_table(data.bacteria.samples_phylum.ra)[, "Phylum"])),
    length(unique(phyloseq::tax_table(data.bacteria.samples_class.ra)[, "Class"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_order.ra)[, "Order"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_family.ra)[, "Family"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_genus.ra)[, "Genus"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_species.ra)[, "Species"]))),
  
  # Unclassified taxa
  "Number of Unclassified Taxa" = c(
    length(taxa_names(data.bacteria.samples_phylum.unclassified.ra)),
    length(taxa_names(data.bacteria.samples_class.unclassified.ra)), 
    length(taxa_names(data.bacteria.samples_order.unclassified.ra)), 
    length(taxa_names(data.bacteria.samples_family.unclassified.ra)), 
    length(taxa_names(data.bacteria.samples_genus.unclassified.ra)), 
    length(taxa_names(data.bacteria.samples_species.unclassified.ra))),
  
  # Unique unclassified taxa
  "Number of Unique Unclassified Taxa" = c(
    length(unique(phyloseq::tax_table(data.bacteria.samples_phylum.unclassified.ra)[, "Phylum"])),
    length(unique(phyloseq::tax_table(data.bacteria.samples_class.unclassified.ra)[, "Class"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_order.unclassified.ra)[, "Order"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_family.unclassified.ra)[, "Family"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_genus.unclassified.ra)[, "Genus"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_species.unclassified.ra)[, "Species"]))),
  
  # Unknown taxa
  "Number of Unknown Taxa" = c(
    length(taxa_names(data.bacteria.samples_phylum.unknown.ra)),
    length(taxa_names(data.bacteria.samples_class.unknown.ra)), 
    length(taxa_names(data.bacteria.samples_order.unknown.ra)), 
    length(taxa_names(data.bacteria.samples_family.unknown.ra)), 
    length(taxa_names(data.bacteria.samples_genus.unknown.ra)), 
    NA),
  
  # Unique unknown taxa (no species level)
  "Number of Unique Unknown Taxa" = c(
    length(unique(phyloseq::tax_table(data.bacteria.samples_phylum.unknown.ra)[, "Phylum"])),
    length(unique(phyloseq::tax_table(data.bacteria.samples_class.unknown.ra)[, "Class"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_order.unknown.ra)[, "Order"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_family.unknown.ra)[, "Family"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_genus.unknown.ra)[, "Genus"])), 
    NA),
  
  # Classified taxa
  "Number of Classified Taxa" = c(
    length(taxa_names(data.bacteria.samples_phylum.classified.ra)),
    length(taxa_names(data.bacteria.samples_class.classified.ra)), 
    length(taxa_names(data.bacteria.samples_order.classified.ra)), 
    length(taxa_names(data.bacteria.samples_family.classified.ra)), 
    length(taxa_names(data.bacteria.samples_genus.classified.ra)),
    length(taxa_names(data.bacteria.samples_species.classified.ra))),
  
  # Unique classified taxa
  "Number of Unique Classified Taxa" = c(
    length(unique(phyloseq::tax_table(data.bacteria.samples_phylum.classified.ra)[, "Phylum"])),
    length(unique(phyloseq::tax_table(data.bacteria.samples_class.classified.ra)[, "Class"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_order.classified.ra)[, "Order"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_family.classified.ra)[, "Family"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_genus.classified.ra)[, "Genus"])), 
    length(unique(phyloseq::tax_table(data.bacteria.samples_species.classified.ra)[, "Species"]))),
  
  # Relative abundances
  "Mean Relative Abundance (%) of Unclassified Taxa Across Samples" = c(
    Unclassified_phylum_abundance$Unclassified_sum, 
    Unclassified_class_abundance$Unclassified_sum,
    Unclassified_order_abundance$Unclassified_sum,
    Unclassified_family_abundance$Unclassified_sum,
    Unclassified_genus_abundance$Unclassified_sum,
    Unclassified_species_abundance$Unclassified_sum
  ),
  "Mean Relative Abundance (%) of Unknown Taxa Across Samples" = c(
    Unknown_phylum_abundance$Unknown_sum, 
    Unknown_class_abundance$Unknown_sum, 
    Unknown_order_abundance$Unknown_sum,
    Unknown_family_abundance$Unknown_sum,
    Unknown_genus_abundance$Unknown_sum,
    NA),
  "Mean Relative Abundance (%) of Classified Taxa Across Samples" = c(
    Classified_phylum_abundance$Classified_sum,
    Classified_class_abundance$Classified_sum, 
    Classified_order_abundance$Classified_sum,
    Classified_family_abundance$Classified_sum,
    Classified_genus_abundance$Classified_sum,
    Classified_species_abundance$Classified_sum),
  check.names = FALSE
) %>%
  mutate(`Percentage of Classified Taxa` =
           (`Number of Unique Classified Taxa` / `Number of Unique Taxa`) * 100)
classified_unknown_table


####SUPPLEMENTARY TABLE 3####
stable3 <- classified_unknown_table%>%
  mutate(`Mean Relative Abundance (%) of Classified Taxa Across Samples` = 
           round(`Mean Relative Abundance (%) of Classified Taxa Across Samples`, 2), 
         `Percentage of Classified Taxa` = 
           round(`Percentage of Classified Taxa`, 2)
         )%>%
  select(`Taxonomic level`, 
         `Number of Taxa`, 
         `Number of Classified Taxa`, 
         `Mean Relative Abundance (%) of Classified Taxa Across Samples`, 
         `Percentage of Classified Taxa`
         )
stable3
write_xlsx(stable3, 
           "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryTable3.xlsx")

#ALPHA DIVERSITY ######
##OTU LEVEL#####
alpha_div1 <- phyloseq::estimate_richness(data.bacteria.samples, measures = c("Observed", "Shannon")) # richness, diversity
alpha_div2 <- microbiome::evenness(data.bacteria.samples, index = "pielou", 
                                   zeroes = FALSE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
                                   detection = 0) ##evenness


# combine metrics with metadata
alpha_div <- cbind(alpha_div1, alpha_div2)
alpha_div
alpha_div_meta <- cbind(data.bacteria.samples@sam_data, alpha_div) 
alpha_div_meta # metadata and div metrics

#Factor variables
alpha_div_meta$feedlot <- factor(alpha_div_meta$feedlot,
                                 levels = c("1", "2", "3", "4", "5"))
alpha_div_meta$sample_type <- factor(alpha_div_meta$sample_type,
                                         levels = c("Feces", "Water"))
alpha_div_meta$gen_material <- factor(alpha_div_meta$gen_material,
                                          levels = c("DNA", "cDNA"))

#Pivot to long format 
alpha_div_meta_long <- 
  alpha_div_meta %>%
  pivot_longer(cols = c(Observed, Shannon, pielou),  
               names_to = "alpha_div_metric", 
               values_to = "alpha_div_value") 
alpha_div_meta_long
##Factoring alpha div metrics
alpha_div_meta_long$alpha_div_metric<- factor(alpha_div_meta_long$alpha_div_metric, levels = c("Observed","pielou", "Shannon"))

##LM MODEL ALPHA DIV INDECES#####
###Observed (Richness)#######
#How does its distribution look?
ggplot(alpha_div_meta, aes (x = Observed))+
  geom_histogram()

#Model (LM)
# do the linear mixed model with random sampleID 
Observed_model_otu <- lmerTest::lmer(Observed~ sample_type * gen_material * feedlot + (1 | original_sample),
                                       data = alpha_div_meta, 
                                       REML = TRUE)# for better estimate of the random-effects variance


#Check model
plot(Observed_model_otu) 
qqnorm(residuals(Observed_model_otu))
qqline(residuals(Observed_model_otu))
shapiro.test(residuals(Observed_model_otu))  

#Summary
summary(Observed_model_otu)
#Confidence Intervals
confint(Observed_model_otu)

#Anova type3 - instead of testing each coefficient individually, Type III ANOVA tests the factor as a whole.
Anova(Observed_model_otu, type = "III") #Interaction between sample_type and feedlot, as well as among all effects (three-way)

#####Emmeans - three way interaction
emmeans(Observed_model_otu, pairwise~sample_type|gen_material + feedlot) #Feedlots 1 and 3 sig differences between feces and water in both DNA and RNA (cDNA). Feedlots 4 and 5 sig differences in cDNA between feces and water
emmeans(Observed_model_otu, pairwise~gen_material|sample_type + feedlot) #In feces across all feedlots, DNA more OTUs than RNA. For water, only in feedlots 2 and 4
emmeans(Observed_model_otu, pairwise~feedlot|gen_material + sample_type) #Feedlots 1, 2, and 3 had more OTUs than feedlot 4 in water


###Evenness#######
#How does its distribution look?
ggplot(alpha_div_meta, aes (x = pielou))+
  geom_histogram()

#Model (LM)
Evenness_model_otu <- lmerTest::lmer(pielou~ sample_type * gen_material * feedlot + (1 | original_sample),
                                       data = alpha_div_meta, 
                                       REML = TRUE)# for better estimate of the random-effects variance
#Check model
plot(Evenness_model_otu) 
qqnorm(residuals(Evenness_model_otu))
qqline(residuals(Evenness_model_otu))
shapiro.test(residuals(Evenness_model_otu)) #Not really normally distributed 

#Summary
summary(Evenness_model_otu)
#Confidence Intervals
confint(Evenness_model_otu)
#####Anova type3
Anova(Evenness_model_otu, type = "III")

#####Emmeans
emmeans(Evenness_model_otu, pairwise~sample_type|gen_material+feedlot) ##Sample type : gen_material : feedlot interaction 
emmeans(Evenness_model_otu, pairwise~gen_material|sample_type+feedlot) ##Sample type : gen_material: feedlot interaction 
emmeans(Evenness_model_otu, pairwise~feedlot|gen_material+sample_type) ##Sample type : gen_material : feedlot interaction 

###Shannon#######
#How does its distribution look?
ggplot(alpha_div_meta, aes (x = Shannon))+
  geom_histogram()
alpha_div_meta
#Model (MIXED EFFECTS LINEAR REGRESSION)
Shannons_model_otu <- lmerTest::lmer(Shannon ~ sample_type * gen_material * feedlot + (1 | original_sample),
                                       data = alpha_div_meta,
                                       REML = TRUE) # for better estimate of the random-effects variances

#Check model
plot(Shannons_model_otu) 
qqnorm(residuals(Shannons_model_otu))
qqline(residuals(Shannons_model_otu))
shapiro.test(residuals(Shannons_model_otu)) #Not really normally distributed 

#Summary
summary(Shannons_model_otu) 
#Confidence Intervals
confint(Shannons_model_otu)

#Anova type3 - instead of testing each coefficient individually, Type III ANOVA tests the factor as a whole#
Anova(Shannons_model_otu, type = "III") #Interaction sample_type:gen_material , and sample_type:gen_material:feedlot
#Emmeans
emmeans(Shannons_model_otu, pairwise~sample_type|gen_material+feedlot) ##Sample type : gen_material : feedlot interaction 
emmeans(Shannons_model_otu, pairwise~gen_material|sample_type+feedlot) ##Sample type : gen_material: feedlot interaction 
emmeans(Shannons_model_otu, pairwise~feedlot|gen_material+sample_type) ##Sample type : gen_material : feedlot interaction 

###ANOVA TABLES TOGETHER ######
####RICHNESS#####
richness_anovaIII <- data.frame(
  Anova(Observed_model_otu, type = "III", test.statistic = "F"),
  check.names = F) %>%
  mutate(Metric = "Richness (Observed)",
         Dataset = "Microbiome (Bacteria - Archaea)",
         `F` = round(`F`, 2),
         Df = round(Df, 2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3)
  )%>%
  rownames_to_column(var = "Fixed Effects")%>% 
  filter(`Fixed Effects` != "(Intercept)") %>%#Don't need intercept
  #Renaming strings
  mutate(`Fixed Effects` = str_replace(`Fixed Effects`, "sample_type", "Sample Type"))%>%
  mutate(`Fixed Effects` = str_replace(`Fixed Effects`, "feedlot", "Feedlot"))%>%
  mutate(`Fixed Effects` = str_replace(`Fixed Effects`, "gen_material", "Library Type"))


####EVENNESS#####
evenness_anovaIII <- data.frame(
  Anova(Evenness_model_otu, type = "III", test.statistic = "F"),
  check.names = F)%>%
  mutate(Metric = "Evenness (pielou)",
         Dataset = "Microbiome (Bacteria - Archaea)",
         `F` = round(`F`, 2),
         Df = round(Df, 2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3)
  )%>%
  rownames_to_column(var = "Fixed Effects")%>% 
  filter(`Fixed Effects` != "(Intercept)") %>%#Don't need intercept
  #Renaming strings
  mutate(`Fixed Effects` = str_replace(`Fixed Effects`, "sample_type", "Sample Type"))%>%
  mutate(`Fixed Effects` = str_replace(`Fixed Effects`, "feedlot", "Feedlot"))%>%
  mutate(`Fixed Effects` = str_replace(`Fixed Effects`, "gen_material", "Library Type"))

####SHANNON#####
shannons_anovaIII <- data.frame(
  Anova(Shannons_model_otu, type = "III", test.statistic = "F"),
  check.names = F)%>%
  mutate(Metric = "Shannon",
         Dataset = "Microbiome (Bacteria - Archaea)",
         `F` = round(`F`, 2),
         Df = round(Df, 2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3)
         )%>%
  rownames_to_column(var = "Fixed Effects")%>% 
  filter(`Fixed Effects` != "(Intercept)") %>%#Don't need intercept
  #Renaming strings
  mutate(`Fixed Effects` = str_replace(`Fixed Effects`, "sample_type", "Sample Type"))%>%
  mutate(`Fixed Effects` = str_replace(`Fixed Effects`, "feedlot", "Feedlot"))%>%
  mutate(`Fixed Effects` = str_replace(`Fixed Effects`, "gen_material", "Library Type"))

alpha_div_anovaIII_microbiome <- bind_rows(richness_anovaIII,
                                          evenness_anovaIII,
                                          shannons_anovaIII)


#####SUPPLEMENTARY TABLE 6_1######
stable6.1 <- alpha_div_anovaIII_microbiome%>%
  select(Dataset, Metric, `Fixed Effects`, `F`, Df, `Pr(>F)`)
write_xlsx(stable6.1, "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryTable6_1.xlsx")

###EMMEANS TABLES#######
####RICHNESS####
#####Pairwise~sample_type|gen_material + feedlot
richness_emmeans_st_gm_feedlot <- emmeans(Observed_model_otu, pairwise~sample_type|gen_material + feedlot)
richness_emmeans_st_gm_feedlot_df <- data.frame(richness_emmeans_st_gm_feedlot$contrasts, check.names = F)%>%
  mutate(Metric = "Richness (Observed)",
         Dataset = "Microbiome (Bacteria - Archaea)")

#####Pairwise~gen_material|sample_type + feedlot
richness_emmeans_gm_st_feedlot <- emmeans(Observed_model_otu, pairwise~gen_material|sample_type + feedlot) 
richness_emmeans_gm_st_feedlot_df <- data.frame(richness_emmeans_gm_st_feedlot$contrasts, check.names = F)%>%
  mutate(Metric = "Richness (Observed)",
         Dataset = "Microbiome (Bacteria - Archaea)")

#####Pairwise~feedlot|gen_material + sample_type
richness_emmeans_feedlot_gm_st <-emmeans(Observed_model_otu, pairwise~feedlot|gen_material + sample_type)
richness_emmeans_feedlot_gm_st_df <- data.frame(richness_emmeans_feedlot_gm_st$contrasts, check.names = F)%>%
  mutate(Metric = "Richness (Observed)",
         Dataset = "Microbiome (Bacteria - Archaea)")

####EVENNESS####
#####Pairwise~sample_type|gen_material + feedlot
evenness_emmeans_st_gm_feedlot <- emmeans(Evenness_model_otu, pairwise~sample_type|gen_material + feedlot)
evenness_emmeans_st_gm_feedlot_df <- data.frame(evenness_emmeans_st_gm_feedlot$contrasts, check.names = F)%>%
  mutate(Metric = "Evenness (pielou)",
         Dataset = "Microbiome (Bacteria - Archaea)")

#####Pairwise~gen_material|sample_type + feedlot
evenness_emmeans_gm_st_feedlot <- emmeans(Evenness_model_otu, pairwise~gen_material|sample_type + feedlot) 
evenness_emmeans_gm_st_feedlot_df <- data.frame(evenness_emmeans_gm_st_feedlot$contrasts, check.names = F)%>%
  mutate(Metric = "Evenness (pielou)",
         Dataset = "Microbiome (Bacteria - Archaea)")

#####Pairwise~feedlot|gen_material + sample_type
evenness_emmeans_feedlot_gm_st <-emmeans(Evenness_model_otu, pairwise~feedlot|gen_material + sample_type)
evenness_emmeans_feedlot_gm_st_df <- data.frame(evenness_emmeans_feedlot_gm_st$contrasts, check.names = F)%>%
  mutate(Metric = "Evenness (pielou)",
         Dataset = "Microbiome (Bacteria - Archaea)")


####SHANNONS#####
#####Pairwise~sample_type|gen_material + feedlot
shannon_emmeans_st_gm_feedlot <- emmeans(Shannons_model_otu, pairwise~sample_type|gen_material + feedlot)
shannon_emmeans_st_gm_feedlot_df <- data.frame(shannon_emmeans_st_gm_feedlot$contrasts, check.names = F)%>%
  mutate(Metric = "Shannon",
         Dataset = "Microbiome (Bacteria - Archaea)")

#####Pairwise~gen_material|sample_type + feedlot
shannon_emmeans_gm_st_feedlot <- emmeans(Shannons_model_otu, pairwise~gen_material|sample_type + feedlot) 
shannon_emmeans_gm_st_feedlot_df <- data.frame(shannon_emmeans_gm_st_feedlot$contrasts, check.names = F)%>%
  mutate(Metric = "Shannon",
         Dataset = "Microbiome (Bacteria - Archaea)")

#####Pairwise~feedlot|gen_material + sample_type
shannon_emmeans_feedlot_gm_st <-emmeans(Shannons_model_otu, pairwise~feedlot|gen_material + sample_type)
shannon_emmeans_feedlot_gm_st_df <- data.frame(shannon_emmeans_feedlot_gm_st$contrasts, check.names = F)%>%
  mutate(Metric = "Shannon",
         Dataset = "Microbiome (Bacteria - Archaea)")

#Bind all together
alpha_div_emmeans_microbiome <- bind_rows(richness_emmeans_st_gm_feedlot_df,
                                          richness_emmeans_gm_st_feedlot_df,
                                          richness_emmeans_feedlot_gm_st_df,
                                          evenness_emmeans_st_gm_feedlot_df,
                                          evenness_emmeans_gm_st_feedlot_df,
                                          evenness_emmeans_feedlot_gm_st_df,
                                          shannon_emmeans_st_gm_feedlot_df,
                                          shannon_emmeans_gm_st_feedlot_df,
                                          shannon_emmeans_feedlot_gm_st_df)
alpha_div_emmeans_microbiome

#Editing
alpha_div_emmeans_microbiome <- alpha_div_emmeans_microbiome%>%
  select(Dataset, contrast, gen_material, sample_type, feedlot, everything())%>% #Order columns
  rename(Contrast = contrast, 
         `Library Type`= gen_material,
         `Sample Type` = sample_type,
         Feedlot = feedlot,
         )%>%
  mutate(across(where(is.factor), as.character)) %>%  # Convert factors to character
  mutate(Estimate = round(estimate, 2), 
         SE = round(SE, 3), 
         Df = round(df, 2), 
         `T ratio` = round(`t.ratio`, 2), 
         `P value` = format(p.value, scientific = TRUE, digits = 3))%>%
  mutate(across(where(is.character),
                ~ str_replace_all(., c("Water" = "Catch Basins",
                                       "cDNA" = "Metatranscriptomic (RNA (cDNA))", 
                                       "(?<!c)DNA" = "Metagenomic (DNA)", 
                                       "feedlot" = "Feedlot"))))%>% #"Match DNA only if it is NOT immediately preceded by c."
  select(Dataset, Metric, `Sample Type`, `Library Type`, Feedlot, Contrast, Estimate, SE, Df, `T ratio`, `P value`)
alpha_div_emmeans_microbiome              

#####SUPPLEMENTARY TABLE 7_1 ###### 
stable7.1 <- alpha_div_emmeans_microbiome
write_xlsx(stable7.1, 
           "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryTable7_1.xlsx")


###EMMEANS PLOTS######
###Plotting data####
Shannon_otu_WvF.cDNAandDNA_emmeans <- emmip(Shannons_model_otu,~sample_type|gen_material,
                                              CIs = T, 
                                              type = "response",
                                              nesting.order = F, 
                                              plotit = F)%>%
  mutate(alpha_div_metric = "Shannon")

Evenness_otu_WvF.cDNAandDNA_emmean <- emmip(Evenness_model_otu,~sample_type|gen_material,
                                              CIs = T, 
                                              type = "response",
                                              nesting.order = F, 
                                              plotit = F)%>%
  mutate(alpha_div_metric = "pielou")
Observed_otu_WvF.cDNAandDNA_emmean <- emmip(Observed_model_otu,~sample_type|gen_material,
                                              CIs = T, 
                                              type = "response",
                                              nesting.order = F, 
                                              plotit = F)%>%
  mutate(alpha_div_metric = "Observed")
##Put them together 
alpha_div_emmeans_data <- bind_rows(Shannon_otu_WvF.cDNAandDNA_emmeans, 
                                    Evenness_otu_WvF.cDNAandDNA_emmean, 
                                    Observed_otu_WvF.cDNAandDNA_emmean)


###Plotting data - facetted by gen material and feedlot####
Shannon_otu_WvF.cDNAandDNA_feedlot_emmeans <- emmip(Shannons_model_otu,~sample_type|gen_material + feedlot,
                                            CIs = T, 
                                            type = "response",
                                            nesting.order = F, 
                                            plotit = F)%>%
  mutate(alpha_div_metric = "Shannon")

Evenness_otu_WvF.cDNAandDNA_feedlot_emmean <- emmip(Evenness_model_otu,~sample_type|gen_material + feedlot,
                                            CIs = T, 
                                            type = "response",
                                            nesting.order = F, 
                                            plotit = F)%>%
  mutate(alpha_div_metric = "pielou")
Observed_otu_WvF.cDNAandDNA_feedlot_emmean <- emmip(Observed_model_otu,~sample_type|gen_material + feedlot,
                                            CIs = T, 
                                            type = "response",
                                            nesting.order = F, 
                                            plotit = F)%>%
  mutate(alpha_div_metric = "Observed")
##Put them together 
alpha_div_emmeans_data_WvsF_cDNAandDNA_feedlot <- bind_rows(Shannon_otu_WvF.cDNAandDNA_feedlot_emmeans, 
                                                            Evenness_otu_WvF.cDNAandDNA_feedlot_emmean, 
                                                            Observed_otu_WvF.cDNAandDNA_feedlot_emmean)

####PLOTS#######
#####Feces vs Water facet by cDNA and DNA#####
alpha_div_emmeans_WvsF_DNAandcDNA <- alpha_div_emmeans_data %>%
  ggplot(aes(x = sample_type, y = yvar, color = sample_type)) +
  geom_jitter(aes(x = sample_type,
                  y = alpha_div_value,
                  shape = feedlot),
              alpha = 0.35,
              width = 0.2,
              size = 3,
              data = alpha_div_meta_long) +#raw data
  geom_point(size = 4, shape = 20) + ##emmean
  geom_errorbar(aes(ymin = LCL, ymax = UCL), 
                position = position_dodge(width = 0.5), 
                width = 0.2,
                linewidth = 1.5) + #error bars for confidence intervals
  scale_color_manual(values = sample.type.palette, 
                     name = "Sample Type",
                     labels = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  scale_shape_manual(name = "Feedlot",
                     values = c(15,17,18,19, 12))+
  facet_grid(alpha_div_metric~ gen_material,
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "DNA" = "DNA", 
                                      "cDNA" = "RNA\n(cDNA)")), 
             scales = "free")+
  theme_bw() +
  labs(title= "MICROBIOME") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22, 
                                face = "bold"),
    legend.box = "vertical",   # stack legend boxes vertically
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text = element_text(colour = "white", size = 28, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 30, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  geom_pwc(method = "emmeans_test",
           label = "p.signif",
           step.increase = 0.1,
           size = 0.5,
           label.size = 5,
           tip.length = 0.02,
           hide.ns = T) 
alpha_div_emmeans_WvsF_DNAandcDNA

#####Feces vs CB facet by cDNA and DNA as well as feedlot#####
alpha_div_emmeans_WvsF_DNAandcDNA_feedlot <- alpha_div_emmeans_data_WvsF_cDNAandDNA_feedlot %>%
  ggplot(aes(x = sample_type, y = yvar, color = sample_type)) +
  geom_jitter(aes(x = sample_type,
                  y = alpha_div_value),
              alpha = 0.35,
              width = 0.2,
              size = 3,
              data = alpha_div_meta_long) +#raw data
  geom_point(size = 4, shape = 20) + ##emmean
  geom_errorbar(aes(ymin = LCL, ymax = UCL), 
                position = position_dodge(width = 0.5), 
                width = 0.2,
                linewidth = 1.5) + #error bars for confidence intervals
  scale_color_manual(values = sample.type.palette, 
                     name = "Sample Type",
                     labels = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  facet_grid(alpha_div_metric ~ feedlot + gen_material ,
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "DNA" = "DNA", 
                                      "cDNA" = "RNA\n(cDNA)", 
                                      "1" = "1",
                                      "2" = "2", 
                                      "3" = "3", 
                                      "4" = "4", 
                                      "5" = "5")), 
             scales = "free")+
  theme_bw() +
  labs(title= "MICROBIOME") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22, 
                                face = "bold"),
    legend.box = "vertical",   # stack legend boxes vertically
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text = element_text(colour = "white", size = 28, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 30, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  geom_pwc(method = "emmeans_test",
           label = "p.signif",
           step.increase = 0.1,
           size = 0.5,
           label.size = 5,
           tip.length = 0.02,
           hide.ns = T) 
alpha_div_emmeans_WvsF_DNAandcDNA_feedlot

######FIGURE 3A - FOCUS ON SAMPLE TYPE#####
figure3A <- plot_grid(alpha_div_emmeans_WvsF_DNAandcDNA_feedlot+
                        theme(plot.title = element_blank()), 
                      labels = c("A"), 
                      label_size = 22)
figure3A
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure3A.png", 
       plot = figure3A, 
       device = "png",
       dpi = 600,
       # width = 8,
       # height =12,
       width = 16,
       height = 10,
       bg = "white")

#####DNA vs cDNA facet by Feces and Water#####
alpha_div_emmeans_cDNAvsDNA_WandF <- alpha_div_emmeans_data %>%
  ggplot(aes(x = gen_material, y = yvar, color = gen_material)) +
  geom_jitter(aes(x = gen_material,
                  y = alpha_div_value,
                  shape = feedlot),
              alpha = 0.35,
              width = 0.2,
              size = 3,
              data = alpha_div_meta_long) +#raw data
  geom_point(size = 4, shape = 20) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), 
                position = position_dodge(width = 0.5), 
                width = 0.2,
                linewidth = 1.5) +
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  scale_shape_manual(name = "Feedlot",
                     values = c(15,17,18,19, 12))+
  theme_bw() +
  labs(
    #title= "MICROBIOME", 
    color = "Library Type", 
    fill = "Library Type") +
  facet_grid(alpha_div_metric~sample_type,
             scales = "free",
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "Water" = "CB", 
                                      "Feces" = "FECES"))) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        legend.box = "vertical",   # stack legend boxes vertically
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  geom_pwc(method = "emmeans_test",
           label = "{p.signif}",
           step.increase = 0.1,
           size = 0.5,
           label.size = 5,
           tip.length = 0.02,
           # hide.ns = T
           ) 
alpha_div_emmeans_cDNAvsDNA_WandF

#####DNA vs cDNA facet by Feces and CB and feedlot#####
alpha_div_emmeans_cDNAvsDNA_WandF_feedlot <- 
  alpha_div_emmeans_data_WvsF_cDNAandDNA_feedlot %>%
  ggplot(aes(x = gen_material, y = yvar, color = gen_material)) +
  geom_jitter(aes(x = gen_material,
                  y = alpha_div_value),
              alpha = 0.5,
              width = 0.2,
              size = 3,
              data = alpha_div_meta_long) +#raw data
  geom_point(size = 4, shape = 20) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), 
                position = position_dodge(width = 0.5), 
                width = 0.2,
                linewidth = 1.5) +
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  theme_bw() +
  labs(title= "MICROBIOME", color = "Library Type") +
  facet_grid(alpha_div_metric ~ feedlot +sample_type ,
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "Feces" = "FECES", 
                                      "Water" = "CB", 
                                      "1" = "1",
                                      "2" = "2", 
                                      "3" = "3", 
                                      "4" = "4", 
                                      "5" = "5")), 
             scales = "free")+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        legend.box = "vertical",   # stack legend boxes vertically
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  geom_pwc(method = "emmeans_test",
           label = "p.signif",
           step.increase = 0.1,
           size = 0.5,
           label.size = 5,
           tip.length = 0.02,
           hide.ns = T) 
alpha_div_emmeans_cDNAvsDNA_WandF_feedlot


######FIGURE 6A - FOCUS ON LIBRARY TYPE#####
figure6A <- plot_grid(alpha_div_emmeans_cDNAvsDNA_WandF_feedlot+
                        theme(plot.title = element_blank()), 
                      labels = c("A"), 
                      label_size = 22)
figure6A
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure6A.png", 
       plot = figure6A, 
       device = "png", 
       dpi = 600,
       # width = 8, 
       # height =12,
       width = 16,
       height = 10,
       bg = "white")


## BOX PLOTS ALPHA DIV ######
###Feces vs Water faceted by cDNA and DNA#####
alpha_div_WvF.cDNAandDNA <- ggplot(alpha_div_meta_long, aes(x = sample_type, y= alpha_div_value, fill= sample_type, colour = sample_type)) +
  theme_bw() +
  labs(title= "ALPHA DIVERSITY", color = "Sample Type", fill = "Sample Type") +
  facet_grid(alpha_div_metric~gen_material,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "DNA" = "DNA", 
                                      "cDNA" = "RNA (cDNA)"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
  scale_color_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  scale_fill_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  #For easy plotting of p values (I know it's not the same as ANOVA III, but keeps the significant results)
  # stat_compare_means(method = "anova", 
  #                    label.y.npc = "top", 
  #                    hide.ns = TRUE,  
  #                    show.legend = F,
  #                    label = "p.signif",
  #                    size = 8)
  geom_pwc (method = "wilcox_test",
            label = "Wilcoxon, p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
alpha_div_WvF.cDNAandDNA 

###cDNA vs DNA faceted by Water and Feces#####
alpha_div_cDNAvsDNA.WandF <- ggplot(alpha_div_meta_long, aes(x = gen_material, y= alpha_div_value, fill= gen_material, colour = gen_material)) +
  theme_bw() +
  labs(title= "ALPHA DIVERSITY", color = "Library Type", fill = "Library Type") +
  facet_grid(alpha_div_metric~sample_type,
             scales = "free",
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "Water" = "CATCH BASINS", 
                                      "Feces" = "FECES"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
  scale_fill_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  # stat_compare_means(method = "anova",
  #                    label.y.npc = "top",
  #                    hide.ns = TRUE,
  #                    show.legend = F,
  #                    label = "p.signif",
  #                    size = 8)
  geom_pwc (method = "wilcox_test",
            label = "p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
alpha_div_cDNAvsDNA.WandF

#####Rare taxa?#######
##PONDS
data.bacteria.samples.tss.water <- subset_samples(data.bacteria.samples.tss, sample_type == "Water")
data.bacteria.samples.tss.water <- prune_taxa(taxa_sums(data.bacteria.samples.tss.water)>0, data.bacteria.samples.tss.water)
data.bacteria.samples.tss.water #29583 otus and 24 samples

###Ponds DNA
data.bacteria.samples.tss.water.DNA <- subset_samples(data.bacteria.samples.tss.water, gen_material == "DNA")
data.bacteria.samples.tss.water.DNA <- prune_taxa(taxa_sums(data.bacteria.samples.tss.water.DNA)>0, data.bacteria.samples.tss.water.DNA)
data.bacteria.samples.tss.water.DNA #27975 OTUs and 12 samples

data.bacteria.samples.tss.water.DNA.melt <- psmelt(data.bacteria.samples.tss.water.DNA)
low_abund_water_DNA <- data.bacteria.samples.tss.water.DNA.melt%>%
  group_by(OTU) %>%
  summarise(mean_abun = mean(Abundance), .groups = "drop") %>%
  filter(mean_abun < 0.15 & mean_abun > 0) 
nrow(low_abund_water_DNA) #27909 low abun (<0.15% RA) OTUs in water DNA 
nrow(low_abund_water_DNA)/nrow(get_taxa(data.bacteria.samples.tss.water.DNA)) # 0.997 (99.7%) of pond DNA OTUs in average are low abun (<0.15% RA)

###Ponds cDNA
data.bacteria.samples.tss.water.cDNA <- subset_samples(data.bacteria.samples.tss.water, gen_material == "cDNA")
data.bacteria.samples.tss.water.cDNA <- prune_taxa(taxa_sums(data.bacteria.samples.tss.water.cDNA)>0, data.bacteria.samples.tss.water.cDNA)
data.bacteria.samples.tss.water.cDNA #26124 OTUs and 12 samples

data.bacteria.samples.tss.water.cDNA.melt <- psmelt(data.bacteria.samples.tss.water.cDNA)
low_abund_water_cDNA <- data.bacteria.samples.tss.water.cDNA.melt%>%
  group_by(OTU) %>%
  summarise(mean_abun = mean(Abundance), 
            .groups = "drop") %>%
  filter(mean_abun < 0.15 & mean_abun > 0) 
nrow(low_abund_water_cDNA) #26048 low abun OTUs in cDNA 
nrow(low_abund_water_cDNA)/nrow(get_taxa(data.bacteria.samples.tss.water.cDNA)) # 0.997 (99.7%) of pond cDNA OTUs in average are low abun (<0.15% RA)


##FECES
data.bacteria.samples.tss.feces <- subset_samples(data.bacteria.samples.tss, sample_type == "Feces")
data.bacteria.samples.tss.feces <- prune_taxa(taxa_sums(data.bacteria.samples.tss.feces)>0, data.bacteria.samples.tss.feces)
data.bacteria.samples.tss.feces #33493 otus and 96 samples 

###Feces DNA
data.bacteria.samples.tss.feces.DNA <- subset_samples(data.bacteria.samples.tss.feces, gen_material == "DNA")
data.bacteria.samples.tss.feces.DNA <- prune_taxa(taxa_sums(data.bacteria.samples.tss.feces.DNA)>0, data.bacteria.samples.tss.feces.DNA)
data.bacteria.samples.tss.feces.DNA #31904 taxa and 48 samples 

data.bacteria.samples.tss.feces.DNA.melt <- psmelt(data.bacteria.samples.tss.feces.DNA)
low_abund_feces_DNA <- data.bacteria.samples.tss.feces.DNA.melt%>%
  group_by(OTU) %>%
  summarise(mean_abun = mean(Abundance), .groups = "drop") %>%
  filter(mean_abun < 0.15 & mean_abun > 0) 
nrow(low_abund_feces_DNA) #31816 low abun OTUs in DNA 
nrow(low_abund_feces_DNA)/nrow(get_taxa(data.bacteria.samples.tss.feces.DNA)) # 0.997 (99.7%) of fecal DNA OTUs in average are low abun

###Feces cDNA
data.bacteria.samples.tss.feces.cDNA <- subset_samples(data.bacteria.samples.tss.feces, gen_material == "cDNA")
data.bacteria.samples.tss.feces.cDNA <- prune_taxa(taxa_sums(data.bacteria.samples.tss.feces.cDNA)>0, data.bacteria.samples.tss.feces.cDNA)
data.bacteria.samples.tss.feces.cDNA #28722 OTUs and 48 samples

data.bacteria.samples.tss.feces.cDNA.melt <- psmelt(data.bacteria.samples.tss.feces.cDNA)
low_abund_feces_cDNA <- data.bacteria.samples.tss.feces.cDNA.melt%>%
  group_by(OTU) %>%
  summarise(mean_abun = mean(Abundance), .groups = "drop") %>%
  filter(mean_abun < 0.3 & mean_abun > 0) 
nrow(low_abund_feces_cDNA) #28682 low abun OTUs in cDNA 
nrow(low_abund_feces_cDNA)/nrow(get_taxa(data.bacteria.samples.tss.feces.cDNA)) # 0.998 (99.8%) of fecal cDNA OTUs in average are low abun


###Reseq vs not reseq (batch effect check)#####
alpha_div_reseq <- ggplot(alpha_div_meta_long, aes(x = re_sequenced, y= alpha_div_value, fill= re_sequenced, colour = re_sequenced)) +
  theme_bw() +
  labs(title= "ALPHA DIVERSITY") +
  facet_grid(alpha_div_metric ~ gen_material,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "cDNA" = "cDNA", 
                                      "DNA" = "DNA"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
  scale_color_manual(values = batch_palette) +
  scale_fill_manual(values = batch_palette) +
  labs(color = "Re-Sequenced", fill = "Re-Sequenced")+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  geom_pwc (method = "wilcox_test",
            label = "p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
alpha_div_reseq

# BETA DIVERSITY ###########
##BRAY CURTIS####
###ONLY cDNA (METATRANSCRIPTOMIC) SAMPLES#############
##Subsetting only cDNA samples
data.bacteria.samples.cDNA.tss <- subset_samples(data.bacteria.samples.tss, gen_material=="cDNA")
data.bacteria.samples.cDNA.tss <- prune_taxa(taxa_sums(data.bacteria.samples.cDNA.tss) > 0, data.bacteria.samples.cDNA.tss) 

##Distance matrix
data.bacteria.samples.cDNA.tss.bray <- vegdist(t(data.bacteria.samples.cDNA.tss@otu_table), method = "bray") 
data.bacteria.samples.cDNA.tss.df <- as(data.bacteria.samples.cDNA.tss@sam_data,"data.frame") %>% # make DF from metadata
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")))

####PERMANOVA########
#Model will include sample type and feedlot
set.seed(87)
cDNA_BC_adonis_sampletype_feedlot  <- adonis2(data.bacteria.samples.cDNA.tss.bray ~ sample_type+feedlot,
                                       by = "margin",
                                       data.bacteria.samples.cDNA.tss.df, permutations = 9999)
cDNA_BC_adonis_sampletype_feedlot #Sample type and feedlot significant

######SUPPLEMENTARY TABLE 5.6 #######
##Using the model with sample type and feedlot
stable5.6 <- data.frame(cDNA_BC_adonis_sampletype_feedlot, check.names = F)%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)",
         `Library Type` = "Metatranscriptomic (RNA (cDNA))",
         SumOfSqs = round(SumOfSqs, 2), 
         `R2` = round(`R2`,3), 
         `F` = round(`F`,2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3)
         )%>%
  rownames_to_column(var = "Fixed Effect")%>%
  mutate(`Fixed Effect` = str_replace(`Fixed Effect`, "sample_type", "Sample Type"),
         `Fixed Effect` = str_replace(`Fixed Effect`, "feedlot", "Feedlot"))%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))
stable5.6

#PERMDISP - Sample type
# Run the betadisper function, average distance to centroid
bray.cDNA.sampletype.disp <- betadisper(data.bacteria.samples.cDNA.tss.bray, 
                                        data.bacteria.samples.cDNA.tss.df$sample_type)
bray.cDNA.sampletype.disp
##Then test by permuting
set.seed(87)
bray.cDNA.sampletype.permdisp <- permutest(bray.cDNA.sampletype.disp, 
                                           permutations = 9999)

bray.cDNA.sampletype.permdisp 
##Feces vs Water p-value 1e-04 (different variances)

#PERMDISP - Feedlot
# Run the betadisper function, average distance to centroid
bray.cDNA.feedlot.disp <- betadisper(data.bacteria.samples.cDNA.tss.bray, data.bacteria.samples.cDNA.tss.df$feedlot)
bray.cDNA.feedlot.disp
##Then test by permuting
set.seed(87)
bray.cDNA.feedlot.permdisp <- permutest(bray.cDNA.feedlot.disp, permutations = 9999)
bray.cDNA.feedlot.permdisp 
##Not significant for feedlots (p = 0.91)


#### ORDINATION
set.seed(87)
data.bacteria.samples.cDNA.tss.bray.ord <- metaMDS(data.bacteria.samples.cDNA.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

##Get ordination coordinates
data.bacteria.samples.cDNA.tss.bray.plot <- ordiplot(data.bacteria.samples.cDNA.tss.bray.ord$points)
#Extracts the scores (coordinates) of the points from the ordination plot object:
data.bacteria.samples.cDNA.tss.bray.scrs <- scores(data.bacteria.samples.cDNA.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata:
data.bacteria.samples.cDNA.tss.bray.scrs <- cbind(as.data.frame(data.bacteria.samples.cDNA.tss.bray.scrs), 
                                  sample_type = data.bacteria.samples.cDNA.tss.df$sample_type, 
                                  feedlot = factor(data.bacteria.samples.cDNA.tss.df$feedlot)) 

####SAMPLE TYPE EFFECT#######
## BC
data.bacteria.samples.cDNA.tss.bray.cent.sample_type <- aggregate(cbind(MDS1,MDS2) ~ sample_type, 
                                                  data = data.bacteria.samples.cDNA.tss.bray.scrs, FUN = mean) ##Centroids according to sample type (water and feces)
data.bacteria.samples.cDNA.tss.bray.segs.sample_type <- merge(data.bacteria.samples.cDNA.tss.bray.scrs, 
                                              setNames(data.bacteria.samples.cDNA.tss.bray.cent.sample_type, 
                                                       c("sample_type", "cMDS1","cMDS2")), by = 'sample_type', sort = F) ##add centroids to main scrs dataframe

##Abbreviated version of sample_type for easier plotting
data.bacteria.samples.cDNA.tss.bray.segs.sample_type  <- data.bacteria.samples.cDNA.tss.bray.segs.sample_type  %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"))

# Extract R2 and p-values - sample type
R2_cDNA_adonis_sample_type <- cDNA_BC_adonis_sampletype_feedlot$R2[1] 
pvalue_cDNA_adonis_sample_type<-  cDNA_BC_adonis_sampletype_feedlot$`Pr(>F)`[1]

#####PLOT #####
cDNA_BC_beta_div_spider_sampletype <- ggplot(data.bacteria.samples.cDNA.tss.bray.segs.sample_type) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "Metatranscriptomic libraries (RNA (cDNA))", 
       color = "Sample Type") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = sample_type), show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = sample_type), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = sample_type), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2,label= sample_type.abbrv), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values = sample.type.palette,
                     labels = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", size = 22, face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  annotate("text", x = 0.3, y = 0.8, ##change coordinates as needed
           label = "Sample type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = 0.3, y = 0.8, ##change coordinates as needed
           label = paste0("R² = ", round(R2_cDNA_adonis_sample_type * 100, 1), "%",
                          "\np = ", round(pvalue_cDNA_adonis_sample_type, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") +# Annotate R² and p-values 
  guides(
    color= guide_legend(override.aes = list(size = 7))
  )
cDNA_BC_beta_div_spider_sampletype


####FEEDLOT EFFECT#######
data.bacteria.samples.cDNA.tss.bray.scrs #Already have the ordination coordinates with metadata
data.bacteria.samples.cDNA.tss.bray.cent.feedlot <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = data.bacteria.samples.cDNA.tss.bray.scrs, 
                                              FUN = mean) ##Centroids according to feedlot
data.bacteria.samples.cDNA.tss.bray.segs.feedlot <- merge(data.bacteria.samples.cDNA.tss.bray.scrs, 
                                          setNames(data.bacteria.samples.cDNA.tss.bray.cent.feedlot, c("feedlot", "cMDS1","cMDS2")), 
                                          by = 'feedlot', 
                                          sort = F) ##add centroids to main scrs dataframe

# Extract R2 and p-values
R2_cDNA_adonis_feedlot <- cDNA_BC_adonis_sampletype_feedlot$R2[2] 
pvalue_cDNA_adonis_feedlot<-  cDNA_BC_adonis_sampletype_feedlot$`Pr(>F)`[2]

#####PLOT#####
cDNA_BC_beta_div_spider_feedlot <- ggplot(data.bacteria.samples.cDNA.tss.bray.segs.feedlot) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "Metatranscriptomic libraries (RNA (cDNA))", 
       color = "Feedlot") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = feedlot), show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = feedlot), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = feedlot), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2,label= feedlot), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values = feedlot_palette) + 
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", size = 22, face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  annotate("text", x = -0.4, y = 0.9, ##change coordinates as needed
           label = "Feedlot",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = -0.4, y = 0.9, ##change coordinates as needed
           label = paste0("R² = ", round(R2_cDNA_adonis_feedlot * 100, 1), "%",
                          "\np = ", round(pvalue_cDNA_adonis_feedlot, 3)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") +# Annotate R² and p-values
  guides(
    color= guide_legend(override.aes = list(size = 7))
  )
cDNA_BC_beta_div_spider_feedlot

####Sequenced vs resequenced#############
##Only cDNA samples
data.bacteria.samples.cDNA.tss 
##Distance matrix
data.bacteria.samples.cDNA.tss.bray 
##Metadata df
data.bacteria.samples.cDNA.tss.df 
## Ordinated object
data.bacteria.samples.cDNA.tss.bray.ord

#### ADDING CENTROIDS FOR PLOTTING 
## BC
data.bacteria.samples.cDNA.tss.bray.plot.reseq <- ordiplot(data.bacteria.samples.cDNA.tss.bray.ord$points)
data.bacteria.samples.cDNA.tss.bray.reseq.scrs <- scores(data.bacteria.samples.cDNA.tss.bray.plot.reseq, display = "sites") 
data.bacteria.samples.cDNA.tss.bray.reseq.scrs <- cbind(as.data.frame(data.bacteria.samples.cDNA.tss.bray.reseq.scrs), 
                                         sample_type = data.bacteria.samples.cDNA.tss.df$sample_type, 
                                         feedlot = data.bacteria.samples.cDNA.tss.df$feedlot,
                                        re_sequenced = data.bacteria.samples.cDNA.tss.df$re_sequenced)
data.bacteria.samples.cDNA.tss.bray.reseq.cent <- aggregate(cbind(MDS1,MDS2) ~ sample_type + re_sequenced, data = data.bacteria.samples.cDNA.tss.bray.reseq.scrs, FUN = mean) ##Centroids according to sample type (water and feces) and reseq (yes or no)
data.bacteria.samples.cDNA.tss.bray.reseq.segs <- merge(
  data.bacteria.samples.cDNA.tss.bray.reseq.scrs,
  setNames(data.bacteria.samples.cDNA.tss.bray.reseq.cent, c("sample_type", "re_sequenced", "cMDS1", "cMDS2")),
  by = c('sample_type', 're_sequenced'),
  sort = FALSE
)
data.bacteria.samples.cDNA.tss.bray.reseq.segs <- data.bacteria.samples.cDNA.tss.bray.reseq.segs%>%
  mutate(reseq.abbrv = if_else(re_sequenced == "yes", "Y", "N"))

#PERMANOVA
set.seed(87)
cDNA_BC_adonis_reseq  <- adonis2(data.bacteria.samples.cDNA.tss.bray ~ sample_type + re_sequenced, 
                                 data.bacteria.samples.cDNA.tss.df, by = "margin", permutations = 9999)
cDNA_BC_adonis_reseq  #1.2% of variation is due to resequencing (yes vs no) p = 0.1733


# Extract R2 and p-values
R2_cDNA_reseq_adonis_sample_type <- cDNA_BC_adonis_reseq$R2[2] 
pvalue_cDNA_reseq_adonis_sample_type<-  cDNA_BC_adonis_reseq$`Pr(>F)`[2]

##PLOT
cDNA_BC_beta_div_spider_reseq <- ggplot(data.bacteria.samples.cDNA.tss.bray.reseq.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", 
       title= "Metabolically active members\nwithin communities (cDNA)", 
       color = "Sample type",
       shape = "Re-sequencing") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = data.bacteria.samples.cDNA.tss.bray.reseq.segs%>% filter(sample_type=="Feces" & re_sequenced== "yes"),
               aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = sample_type), show.legend = F)+
  geom_segment(data = data.bacteria.samples.cDNA.tss.bray.reseq.segs%>% filter(sample_type=="Feces" & re_sequenced== "no"),
               aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = sample_type), show.legend = F)+
  geom_segment(data = data.bacteria.samples.cDNA.tss.bray.reseq.segs%>% filter(sample_type=="Water" & re_sequenced== "yes"),
               aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = sample_type), show.legend = F)+
  geom_segment(data = data.bacteria.samples.cDNA.tss.bray.reseq.segs%>% filter(sample_type=="Water" & re_sequenced== "no"),
               aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = sample_type), show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = sample_type, shape = re_sequenced), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = sample_type, shape = re_sequenced), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2,label= reseq.abbrv), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values = sample.type.palette) + 
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", size = 24, face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  annotate("text", x = 1.65, y = 0.9, ##change coordinates as needed
           label = "Re-sequencing", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = 1.65, y = 0.9, ##change coordinates as needed
           label = paste0("R² = ", round(R2_cDNA_reseq_adonis_sample_type * 100, 1), "%",
                          "\np = ", round(pvalue_cDNA_reseq_adonis_sample_type, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") +# Annotate R² and p-values 
  guides(
    color= guide_legend(override.aes = list(size = 7))
  )
cDNA_BC_beta_div_spider_reseq

#PERMDISP
# Run the betadisper function, average distance to centroid
sampletype_reseq_interaction <- interaction(data.bacteria.samples.cDNA.tss.df$sample_type, data.bacteria.samples.cDNA.tss.df$re_sequenced)
bray.cDNA.sampletype.disp.reseq <- betadisper(data.bacteria.samples.cDNA.tss.bray, sampletype_reseq_interaction)
bray.cDNA.sampletype.disp.reseq
##Then test by permuting
set.seed(87)
bray.cDNA.sampletype.permdisp.reseq <- permutest(bray.cDNA.sampletype.disp.reseq, permutations = 9999, pairwise = 1)
bray.cDNA.sampletype.permdisp.reseq
##Water reseq vs not reseq p =0.8933 (n.s.)
##Feces reseq vs not reseq p=0.5 (n.s.)


###ONLY DNA (METAGENOMIC) SAMPLES#######
##Subsetting only DNA samples
data.bacteria.samples.DNA.tss <- subset_samples(data.bacteria.samples.tss, gen_material=="DNA")
data.bacteria.samples.DNA.tss <- prune_taxa(taxa_sums(data.bacteria.samples.DNA.tss) > 0, data.bacteria.samples.DNA.tss) 

##Distance matrix
data.bacteria.samples.DNA.tss.bray <- vegdist(t(data.bacteria.samples.DNA.tss@otu_table), method = "bray") 
data.bacteria.samples.DNA.tss.df <- as(data.bacteria.samples.DNA.tss@sam_data,"data.frame") %>% # make DF from metadata
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")))

####PERMANOVA########
#Modelling sample type and feedlot
set.seed(87)
DNA_BC_adonis_sampletype_feedlot  <- adonis2(data.bacteria.samples.DNA.tss.bray ~ sample_type+feedlot,
                                       by = "margin",
                                       data.bacteria.samples.DNA.tss.df, permutations = 9999)
DNA_BC_adonis_sampletype_feedlot #Yes

######SUPPLEMENTARY TABLE 5.5#######
stable5.5 <- data.frame(DNA_BC_adonis_sampletype_feedlot, check.names = F)%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)",
         `Library Type` = "Metagenomic (DNA)",
         SumOfSqs = round(SumOfSqs, 2), 
         `R2` = round(`R2`,3), 
         `F` = round(`F`,2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3)
         )%>%
  rownames_to_column(var = "Fixed Effect")%>%
  mutate(`Fixed Effect` = str_replace(`Fixed Effect`, "sample_type", "Sample Type"), 
         `Fixed Effect` = str_replace(`Fixed Effect`, "feedlot", "Feedlot"))%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))
stable5.5

#PERMDISP - SAMPLE TYPE
# Run the betadisper function, average distance to centroid
bray.DNA.sampletype.disp <- betadisper(data.bacteria.samples.DNA.tss.bray, data.bacteria.samples.DNA.tss.df$sample_type)
bray.DNA.sampletype.disp
##Then test by permuting
set.seed(87)
bray.DNA.sampletype.permdisp <- permutest(bray.DNA.sampletype.disp, permutations = 9999)
bray.DNA.sampletype.permdisp 
##Different variances Feces vs Water p-value 1e-04

#PERMDISP - FEEDLOT
# Run the betadisper function, average distance to centroid
bray.DNA.feedlot.disp <- betadisper(data.bacteria.samples.DNA.tss.bray, data.bacteria.samples.DNA.tss.df$feedlot)

bray.DNA.feedlot.disp
##Then test by permuting
set.seed(87)
bray.DNA.feedlot.permdisp <- permutest(bray.DNA.feedlot.disp, permutations = 9999)
bray.DNA.feedlot.permdisp 
##No different variances between feedlots (p = 0.9)

#### ORDINATION
set.seed(87)
data.bacteria.samples.DNA.tss.bray.ord <- metaMDS(data.bacteria.samples.DNA.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

##Get ordination coordinates
data.bacteria.samples.DNA.tss.bray.plot <- ordiplot(data.bacteria.samples.DNA.tss.bray.ord$points)
#Extracts the scores (coordinates) of the points from the ordination plot object:
data.bacteria.samples.DNA.tss.bray.scrs <- scores(data.bacteria.samples.DNA.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata:
data.bacteria.samples.DNA.tss.bray.scrs <- cbind(as.data.frame(data.bacteria.samples.DNA.tss.bray.scrs), 
                                  sample_type = data.bacteria.samples.DNA.tss.df$sample_type, 
                                  feedlot = factor(data.bacteria.samples.DNA.tss.df$feedlot)) 

####SAMPLE TYPE EFFECT#######
## BC
data.bacteria.samples.DNA.tss.bray.cent.sample_type <- aggregate(cbind(MDS1,MDS2) ~ sample_type, 
                                                 data = data.bacteria.samples.DNA.tss.bray.scrs, FUN = mean) ##Centroids according to sample type (water and feces)
data.bacteria.samples.DNA.tss.bray.segs.sample_type <- merge(data.bacteria.samples.DNA.tss.bray.scrs, 
                                             setNames(data.bacteria.samples.DNA.tss.bray.cent.sample_type, 
                                                                               c("sample_type", "cMDS1","cMDS2")), 
                                             by = 'sample_type', sort = F) ##add centroids to main scrs dataframe

##Abbreviated version of sample_type for easier plotting
data.bacteria.samples.DNA.tss.bray.segs.sample_type  <- data.bacteria.samples.DNA.tss.bray.segs.sample_type  %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"))

# Extract R2 and p-values
R2_DNA_adonis_sample_type <- DNA_BC_adonis_sampletype_feedlot$R2[1] 
pvalue_DNA_adonis_sample_type<-  DNA_BC_adonis_sampletype_feedlot$`Pr(>F)`[1]

#####PLOT#####
DNA_BC_beta_div_spider_sampletype <- ggplot(data.bacteria.samples.DNA.tss.bray.segs.sample_type) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", 
       title= "Metagenomic libraries (DNA)", 
       color = "Sample Type") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = sample_type), show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = sample_type), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = sample_type), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2,label= sample_type.abbrv), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values = sample.type.palette,
                     labels = c("Feces"= "Feces", "Water" = "Catch Basins"))+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", size = 22, face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"))+
  annotate("text", x = 0.3, y = 0.4, ##change coordinates as needed
           label = "Sample type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = 0.3, y = 0.4, ##change coordinates as needed
           label = paste0("R² = ", round(R2_DNA_adonis_sample_type * 100, 1), "%",
                          "\np = ", round(pvalue_DNA_adonis_sample_type, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") +# Annotate R² and p-values 
  guides(
    color= guide_legend(override.aes = list(size = 7))
  )
DNA_BC_beta_div_spider_sampletype

####FEEDLOT EFFECT#######
data.bacteria.samples.DNA.tss.bray.scrs #Already have the ordination coordinates with metadata
data.bacteria.samples.DNA.tss.bray.cent.feedlot <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = data.bacteria.samples.DNA.tss.bray.scrs, 
                                              FUN = mean) ##Centroids according to feedlot
data.bacteria.samples.DNA.tss.bray.segs.feedlot <- merge(data.bacteria.samples.DNA.tss.bray.scrs, 
                                          setNames(data.bacteria.samples.DNA.tss.bray.cent.feedlot, c("feedlot", "cMDS1","cMDS2")), 
                                          by = 'feedlot', 
                                          sort = F) ##add centroids to main scrs dataframe

# Extract R2 and p-values
R2_DNA_adonis_feedlot <- DNA_BC_adonis_sampletype_feedlot$R2[2] 
pvalue_DNA_adonis_feedlot<-  DNA_BC_adonis_sampletype_feedlot$`Pr(>F)`[2]

#####PLOT#####
DNA_BC_beta_div_spider_feedlot <- ggplot(data.bacteria.samples.DNA.tss.bray.segs.feedlot) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", 
       title= "Metagenomic libraries (DNA)",
       color = "Feedlot") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = feedlot), show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = feedlot), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = feedlot), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2,label= feedlot), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values = feedlot_palette) + 
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", size = 22, face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  annotate("text", x = -0.4, y = 0.4, ##change coordinates as needed
           label = "Feedlot",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = -0.4, y = 0.4, ##change coordinates as needed
           label = paste0("R² = ", round(R2_DNA_adonis_feedlot * 100, 1), "%",
                          "\np = ", round(pvalue_DNA_adonis_feedlot, 4)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") +# Annotate R² and p-values
  guides(
    color= guide_legend(override.aes = list(size = 7))
  )
DNA_BC_beta_div_spider_feedlot

######SAMPLE TYPE COMPARISON PLOTS - DNA AND CDNA SAMPLES######
DNAandcDNA_BC_beta_div_spider_sampletype <- plot_grid(DNA_BC_beta_div_spider_sampletype+ 
                                                 theme(axis.title.x = element_blank(), legend.position = "none"), 
                                                 cDNA_BC_beta_div_spider_sampletype + 
                                                   theme(axis.title = element_blank(), legend.position = "none"),
                                               labels = "AUTO", label_size = 22, 
                                               rel_widths = c(1.05, 1))
DNAandcDNA_BC_beta_div_spider_sampletype

##Extracting legend and x axis
leg_WaterandFeces_fill <- get_plot_component(DNA_BC_beta_div_spider_sampletype, "guide-box", return_all = TRUE)
leg_WaterandFeces_fill <- ggdraw(leg_WaterandFeces_fill[[3]])##extract the legend, make into ggplot object
x_axis_label_NMDS1 <- ggdraw() + draw_label("NMDS1", size = 28, hjust = 0.5) ##Make an xaxis title 


##Just want one legend for the main plot and one x axis
DNAandcDNA_BC_beta_div_spider_sampletype <- plot_grid(
  DNAandcDNA_BC_beta_div_spider_sampletype, # Main plot
  x_axis_label_NMDS1,# Common x-axis label
  leg_WaterandFeces_fill, # Legend
  ncol = 1, 
  rel_heights = c(1, 0.05, 0.1) # Adjust heights as needed
)
DNAandcDNA_BC_beta_div_spider_sampletype


######FEEDLOT COMPARISON PLOTS - DNA AND CDNA SAMPLES######
DNAandcDNA_BC_beta_div_spider_feedlot <- plot_grid(DNA_BC_beta_div_spider_feedlot+ 
                                                        theme(axis.title.x = element_blank(), legend.position = "none"), 
                                                      cDNA_BC_beta_div_spider_feedlot + 
                                                        theme(axis.title = element_blank(), legend.position = "none"),
                                                      labels = "AUTO", label_size = 22, 
                                                      rel_widths = c(1.05, 1))
DNAandcDNA_BC_beta_div_spider_feedlot

##Extracting legend and x axis
leg_WaterandFeces_fill <- get_plot_component(DNA_BC_beta_div_spider_feedlot, "guide-box", return_all = TRUE)
leg_WaterandFeces_fill <- ggdraw(leg_WaterandFeces_fill[[3]])##extract the legend, make into ggplot object
x_axis_label_NMDS1 <- ggdraw() + draw_label("NMDS1", size = 28, hjust = 0.5) ##Make an xaxis title 

##Just want one legend for the main plot and one x axis
DNAandcDNA_BC_beta_div_spider_feedlot <- plot_grid(
  DNAandcDNA_BC_beta_div_spider_feedlot, # Main plot
  x_axis_label_NMDS1,# Common x-axis label
  leg_WaterandFeces_fill, # Legend
  ncol = 1, 
  rel_heights = c(1, 0.05, 0.1) # Adjust heights as needed
)
DNAandcDNA_BC_beta_div_spider_feedlot


###CATCH BASIN SAMPLES#######
##Subsetting only catch basin (wastewater) samples
data.bacteria.samples.water.tss <- subset_samples(data.bacteria.samples.tss, sample_type=="Water")
data.bacteria.samples.water.tss <- prune_taxa(taxa_sums(data.bacteria.samples.water.tss) > 0, data.bacteria.samples.water.tss) 

##Distance matrix
data.bacteria.samples.water.tss.bray <- vegdist(t(data.bacteria.samples.water.tss@otu_table), method = "bray") 
data.bacteria.samples.water.tss.df <- as(data.bacteria.samples.water.tss@sam_data,"data.frame") # make DF from metadata

##Adding a column concatenating "gen_material" and "feedlot", then making feedlot a factor
data.bacteria.samples.water.tss.df<- data.bacteria.samples.water.tss.df %>%
  mutate (gen_material_spec_2 = paste(gen_material, feedlot, sep = '_'),
          feedlot = factor(feedlot))

####PERMANOVA#####
#Is there an interaction betwwen library type and feedlot? 
set.seed(87)
water_BC_adonis_interaction <- adonis2(data.bacteria.samples.water.tss.bray ~ gen_material * feedlot,
                                       by = "margin",
                                       data.bacteria.samples.water.tss.df, 
                                       p.adjust.methods = "BH", permutations = 9999)
water_BC_adonis_interaction #No

##Model Library Type and feedlot, stratified by original sample from which DNA and RNA were extracted 
set.seed(87)
water_BC_adonis <- adonis2(data.bacteria.samples.water.tss.bray ~ gen_material + feedlot, 
                           strata = data.bacteria.samples.water.tss.df$original_sample,
                           by = "margin",data.bacteria.samples.water.tss.df, 
                           p.adjust.methods = "BH", permutations = 9999)
water_BC_adonis #4.1% of the variation is due to Library Type, p = 0.0005
#62.9% of the variation is due to feedlot, p = 0.0005

######SUPPLEMENTARY TABLE 5.8 #######
stable5.8 <- data.frame(water_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)",
         `Sample Type` = "Catch Basins",
         SumOfSqs = round(SumOfSqs, 2), 
         `R2` = round(`R2`,3), 
         `F` = round(`F`,2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3)
  )%>%
  rownames_to_column(var = "Fixed Effect")%>%
  mutate(`Fixed Effect` = str_replace(`Fixed Effect`, "gen_material", "Library Type"), 
         `Fixed Effect` = str_replace(`Fixed Effect`, "feedlot", "Feedlot"))%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))
stable5.8

#PERMDISP- Library Type
# Run the betadisper function, average distance to centroid
bray.water.genmat.disp <- betadisper(data.bacteria.samples.water.tss.bray, data.bacteria.samples.water.tss.df$gen_material)
bray.water.genmat.disp
##Then test by permuting
set.seed(87)
bray.water.genmat.permdisp <- permutest(bray.water.genmat.disp, permutations = 9999)
bray.water.genmat.permdisp  #No different variances

#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.water.feedlot.disp <- betadisper(data.bacteria.samples.water.tss.bray, 
                                      data.bacteria.samples.water.tss.df$feedlot)
bray.water.feedlot.disp
##Then test by permuting
set.seed(87)
bray.water.feedlot.permdisp <- permutest(bray.water.feedlot.disp, permutations = 9999, pairwise = 1)
bray.water.feedlot.permdisp #Different variances p = 0.02

#### ORDINATION
set.seed(87)
data.bacteria.samples.water.tss.bray.ord <- metaMDS(data.bacteria.samples.water.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

#Extract ordination coordinates, add metadata
data.bacteria.samples.water.tss.bray.plot <- ordiplot(data.bacteria.samples.water.tss.bray.ord$points)
#Extract the scores (coordinates) of the points from the ordination plot object 
data.bacteria.samples.water.tss.bray.scrs <- scores(data.bacteria.samples.water.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
data.bacteria.samples.water.tss.bray.scrs <- cbind(as.data.frame(data.bacteria.samples.water.tss.bray.scrs), 
                                   gen_material = data.bacteria.samples.water.tss.df$gen_material, 
                                   feedlot = factor(data.bacteria.samples.water.tss.df$feedlot), 
                                   gen_material.recode = ifelse(data.bacteria.samples.water.tss.df$gen_material == "DNA", "DNA", "RNA"),
                                   gen_material_spec_2 = data.bacteria.samples.water.tss.df$gen_material_spec_2,
                                   sampleID = rownames(data.bacteria.samples.water.tss.df))

####FEEDLOT EFFECT#####
## BC
data.bacteria.samples.water.tss.bray.scrs #have coodrinates and metadata
##Centroids according to feedlot
data.bacteria.samples.water.feedlot.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = data.bacteria.samples.water.tss.bray.scrs, FUN = mean) 
data.bacteria.samples.water.feedlot.tss.bray.segs <- merge(data.bacteria.samples.water.tss.bray.scrs, 
                                           setNames(data.bacteria.samples.water.feedlot.tss.bray.cent, c("feedlot","cMDS1","cMDS2")),
                                           by = 'feedlot', sort = F)

# Extract R2 and p-values
##Model gen_material and feedlot
water_BC_adonis #have the model
R2_water_feedlot_BC_adonis <- water_BC_adonis$R2[2] 
pvalue_water_feedlot_BC_adonis<-water_BC_adonis$`Pr(>F)`[2]

#### PLOT
water_feedlot_BC_beta_div <- ggplot(data.bacteria.samples.water.feedlot.tss.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "CATCH BASINS", color = "Feedlot") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  # stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= color_gen_material, colour = color_gen_material), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = feedlot), show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = feedlot), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = feedlot), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2, label= feedlot), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values = feedlot_palette) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22, face = "bold"),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  annotate("text", x = 0.54, y = 0.8, ##change coordinates as needed
           label = "Feedlot", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (feedlot)
  annotate("text", x = 0.54, y = 0.8, ##change coordinates as needed
           label = paste("R² = ", round(R2_water_feedlot_BC_adonis * 100, 1), "%",
                         "\np = ", round(pvalue_water_feedlot_BC_adonis, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black")+ # Annotate R² and p-values
  guides(color = guide_legend(override.aes = list(size = 7)))
water_feedlot_BC_beta_div 

####LIBRARY TYPE EFFECT#####
## BC
data.bacteria.samples.water.tss.bray.scrs #have coordinates and metadata
##Centroids according to gen_material (cDNA/DNA)
data.bacteria.samples.water.genmat.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ gen_material, 
                                              data = data.bacteria.samples.water.tss.bray.scrs, FUN = mean) 
data.bacteria.samples.water.genmat.tss.bray.segs <- merge(data.bacteria.samples.water.tss.bray.scrs, 
                                          setNames(data.bacteria.samples.water.genmat.tss.bray.cent, c("gen_material","cMDS1","cMDS2")), by = 'gen_material', sort = F)

# Extract R2 and p-values for genmat
water_BC_adonis #have the model
R2_water_BC_adonis_genmat <- water_BC_adonis$R2[1] 
pvalue_water_BC_adonis_genmat<-water_BC_adonis$`Pr(>F)`[1]

#### PLOT
water_genmat_BC_beta_div <- ggplot(data.bacteria.samples.water.genmat.tss.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "CATCH BASINS",
       color = "Library Type") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  # stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= color_gen_material, colour = color_gen_material), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = gen_material), show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = gen_material), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = gen_material), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2,label= gen_material.recode), colour= "white", size = 3) +
  scale_color_manual(values = gen.material.palette, 
                     labels = c("cDNA" = "RNA (cDNA)", "DNA" = "DNA")) +
  theme(legend.position = "bottom",
        legend.title= element_text(colour = "black", size = 22, face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 16, colour = "black")) +
  guides(
    color = guide_legend (override.aes = list(size = 7))
  )+
  annotate("text", x = 0.6, y = 0.7, ##change coordinates as needed
           label = "Library Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Library Type)
  annotate("text", x = 0.6, y = 0.7, ##change coordinates as needed
           label = paste("R² = ", round(R2_water_BC_adonis_genmat * 100, 1), "%",
                         "\np = ", round(pvalue_water_BC_adonis_genmat, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
water_genmat_BC_beta_div 

###FECAL SAMPLES#######
##Subsetting only feces samples
data.bacteria.samples.feces.tss <- subset_samples(data.bacteria.samples.tss, sample_type=="Feces")
data.bacteria.samples.feces.tss <- prune_taxa(taxa_sums(data.bacteria.samples.feces.tss) > 0, data.bacteria.samples.feces.tss) 

##Distance matrix
data.bacteria.samples.feces.tss.bray <- vegdist(t(data.bacteria.samples.feces.tss@otu_table), method = "bray") 
data.bacteria.samples.feces.tss.df <- as(data.bacteria.samples.feces.tss@sam_data,"data.frame") # make DF from metadata

##Adding a column concatenating "gen_material" and "feedlot", then making feedlot a factor
data.bacteria.samples.feces.tss.df<- data.bacteria.samples.feces.tss.df %>%
  mutate (gen_material_spec_2 = paste(gen_material, feedlot, sep = '_'),
          feedlot = factor(feedlot))

####PERMANOVA#####
#Is there an interaction between library type and feedlot? 
set.seed(87)
feces_BC_adonis_interaction <- adonis2(data.bacteria.samples.feces.tss.bray ~ gen_material * feedlot, 
                                       by = "margin",
                                       data.bacteria.samples.feces.tss.df, 
                                       p.adjust.methods = "BH", permutations = 9999)
feces_BC_adonis_interaction #No

##Model gen_material and feedlot, stratify by original sample from which DNA and RNA were extracted from 
set.seed(87)
feces_BC_adonis <- adonis2(data.bacteria.samples.feces.tss.bray ~ gen_material + feedlot,
                           strata = data.bacteria.samples.feces.tss.df$original_sample,
                           by = "margin",
                           data.bacteria.samples.feces.tss.df, 
                           p.adjust.methods = "BH", permutations = 9999)
feces_BC_adonis #36.8% if variation explained by gen_material, p = 1e-04
#12.9% of the varaition explained by feedlot, p = 1e-04

######SUPPLEMENTARY TABLE 5.7 #######
stable5.7 <- data.frame(feces_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)",
         `Sample Type` = "Feces",
         SumOfSqs = round(SumOfSqs, 2), 
         `R2` = round(`R2`,3), 
         `F` = round(`F`,2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3)
  )%>%
  rownames_to_column(var = "Fixed Effect")%>%
  mutate(`Fixed Effect` = str_replace(`Fixed Effect`, "gen_material", "Library Type"), 
         `Fixed Effect` = str_replace(`Fixed Effect`, "feedlot", "Feedlot"))%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))
stable5.7

#PERMDISP- Gen_material
# Run the betadisper function, average distance to centroid
bray.feces.genmat.disp <- betadisper(data.bacteria.samples.feces.tss.bray, data.bacteria.samples.feces.tss.df$gen_material)
bray.feces.genmat.disp
##Then test by permuting
set.seed(87)
bray.feces.genmat.permdisp <- permutest(bray.feces.genmat.disp, permutations = 9999)
bray.feces.genmat.permdisp #No different variances p = 0.108

#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.feces.feedlot.disp <- betadisper(data.bacteria.samples.feces.tss.bray, data.bacteria.samples.feces.tss.df$feedlot)
bray.feces.feedlot.disp
##Then test by permuting
set.seed(87)
bray.feces.feedlot.permdisp <- permutest(bray.feces.feedlot.disp, permutations = 9999, pairwise = 1)
bray.feces.feedlot.permdisp #Different variances  p = 0.0043

#### ORDINATION
set.seed(87)
data.bacteria.samples.feces.tss.bray.ord <- metaMDS(data.bacteria.samples.feces.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

#Extract ordination coordinates, add metadata
data.bacteria.samples.feces.tss.bray.plot <- ordiplot(data.bacteria.samples.feces.tss.bray.ord$points)
#Extract the scores (coordinates) of the points from the ordination plot object 
data.bacteria.samples.feces.tss.bray.scrs <- scores(data.bacteria.samples.feces.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
data.bacteria.samples.feces.tss.bray.scrs <- cbind(as.data.frame(data.bacteria.samples.feces.tss.bray.scrs), 
                                   gen_material = data.bacteria.samples.feces.tss.df$gen_material, 
                                   gen_material.recode = ifelse(data.bacteria.samples.feces.tss.df$gen_material == "DNA", "DNA", "RNA"),
                                   feedlot = factor(data.bacteria.samples.feces.tss.df$feedlot), 
                                   gen_material_spec_2 = data.bacteria.samples.feces.tss.df$gen_material_spec_2,
                                   sampleID = rownames(data.bacteria.samples.feces.tss.df))

####FEEDLOT EFFECT #####
## BC
data.bacteria.samples.feces.tss.bray.scrs #have coodrinates and metadata
##Centroids according to feedlot
data.bacteria.samples.feces.feedlot.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = data.bacteria.samples.feces.tss.bray.scrs, FUN = mean) 
data.bacteria.samples.feces.feedlot.tss.bray.segs <- merge(data.bacteria.samples.feces.tss.bray.scrs, 
                                           setNames(data.bacteria.samples.feces.feedlot.tss.bray.cent, c("feedlot","cMDS1","cMDS2")),
                                           by = 'feedlot', sort = F)

# Extract R2 and p-values
##Model gen_material and feedlot
feces_BC_adonis #have the model
R2_feces_feedlot_BC_adonis <- feces_BC_adonis$R2[2] 
pvalue_feces_feedlot_BC_adonis<-feces_BC_adonis$`Pr(>F)`[2]

#### PLOT
feces_feedlot_BC_beta_div <- ggplot(data.bacteria.samples.feces.feedlot.tss.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "FECES", color = "Feedlot") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  # stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= color_gen_material, colour = color_gen_material), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = feedlot), show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = feedlot), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = feedlot), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2, label= feedlot), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values = feedlot_palette) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22, face = "bold"),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  annotate("text", x = -0.5, y = 0.4, ##change coordinates as needed
           label = "Feedlot", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (feedlot)
  annotate("text", x = -0.5, y = 0.4, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces_feedlot_BC_adonis * 100, 1), "%",
                         "\np = ", round(pvalue_feces_feedlot_BC_adonis, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black")+ # Annotate R² and p-values
  guides(color = guide_legend(override.aes = list(size = 7)))
feces_feedlot_BC_beta_div 

####LIBRARY TYPE EFFECT#####
## BC
data.bacteria.samples.feces.tss.bray.scrs #have coordinates and metadata
##Centroids according to gen_material (cDNA/DNA)
data.bacteria.samples.feces.genmat.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ gen_material, 
                                              data = data.bacteria.samples.feces.tss.bray.scrs, FUN = mean) 
data.bacteria.samples.feces.genmat.tss.bray.segs <- merge(data.bacteria.samples.feces.tss.bray.scrs, 
                                          setNames(data.bacteria.samples.feces.genmat.tss.bray.cent, c("gen_material","cMDS1","cMDS2")), 
                                          by = 'gen_material', sort = F)

# Extract R2 and p-values for genmat
feces_BC_adonis #have the model
R2_feces_BC_adonis_genmat <- feces_BC_adonis$R2[1] 
pvalue_feces_BC_adonis_genmat<-feces_BC_adonis$`Pr(>F)`[1]

#### PLOT
feces_genmat_BC_beta_div <- ggplot(data.bacteria.samples.feces.genmat.tss.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "FECES",
       color = "Library Type") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  # stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= color_gen_material, colour = color_gen_material), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = gen_material), show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = gen_material), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = gen_material), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2, label= gen_material.recode), colour= "white", size = 3) +
  scale_color_manual(values = gen.material.palette, 
                     labels = c("cDNA" = "RNA (cDNA)", "DNA" = "DNA")) +
  theme(legend.position = "bottom",
        legend.title= element_text(colour = "black", size = 22, face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 16, colour = "black")) +
  guides(
    color = guide_legend (override.aes = list(size = 7))
  )+
  annotate("text", x = 0.4, y = -0.05, ##change coordinates as needed
           label = "Library Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Library Type)
  annotate("text", x = 0.4, y = -0.05, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces_BC_adonis_genmat * 100, 1), "%",
                         "\np = ", round(pvalue_feces_BC_adonis_genmat, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
feces_genmat_BC_beta_div 

###### PLOT OF LIBRARY EFFECT IN FECES AND CATCH BASIN SAMPLES#######
WaterandFeces_genmat_BC_beta_div <- plot_grid(feces_genmat_BC_beta_div+theme(axis.title.x = element_blank(), legend.position = "none"), 
                                              water_genmat_BC_beta_div+theme(axis.title = element_blank(), legend.position = "none"), 
                                              labels = "AUTO", label_size = 22,
                                              rel_widths = c(1, 0.8))
WaterandFeces_genmat_BC_beta_div
##Extracting legend and x axis
leg_WaterandFeces_genmat <- get_plot_component(feces_genmat_BC_beta_div, "guide-box", return_all = TRUE)
leg_WaterandFeces_genmat <- ggdraw(leg_WaterandFeces_genmat[[3]])##extract the legend, make into ggplot object


x_axis_label_NMDS1 <- ggdraw() + draw_label("NMDS1", size = 28, hjust = 0.5) ##Make an xaxis title 

##Just want one legend for the main plot and one x axis
WaterandFeces_genmat_BC_beta_div  <- plot_grid(
  WaterandFeces_genmat_BC_beta_div , # Main plot
  x_axis_label_NMDS1,# Common x-axis label
  leg_WaterandFeces_genmat, # Legend
  ncol = 1, 
  rel_heights = c(1, 0.05, 0.1) # Adjust heights as needed
)
WaterandFeces_genmat_BC_beta_div 

###### FEEDLOT EFFECT IN FECES AND CATCH BASINS#######
WaterandFeces_feedlot_BC_beta_div <- plot_grid(feces_feedlot_BC_beta_div + theme(axis.title.x = element_blank(), legend.position = "none"), 
                                               water_feedlot_BC_beta_div +theme(axis.title = element_blank(), legend.position = "none"), 
                                               labels = c("A", "B"), label_size = 22,
                                               rel_widths = c(1, 0.9))
WaterandFeces_feedlot_BC_beta_div
##Extracting legend and x axis
leg_WaterandFeces_feedlot <- get_plot_component(water_feedlot_BC_beta_div, "guide-box", return_all = TRUE)
leg_WaterandFeces_feedlot <- ggdraw(leg_WaterandFeces_feedlot[[3]])##extract the legend, make into ggplot object


x_axis_label_NMDS1 <- ggdraw() + draw_label("NMDS1", size = 28, hjust = 0.5) ##Make an xaxis title 

##Just want one legend for the main plot and one x axis
WaterandFeces_feedlot_BC_beta_div <- plot_grid(
  WaterandFeces_feedlot_BC_beta_div, # Main plot
  x_axis_label_NMDS1,# Common x-axis label
  leg_WaterandFeces_feedlot, # Legend
  ncol = 1, 
  rel_heights = c(1, 0.05, 0.1) # Adjust heights as needed
)
WaterandFeces_feedlot_BC_beta_div


####SUPPLEMENTARY FIGURE 6AB#######
#####Effect of feedlot in DNA and cDNA samples#####
sfigure6AandB <- ggarrange(DNA_BC_beta_div_spider_feedlot+
                             theme(plot.title = element_blank()), 
                           cDNA_BC_beta_div_spider_feedlot+
                             theme(plot.title = element_blank()),
                           labels = "AUTO", 
                           font.label = list(size = 22),
                           legend = "none",
                           common.legend = T)
  # labs(title = "MICROBIOME")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
sfigure6AandB
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryFigure6AB.svg", 
       plot = sfigure6AandB, device = "svg", 
       #dpi = 600, 
       bg = "white",
       width = 18, height = 6)

####SUPPLEMENTARY FIGURE 7AB#######
#####Effect of feedlot On fecal and catch basin samples#####
sfigure7AandB <- ggarrange(feces_feedlot_BC_beta_div+
                             theme(plot.title = element_blank()), 
                           water_feedlot_BC_beta_div+
                             theme(plot.title = element_blank()),
                           labels = c("A", "B"), 
                           font.label = list(size = 22),
                           legend = "none",
                           common.legend = T)
  # labs(title = "MICROBIOME")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
sfigure7AandB
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryFigure7AB.svg", 
       plot = sfigure7AandB, 
       device = "svg", 
       width = 18,
       height = 6, 
       dpi = 600)

###SALMONELLA CULTURE STATUS #######
####FECES#########
data.bacteria.samples.feces.tss ##already have this object with normalized counts
data.bacteria.samples.feces.tss.bray ##Distance matrix
data.bacteria.samples.feces.tss.df #metadata dataframe
data.bacteria.samples.feces.tss.bray.ord ##ordination

data.bacteria.samples.feces.salmonella.tss.df<- data.bacteria.samples.feces.tss.df %>%
  mutate (feedlot_salmonella = paste(feedlot, salmonella_culture_status, sep = '_'))%>% ##adding a column concatenating "feedlot" and "Salmonella"
  mutate (salmonella_stat.abbrv = dplyr::recode(salmonella_culture_status, "positive"= "P", "negative"= "N"),
          feedlot = factor(feedlot)) ##making feedlot a factor

#PERMANOVA
set.seed(87)
feces_salmonella_BC_adonis <- adonis2(data.bacteria.samples.feces.tss.bray ~ salmonella_culture_status + feedlot, 
                                      by = "margin", data.bacteria.samples.feces.salmonella.tss.df,  
                                      permutations = 9999)
feces_salmonella_BC_adonis 
#2.8% of variation is due to salmonella_culture_status (positive vs negative) p= 0.039 (s)
#12.2% of variation is due to feedlot p= 0.001 (S)

#PERMANOVA- Salmonella status interaction with feedlot 
set.seed(87)
feces_salmonella_feedlot_adonis_interaction <- adonis2(data.bacteria.samples.feces.tss.bray ~ salmonella_culture_status*feedlot, 
                                                       by = "margin", 
                                             data.bacteria.samples.feces.salmonella.tss.df, permutations = 9999)
feces_salmonella_feedlot_adonis_interaction
#5.3% of variation is due to salmonella_culture_status interaction with feedlot (p= 0.052) 

#PERMANOVA- Salmonella status stratified by feedlot 
set.seed(87)
feces_salmonella_BC_adonis_status <- adonis2(data.bacteria.samples.feces.tss.bray ~ salmonella_culture_status, 
                                             by = "margin", 
                                             strata = data.bacteria.samples.feces.salmonella.tss.df$feedlot, 
                                             data.bacteria.samples.feces.salmonella.tss.df, permutations = 9999)
feces_salmonella_BC_adonis_status
#2.7% of variation is due to salmonella_culture_status (positive vs negative) p= 0.0523 

#PERMANOVA- Feedlot stratified by salmonella status
set.seed(87)
feces_salmonella_BC_adonis_feedlot <- adonis2(data.bacteria.samples.feces.tss.bray ~ feedlot, 
                                              by = "margin", 
                                              strata = data.bacteria.samples.feces.salmonella.tss.df$salmonella_culture_status, 
                                              data.bacteria.samples.feces.salmonella.tss.df, p.adjust.methods = "BH", permutations = 9999)
feces_salmonella_BC_adonis_feedlot
#13% of variation is due to feedlot p= 0.0014

#####FECES DNA#######
data.bacteria.samples.feces.DNA.tss <- subset_samples(data.bacteria.samples.feces.tss, gen_material =="DNA")
data.bacteria.samples.feces.DNA.tss <- prune_taxa(taxa_sums(data.bacteria.samples.feces.DNA.tss) > 0, 
                                                  data.bacteria.samples.feces.DNA.tss) 

data.bacteria.samples.feces.DNA.tss #34167 taxa and 48 samples

#metadata 
data.bacteria.samples.feces.DNA.tss.df <- data.frame(data.bacteria.samples.feces.DNA.tss@sam_data)

###### Including only feedlots that have both positive and negative samples#######
## Taking out feedlot 5 since it did not have any positive samples###
data.bacteria.samples.feces.DNA.tss #Only feces DNA samples 
data.bacteria.samples.feces.DNA.f14.tss <- subset_samples(data.bacteria.samples.feces.DNA.tss, feedlot !="5") #Subsetting out feedlot 5 - doesn't have salmonella positive samples
data.bacteria.samples.feces.DNA.f14.tss <- prune_taxa(taxa_sums(data.bacteria.samples.feces.DNA.f14.tss) > 0, data.bacteria.samples.feces.DNA.f14.tss) 
data.bacteria.samples.feces.DNA.f14.tss ##33683 taxa and 44 samples

##Distance matrix
data.bacteria.samples.feces.DNA.f14.tss.bray <- vegdist(t(data.bacteria.samples.feces.DNA.f14.tss@otu_table), method = "bray") 
#Ordination 
set.seed(87)
data.bacteria.samples.feces.DNA.f14.tss.bray.ord <- metaMDS(data.bacteria.samples.feces.DNA.f14.tss.bray, k = 2, try= 20, trymax= 1000, autotransform = F)

##metadata
data.bacteria.samples.feces.DNA.f14.salmonella.tss.df<- data.bacteria.samples.feces.DNA.tss.df %>%
  filter(feedlot != "5")%>%
  mutate (feedlot_salmonella = paste(feedlot, salmonella_culture_status, sep = '_'))%>% ##adding a column concatenating "feedlot" and "Salmonella"
  mutate (salmonella_stat.abbrv = dplyr::recode(salmonella_culture_status, "positive"= "P", "negative"= "N"),
          feedlot = factor(feedlot))

#### ADDING CENTROIDS FOR PLOTTING
## BC- centroids on only salmonella status
data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.plot <- ordiplot(data.bacteria.samples.feces.DNA.f14.tss.bray.ord$points)
data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.scrs <- scores(data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.plot, display = "sites") 
data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.scrs <- cbind(as.data.frame(data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.scrs), 
                                                       salmonella_culture_status = data.bacteria.samples.feces.DNA.f14.salmonella.tss.df$salmonella_culture_status,
                                                       salmonella_stat.abbrv = data.bacteria.samples.feces.DNA.f14.salmonella.tss.df$salmonella_stat.abbrv ) #This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ salmonella_culture_status, data = data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.scrs, FUN = mean) ##Centroids according to gen_material_spec_2 (cDNA/DNA_feedlot)
data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.segs <- merge(data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.scrs, setNames(data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.cent, c("salmonella_culture_status","cMDS1","cMDS2")), by = 'salmonella_culture_status', sort = F)

######PERMANOVA########
set.seed(87)
feces.DNA_salmonellaonly_BC_adonis <- adonis2(data.bacteria.samples.feces.DNA.f14.tss.bray ~ salmonella_culture_status + feedlot,
                                              data.bacteria.samples.feces.DNA.f14.salmonella.tss.df,
                                              by = "margin",
                                              p.adjust.methods = "BH", 
                                              permutations = 9999)
feces.DNA_salmonellaonly_BC_adonis #5.5% of variation is due to Salmonella Status (p= 0.0486)

######SUPPLEMENTARY TABLE 5.1#######
stable5.1 <- data.frame(feces.DNA_salmonellaonly_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)",
         `Library Type` = "Metagenomic (DNA)", 
         `Sample Type` = "Feces",
         SumOfSqs = round(SumOfSqs, 2), 
         `R2` = round(`R2`,3), 
         `F` = round(`F`,2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3))%>%
  rownames_to_column(var = "Fixed Effect")%>%
  mutate(`Fixed Effect` = str_replace(`Fixed Effect`, "salmonella_culture_status", "Salmonella Culture Status"), 
         `Fixed Effect` = str_replace(`Fixed Effect`, "feedlot", "Feedlot"))%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))
stable5.1

#PERMDISP- Salmonella status
# Run the betadisper function, average distance to centroid
bray.feces.DNA.f14.salmonella.disp <- betadisper(data.bacteria.samples.feces.DNA.f14.tss.bray, data.bacteria.samples.feces.DNA.f14.salmonella.tss.df$salmonella_culture_status)
bray.feces.DNA.f14.salmonella.disp
##Then test by permuting
set.seed(87)
bray.feces.DNA.f14.salmonella.permdisp <- permutest(bray.feces.DNA.f14.salmonella.disp, permutations = 9999)
bray.feces.DNA.f14.salmonella.permdisp 
##No different dispersions of variance between salmonella status (p = 0.38)

#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.feces.DNA.f14.feedlot.disp <- betadisper(data.bacteria.samples.feces.DNA.f14.tss.bray, 
                                               data.bacteria.samples.feces.DNA.f14.salmonella.tss.df$feedlot)
bray.feces.DNA.f14.feedlot.disp
##Then test by permuting
set.seed(87)
bray.feces.DNA.f14.feedlot.permdisp <- permutest(bray.feces.DNA.f14.salmonella.disp, permutations = 9999)
bray.feces.DNA.f14.feedlot.permdisp 
##No different dispersions of variance between feedlots (p = 0.38)

# Extract R2 and p-values
#Salmonella 
R2_feces.DNA_salmonella_BC_adonis_status <- feces.DNA_salmonellaonly_BC_adonis$R2[1] 
pvalue_feces.DNA_salmonella_BC_adonis_status<-feces.DNA_salmonellaonly_BC_adonis$`Pr(>F)`[1]

##Plot
feces.DNA_salmonella_only_BC_beta_div <- ggplot(data.bacteria.samples.feces.DNA.f14.salmonella.only.tss.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "FECES (DNA)", color = "Salmonella Status") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = salmonella_culture_status))+
  geom_point(aes(x=MDS1, y=MDS2, colour = salmonella_culture_status), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = salmonella_culture_status), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2, label= salmonella_stat.abbrv), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values= salmonella.palette)+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 22),
        legend.text = element_text(colour = "black", size = 22),
        legend.title.position = "top",
        legend.title.align = 0.5,
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+ 
  guides(
    color = guide_legend(override.aes = list(size = 7)))+
  annotate("text", x = 0.2, y = -0.15, ##change coordinates as needed
           label = "Salmonella Status",
           hjust = 0.5, vjust = -0.5, size = 5, colour = "black", fontface = "bold") + ##annotate variable (SAlmonella status)
  annotate("text", x =0.2, y = -0.15, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces.DNA_salmonella_BC_adonis_status * 100, 1), "%",
                         "\np = ", round(pvalue_feces.DNA_salmonella_BC_adonis_status, 3)),
           hjust = 0.5, vjust = 1.1, size = 5, colour = "black") # Annotate R² and p-values
feces.DNA_salmonella_only_BC_beta_div

#####FECES cDNA#######
data.bacteria.samples.feces.cDNA.tss <- subset_samples(data.bacteria.samples.feces.tss, gen_material =="cDNA")
data.bacteria.samples.feces.cDNA.tss <- prune_taxa(taxa_sums(data.bacteria.samples.feces.cDNA.tss) > 0, 
                                                  data.bacteria.samples.feces.cDNA.tss) 

data.bacteria.samples.feces.cDNA.tss #29803 taxa and 48 samples

#metadata
data.bacteria.samples.feces.cDNA.tss.df <- data.frame(data.bacteria.samples.feces.cDNA.tss@sam_data)

###### Including only feedlots with both positive and negative samples##########
data.bacteria.samples.feces.cDNA.tss #Only feces cDNA samples 
data.bacteria.samples.feces.cDNA.f14.tss <- subset_samples(data.bacteria.samples.feces.cDNA.tss, feedlot !="5") #Subsetting out feedlot 5 - doesn't have salmonella positive samples
data.bacteria.samples.feces.cDNA.f14.tss <- prune_taxa(taxa_sums(data.bacteria.samples.feces.cDNA.f14.tss) > 0, data.bacteria.samples.feces.cDNA.f14.tss) 
data.bacteria.samples.feces.cDNA.f14.tss ##29437 taxa and 44 samples

##Distance matrix
data.bacteria.samples.feces.cDNA.f14.tss.bray <- vegdist(t(data.bacteria.samples.feces.cDNA.f14.tss@otu_table), method = "bray") 
#Ordination 
set.seed(87)
data.bacteria.samples.feces.cDNA.f14.tss.bray.ord <- metaMDS(data.bacteria.samples.feces.cDNA.f14.tss.bray, k = 2, try= 20, trymax= 1000, autotransform = F)

##metadata
data.bacteria.samples.feces.cDNA.f14.salmonella.tss.df<- data.bacteria.samples.feces.cDNA.tss.df %>%
  filter(feedlot != "5")%>%
  mutate (feedlot_salmonella = paste(feedlot, salmonella_culture_status, sep = '_'))%>% ##adding a column concatenating "feedlot" and "Salmonella"
  mutate (salmonella_stat.abbrv = dplyr::recode(salmonella_culture_status, "positive"= "P", "negative"= "N"),
          feedlot = factor(feedlot))

#### ADDING CENTROIDS FOR PLOTTING
## BC- centroids on only salmonella status
data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.plot <- ordiplot(data.bacteria.samples.feces.cDNA.f14.tss.bray.ord$points)
data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.scrs <- scores(data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.plot, display = "sites") 
data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.scrs <- cbind(as.data.frame(data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.scrs), 
                                                           salmonella_culture_status = data.bacteria.samples.feces.cDNA.f14.salmonella.tss.df$salmonella_culture_status,
                                                           salmonella_stat.abbrv = data.bacteria.samples.feces.cDNA.f14.salmonella.tss.df$salmonella_stat.abbrv ) #This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ salmonella_culture_status, 
                                                                data = data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.scrs, FUN = mean) ##Centroids according to gen_material_spec_2 (ccDNA/cDNA_feedlot)
data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.segs <- merge(data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.scrs, 
                                                            setNames(data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.cent, c("salmonella_culture_status","cMDS1","cMDS2")), by = 'salmonella_culture_status', sort = F)

#PERMANOVA 
set.seed(87)
feces.cDNA_salmonellaonly_BC_adonis <- adonis2(data.bacteria.samples.feces.cDNA.f14.tss.bray ~ salmonella_culture_status + feedlot, 
                                              data.bacteria.samples.feces.cDNA.f14.salmonella.tss.df,
                                              by = "margin",
                                              p.adjust.methods = "BH", 
                                              permutations = 9999)
feces.cDNA_salmonellaonly_BC_adonis #4.5% of variation is due to Salmonella Status (p= 0.0908)
#######SUPPLEMENTARY TABLE 5.2#######
stable5.2 <- data.frame(feces.cDNA_salmonellaonly_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)",
         `Library Type` = "Metatranscriptomic (RNA (cDNA))", 
         `Sample Type` = "Feces",
         SumOfSqs = round(SumOfSqs, 2), 
         `R2` = round(`R2`,3), 
         `F` = round(`F`,2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3))%>%
  rownames_to_column(var = "Fixed Effect")%>%
  mutate(`Fixed Effect` = str_replace(`Fixed Effect`, "salmonella_culture_status", "Salmonella Culture Status"), 
         `Fixed Effect` = str_replace(`Fixed Effect`, "feedlot", "Feedlot"))%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))
stable5.2

#PERMDISP- Salmonella status
# Run the betadisper function, average distance to centroid
bray.feces.cDNA.f14.salmonella.disp <- betadisper(data.bacteria.samples.feces.cDNA.f14.tss.bray, data.bacteria.samples.feces.cDNA.f14.salmonella.tss.df$salmonella_culture_status)
bray.feces.cDNA.f14.salmonella.disp
##Then test by permuting
set.seed(87)
bray.feces.cDNA.f14.salmonella.permdisp <- permutest(bray.feces.cDNA.f14.salmonella.disp, permutations = 9999)
bray.feces.cDNA.f14.salmonella.permdisp 
##No different dispersions of variance between salmonella status (p = 0.94)

#PERMDISP - Feedlot
# Run the betadisper function, average distance to centroid
bray.feces.cDNA.f14.feedlot.disp <- betadisper(data.bacteria.samples.feces.cDNA.f14.tss.bray, 
                                               data.bacteria.samples.feces.cDNA.f14.salmonella.tss.df$feedlot)
bray.feces.cDNA.f14.feedlot.disp
##Then test by permuting
set.seed(87)
bray.feces.cDNA.f14.feedlot.permdisp <- permutest(bray.feces.cDNA.f14.salmonella.disp, permutations = 9999)
bray.feces.cDNA.f14.feedlot.permdisp 
##No different dispersions of variance between feedlots (p = 0.94)


# Extract R2 and p-values
#Salmonella 
R2_feces.cDNA_salmonella_BC_adonis_status <- feces.cDNA_salmonellaonly_BC_adonis$R2[1] 
pvalue_feces.cDNA_salmonella_BC_adonis_status<-feces.cDNA_salmonellaonly_BC_adonis$`Pr(>F)`[1]

##Plot
feces.cDNA_salmonella_only_BC_beta_div <- ggplot(data.bacteria.samples.feces.cDNA.f14.salmonella.only.tss.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "FECES (RNA(cDNA))", color = "Salmonella Status") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = salmonella_culture_status))+
  geom_point(aes(x=MDS1, y=MDS2, colour = salmonella_culture_status), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = salmonella_culture_status), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2, label= salmonella_stat.abbrv), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values= salmonella.palette)+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 22),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+ 
  guides(
    color = guide_legend (override.aes = list(size = 7)))+
  annotate("text", x = 0.3, y = -0.2, ##change coordinates as needed
           label = "Salmonella Status",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (SAlmonella status)
  annotate("text", x =0.3, y = -0.2, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces.cDNA_salmonella_BC_adonis_status * 100, 1), "%",
                         "\np = ", round(pvalue_feces.cDNA_salmonella_BC_adonis_status, 3)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
feces.cDNA_salmonella_only_BC_beta_div

####CATCH BASIN SAMLPES#######
data.bacteria.samples.water.tss ##already have this object with normalized counts
data.bacteria.samples.water.tss.bray ##Distance matrix
data.bacteria.samples.water.tss.df #metadata dataframe
data.bacteria.samples.water.tss.bray.ord ##ordination

##metadata
data.bacteria.samples.water.salmonella.tss.df<- data.bacteria.samples.water.tss.df %>%
  mutate (feedlot_salmonella = paste(feedlot, salmonella_culture_status, sep = '_'), #adding a column concatenating "feedlot" and "Salmonella"
          salmonella_stat.abbrv = dplyr::recode(salmonella_culture_status, "positive"= "P", "negative"= "N"),
          feedlot = factor(feedlot)) ##making feedlot a factor

#PERMANOVA
set.seed(87)
water_salmonella_BC_adonis <- adonis2(data.bacteria.samples.water.tss.bray ~ 
                                        salmonella_culture_status + feedlot, by = "margin", 
                                      data.bacteria.samples.water.salmonella.tss.df, p.adjust.methods = "BH", 
                                      permutations = 9999)
water_salmonella_BC_adonis 
#8.6% of variation is due to salmonella_culture_status (positive vs negative) p=  0.0019
#63.7% of variation is due to feedlot p= 0.0001

#PERMANOVA- Salmonella status interaction with feedlot 
set.seed(87)
water_salmonella_feedlot_adonis_interaction <- adonis2(data.bacteria.samples.water.tss.bray ~ salmonella_culture_status*feedlot, 
                                                       by = "margin", 
                                                       data.bacteria.samples.water.salmonella.tss.df, permutations = 9999)
water_salmonella_feedlot_adonis_interaction
#6.4% of variation is due to salmonella_culture_status interaction with feedlot (p= 1e-04)

#PERMANOVA- Salmonella status stratified by feedlot 
set.seed(87)
water_salmonella_BC_adonis_status <- adonis2(data.bacteria.samples.water.tss.bray ~ salmonella_culture_status, 
                                             by = "margin", 
                                             strata = data.bacteria.samples.water.salmonella.tss.df$feedlot, 
                                             data.bacteria.samples.water.salmonella.tss.df, p.adjust.methods = "BH", permutations = 9999)
water_salmonella_BC_adonis_status
#7.7% of variation is due to salmonella_culture_status (positive vs negative) p= 0.12 (ns)

#PERMANOVA- Feedlot stratified by salmonella status
set.seed(87)
water_salmonella_BC_adonis_feedlot <- adonis2(data.bacteria.samples.water.tss.bray ~ feedlot, by = "margin", strata = data.bacteria.samples.water.salmonella.tss.df$salmonella_culture_status, data.bacteria.samples.water.salmonella.tss.df, p.adjust.methods = "BH", permutations = 9999)
water_salmonella_BC_adonis_feedlot
#62.6% of variation is due to feedlot p= 1e-04

#####CATCH BASIN DNA#######
data.bacteria.samples.water.DNA.tss <- subset_samples(data.bacteria.samples.water.tss, gen_material =="DNA")
data.bacteria.samples.water.DNA.tss <- prune_taxa(taxa_sums(data.bacteria.samples.water.DNA.tss) > 0, 
                                                   data.bacteria.samples.water.DNA.tss) 

data.bacteria.samples.water.DNA.tss #30235 taxa and 12 samples

#Metadata
data.bacteria.samples.water.DNA.tss.df <- data.frame(data.bacteria.samples.water.DNA.tss@sam_data)

###### Including inly feedlots with positive and negative samples (feedlots 3 and 5 did not have any positive samples)#######
data.bacteria.samples.water.DNA.tss #Only water DNA samples 
data.bacteria.samples.water.DNA.f14.tss <- subset_samples(data.bacteria.samples.water.DNA.tss, !feedlot %in% c("3", "5"))#Subsetting out feedlot 3 and 5 - doesn't have salmonella positive samples
data.bacteria.samples.water.DNA.f14.tss <- prune_taxa(taxa_sums(data.bacteria.samples.water.DNA.f14.tss) > 0, data.bacteria.samples.water.DNA.f14.tss) 
data.bacteria.samples.water.DNA.f14.tss ##28507 taxa and 9 samples

##Distance matrix
data.bacteria.samples.water.DNA.f14.tss.bray <- vegdist(t(data.bacteria.samples.water.DNA.f14.tss@otu_table), method = "bray") 
#Ordination 
set.seed(87)
data.bacteria.samples.water.DNA.f14.tss.bray.ord <- metaMDS(data.bacteria.samples.water.DNA.f14.tss.bray, k = 2, try= 20, trymax= 1000, autotransform = F)

##metadata
data.bacteria.samples.water.DNA.f14.salmonella.tss.df<- data.bacteria.samples.water.DNA.tss.df %>%
  filter(!feedlot %in% c("3", "5"))%>%
  mutate (feedlot_salmonella = paste(feedlot, salmonella_culture_status, sep = '_'))%>% ##adding a column concatenating "feedlot" and "Salmonella"
  mutate (salmonella_stat.abbrv = dplyr::recode(salmonella_culture_status, "positive"= "P", "negative"= "N"),
          feedlot = factor(feedlot))

#### ADDING CENTROIDS FOR PLOTTING
## BC- centroids on only salmonella status
data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.plot <- ordiplot(data.bacteria.samples.water.DNA.f14.tss.bray.ord$points)
data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.scrs <- scores(data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.plot, display = "sites") 
data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.scrs <- cbind(as.data.frame(data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.scrs), 
                                                           salmonella_culture_status = data.bacteria.samples.water.DNA.f14.salmonella.tss.df$salmonella_culture_status,
                                                           salmonella_stat.abbrv = data.bacteria.samples.water.DNA.f14.salmonella.tss.df$salmonella_stat.abbrv ) #This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ salmonella_culture_status, data = data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.scrs, FUN = mean) ##Centroids according to gen_material_spec_2 (cDNA/DNA_feedlot)
data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.segs <- merge(data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.scrs, setNames(data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.cent, c("salmonella_culture_status","cMDS1","cMDS2")), by = 'salmonella_culture_status', sort = F)

#PERMANOVA 
set.seed(87)
water.DNA_salmonellaonly_BC_adonis <- adonis2(data.bacteria.samples.water.DNA.f14.tss.bray ~ salmonella_culture_status + feedlot, 
                                              #strata = data.bacteria.samples.water.DNA.f14.salmonella.tss.df$feedlot,
                                              data.bacteria.samples.water.DNA.f14.salmonella.tss.df,
                                              by = "margin",
                                              p.adjust.methods = "BH", 
                                              permutations = 9999)
water.DNA_salmonellaonly_BC_adonis #10.8% of variation is due to Salmonella Status (p= 0.145)

######SUPPLEMENTARY TABLE 5.3#######
stable5.3 <- data.frame(water.DNA_salmonellaonly_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)",
         `Library Type` = "Metagenomic (DNA)", 
         `Sample Type` = "Catch Basins",
         SumOfSqs = round(SumOfSqs, 2), 
         `R2` = round(`R2`,3), 
         `F` = round(`F`,2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3))%>%
  rownames_to_column(var = "Fixed Effect")%>%
  mutate(`Fixed Effect` = str_replace(`Fixed Effect`, "salmonella_culture_status", "Salmonella Culture Status"), 
         `Fixed Effect` = str_replace(`Fixed Effect`, "feedlot", "Feedlot"))%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))
stable5.3

#PERMDISP- Salmonella status
# Run the betadisper function, average distance to centroid
bray.water.DNA.f14.salmonella.disp <- betadisper(data.bacteria.samples.water.DNA.f14.tss.bray, data.bacteria.samples.water.DNA.f14.salmonella.tss.df$salmonella_culture_status)
bray.water.DNA.f14.salmonella.disp
##Then test by permuting
set.seed(87)
bray.water.DNA.f14.salmonella.permdisp <- permutest(bray.water.DNA.f14.salmonella.disp, permutations = 9999)
bray.water.DNA.f14.salmonella.permdisp 
##No different dispersions of variance between salmonella status (p = 0.75)


#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.water.DNA.f14.feedlot.disp <- betadisper(data.bacteria.samples.water.DNA.f14.tss.bray, 
                                              data.bacteria.samples.water.DNA.f14.salmonella.tss.df$feedlot)
bray.water.DNA.f14.feedlot.disp
##Then test by permuting
set.seed(87)
bray.water.DNA.f14.feedlot.permdisp <- permutest(bray.water.DNA.f14.feedlot.disp, permutations = 9999)
bray.water.DNA.f14.feedlot.permdisp
##No different dispersions of variance between feedlots (p = 0.45)

# Extract R2 and p-values
#Salmonella 
R2_water.DNA_salmonella_BC_adonis_status <- water.DNA_salmonellaonly_BC_adonis$R2[1] 
pvalue_water.DNA_salmonella_BC_adonis_status<-water.DNA_salmonellaonly_BC_adonis$`Pr(>F)`[1]

#Plot
water.DNA_salmonella_only_BC_beta_div <- ggplot(data.bacteria.samples.water.DNA.f14.salmonella.only.tss.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "CB (DNA)", color = "Salmonella Status") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = salmonella_culture_status))+
  geom_point(aes(x=MDS1, y=MDS2, colour = salmonella_culture_status), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = salmonella_culture_status), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2, label= salmonella_stat.abbrv), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values= salmonella.palette)+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 22),
        legend.text = element_text(colour = "black", size = 22),
        legend.title.position = "top",
        legend.title.align = 0.5,
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+ 
  guides(
    color = guide_legend (override.aes = list(size = 7)))+
  annotate("text", x = -0.3, y = -0.5, ##change coordinates as needed
           label = "Salmonella Status",
           hjust = 0.5, vjust = -0.5, size = 5, colour = "black", fontface = "bold") + ##annotate variable (SAlmonella status)
  annotate("text", x =-0.3, y = -0.5, ##change coordinates as needed
           label = paste("R² = ", round(R2_water.DNA_salmonella_BC_adonis_status * 100, 1), "%",
                         "\np = ", round(pvalue_water.DNA_salmonella_BC_adonis_status, 3)),
           hjust = 0.5, vjust = 1.1, size = 5, colour = "black") # Annotate R² and p-values
water.DNA_salmonella_only_BC_beta_div


####CATCH BASIN cDNA#######
data.bacteria.samples.water.cDNA.tss <- subset_samples(data.bacteria.samples.water.tss, gen_material =="DNA")
data.bacteria.samples.water.cDNA.tss <- prune_taxa(taxa_sums(data.bacteria.samples.water.cDNA.tss) > 0, 
                                                  data.bacteria.samples.water.cDNA.tss) 

data.bacteria.samples.water.cDNA.tss #30235 taxa and 12 samples

#Metadata 
data.bacteria.samples.water.cDNA.tss.df <- data.frame(data.bacteria.samples.water.cDNA.tss@sam_data)
###### Taking out feedlots 3 and 5 since did not have any positive samples#######
data.bacteria.samples.water.cDNA.tss #Only water cDNA samples 
data.bacteria.samples.water.cDNA.f14.tss <- subset_samples(data.bacteria.samples.water.cDNA.tss, !feedlot %in% c("3", "5")) #Subsetting out feedlot 5 - doesn't have salmonella positive samples
data.bacteria.samples.water.cDNA.f14.tss <- prune_taxa(taxa_sums(data.bacteria.samples.water.cDNA.f14.tss) > 0, data.bacteria.samples.water.cDNA.f14.tss) 
data.bacteria.samples.water.cDNA.f14.tss ##28507 taxa and 9 samples

##Distance matrix
data.bacteria.samples.water.cDNA.f14.tss.bray <- vegdist(t(data.bacteria.samples.water.cDNA.f14.tss@otu_table), method = "bray") 
#Ordination 
set.seed(87)
data.bacteria.samples.water.cDNA.f14.tss.bray.ord <- metaMDS(data.bacteria.samples.water.cDNA.f14.tss.bray, k = 2, try= 20, trymax= 1000, autotransform = F)

##metadata
data.bacteria.samples.water.cDNA.f14.salmonella.tss.df<- data.bacteria.samples.water.cDNA.tss.df %>%
  filter(!feedlot %in% c("3", "5"))%>%
  mutate (feedlot_salmonella = paste(feedlot, salmonella_culture_status, sep = '_'))%>% ##adding a column concatenating "feedlot" and "Salmonella"
  mutate (salmonella_stat.abbrv = dplyr::recode(salmonella_culture_status, "positive"= "P", "negative"= "N"),
          feedlot = factor(feedlot))

#### ADDING CENTROIDS FOR PLOTTING
## BC- centroids on only salmonella status
data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.plot <- ordiplot(data.bacteria.samples.water.cDNA.f14.tss.bray.ord$points)
data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.scrs <- scores(data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.plot, display = "sites") 
data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.scrs <- cbind(as.data.frame(data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.scrs), 
                                                            salmonella_culture_status = data.bacteria.samples.water.cDNA.f14.salmonella.tss.df$salmonella_culture_status,
                                                            salmonella_stat.abbrv = data.bacteria.samples.water.cDNA.f14.salmonella.tss.df$salmonella_stat.abbrv ) #This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ salmonella_culture_status, data = data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.scrs, FUN = mean) ##Centroids according to gen_material_spec_2 (ccDNA/cDNA_feedlot)
data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.segs <- merge(data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.scrs, setNames(data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.cent, c("salmonella_culture_status","cMDS1","cMDS2")), by = 'salmonella_culture_status', sort = F)

#PERMANOVA 
set.seed(87)
water.cDNA_salmonellaonly_BC_adonis <- adonis2(data.bacteria.samples.water.cDNA.f14.tss.bray ~ salmonella_culture_status + feedlot, 
                                               data.bacteria.samples.water.cDNA.f14.salmonella.tss.df,
                                               by = "margin",
                                               p.adjust.methods = "BH", 
                                               permutations = 9999)
water.cDNA_salmonellaonly_BC_adonis #10.8% of variation is due to Salmonella Status (p= 0.14)

######SUPPLEMENTARY TABLE 5.4#######
stable5.4 <- data.frame(water.cDNA_salmonellaonly_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)",
         `Library Type` = "Metatranscriptomic (RNA (cDNA))", 
         `Sample Type` = "Catch Basins",
         SumOfSqs = round(SumOfSqs, 2), 
         `R2` = round(`R2`,3), 
         `F` = round(`F`,2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3))%>%
  rownames_to_column(var = "Fixed Effect")%>%
  mutate(`Fixed Effect` = str_replace(`Fixed Effect`, "salmonella_culture_status", "Salmonella Culture Status"), 
         `Fixed Effect` = str_replace(`Fixed Effect`, "feedlot", "Feedlot"))%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))
stable5.4      

#PERMDISP- Salmonella status
# Run the betadisper function, average distance to centroid
bray.water.cDNA.f14.salmonella.disp <- betadisper(data.bacteria.samples.water.cDNA.f14.tss.bray, data.bacteria.samples.water.cDNA.f14.salmonella.tss.df$salmonella_culture_status)
bray.water.cDNA.f14.salmonella.disp
##Then test by permuting
set.seed(87)
bray.water.cDNA.f14.salmonella.permdisp <- permutest(bray.water.cDNA.f14.salmonella.disp, permutations = 9999)
bray.water.cDNA.f14.salmonella.permdisp 
##No different dispersions of variance between salmonella status (p = 0.75)

#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.water.cDNA.f14.feedlot.disp <- betadisper(data.bacteria.samples.water.cDNA.f14.tss.bray, 
                                              data.bacteria.samples.water.cDNA.f14.salmonella.tss.df$feedlot)
bray.water.cDNA.f14.feedlot.disp
##Then test by permuting
set.seed(87)
bray.water.cDNA.f14.feedlot.permdisp <- permutest(bray.water.cDNA.f14.feedlot.disp, permutations = 9999)
bray.water.cDNA.f14.feedlot.permdisp
##No different dispersions of variance between feedlots (p = 0.45)

# Extract R2 and p-values
#Salmonella 
R2_water.cDNA_salmonella_BC_adonis_status <- water.cDNA_salmonellaonly_BC_adonis$R2[1] 
pvalue_water.cDNA_salmonella_BC_adonis_status<-water.cDNA_salmonellaonly_BC_adonis$`Pr(>F)`[1]

##Plot 
water.cDNA_salmonella_only_BC_beta_div <- ggplot(data.bacteria.samples.water.cDNA.f14.salmonella.only.tss.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "CB (RNA(cDNA))", color = "Salmonella Status") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = salmonella_culture_status))+
  geom_point(aes(x=MDS1, y=MDS2, colour = salmonella_culture_status), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = salmonella_culture_status), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2, label= salmonella_stat.abbrv), colour= "white", size = 5, fontface = "bold") +
  scale_color_manual(values= salmonella.palette)+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 22),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+ 
  guides(
    color = guide_legend (override.aes = list(size = 7)))+
  annotate("text", x =-0.5, y = -0.7, ##change coordinates as needed
           label = "Salmonella Status",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (SAlmonella status)
  annotate("text", x =-0.5, y = -0.7, ##change coordinates as needed
           label = paste("R² = ", round(R2_water.cDNA_salmonella_BC_adonis_status * 100, 1), "%",
                         "\np = ", round(pvalue_water.cDNA_salmonella_BC_adonis_status, 3)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
water.cDNA_salmonella_only_BC_beta_div

#PLOT - Salmonella status + feedlot - Used for Figure 2 
WaterandFeces.DNA_salmonella_only_BC_beta_div <- ggarrange(feces.DNA_salmonella_only_BC_beta_div+ 
                                                             theme(legend.position = "none", 
                                                                   axis.title = element_text(size = 20),
                                                                   plot.title = element_text(size = 30),),
                                                           water.DNA_salmonella_only_BC_beta_div+ 
                                                             theme(legend.position = "none", 
                                                                   axis.title = element_text(size = 20),
                                                                   plot.title = element_text(size = 30)), 
                                                           labels = c("B", "C"), 
                                                           font.label = list(size = 22, color = "black", face = "bold", family = NULL),
                                                           ncol = 1,
                                                           common.legend = T,
                                                           legend = "bottom")
WaterandFeces.DNA_salmonella_only_BC_beta_div

###SUPPLEMENTARY TABLE 5 - MICROBIOME SECTION#######
stable5_microbiome <- bind_rows(stable5.1,
                                stable5.2,
                                stable5.3, 
                                stable5.4,
                                stable5.5,
                                stable5.6,
                                stable5.7,
                                stable5.8)%>%
  select(Dataset, `Library Type`, `Sample Type`, `Fixed Effect`, Df, SumOfSqs, R2, `F`, `Pr(>F)`)
stable5_microbiome 
write_xlsx(stable5_microbiome, 
           "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryTable5_1.xlsx")

#RELATIVE ABUNDANCE AND DENDROGRAMS#####
##All samples####
data.bacteria.samples.tss.bray <- vegdist(t(data.bacteria.samples.tss@otu_table), method = "bray") 
data.bacteria.samples.bray.hclust <- hclust(data.bacteria.samples.tss.bray,
                            method = "ward.D2")
plot(data.bacteria.samples.bray.hclust, hang = -1)
data.bacteria.samples.bray.dendro <- as.dendrogram(data.bacteria.samples.bray.hclust) # Build dendrogram object from hclust results
data.bacteria.samples.bray.dendro.data <- dendro_data(data.bacteria.samples.bray.dendro, type = "rectangle") # Extract the dendrogram plot data
data.bacteria.samples.tss@sam_data$sampleID <- rownames(data.bacteria.samples.tss@sam_data) ##need a column with sample_ID
data.bacteria.samples.bray.dendro.metadata <- 
  data.frame(data.bacteria.samples.tss@sam_data) %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"),
         gen_material.abbrv = dplyr::recode(gen_material, "cDNA"= "T", "DNA"= "G"))
data.bacteria.samples.bray.dendro.data$labels <- data.bacteria.samples.bray.dendro.data$labels %>%
  left_join(data.bacteria.samples.bray.dendro.metadata, by = c("label" = "sampleID")) 

###Dendrogram#######
dendro.bray.plot <- ggplot(data.bacteria.samples.bray.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  #Sample type
  geom_point(data = data.bacteria.samples.bray.dendro.data$labels, 
    aes(x = x, y = y, colour = sample_type),size = 5, shape = 15, 
    position = position_nudge(y = -0.2))+
  scale_color_manual(name = "Sample Type", values = sample.type.palette,
                     label = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  guides(color=guide_legend(title.position="top"))+
  new_scale_color()+ 
  #Feedlots
  geom_point(
    data = data.bacteria.samples.bray.dendro.data$labels, 
    aes(x = x, y = y, colour = factor(feedlot)),size = 5, shape = 15, 
    position = position_nudge(y = -0.5))+
  scale_color_manual(name = "Feedlot", values = feedlot_palette) +
  guides(color=guide_legend(title.position="top"))+
  new_scale_color()+
  #Library Type
  geom_point(data = data.bacteria.samples.bray.dendro.data$labels,  
             aes(x = x, y = y, color = gen_material),
             size = 5, shape = 15, 
             position = position_nudge(y = -0.8)) +
  scale_color_manual(name = "Library Type",values = gen.material.palette,  
                     label = c("DNA" = "DNA", "cDNA" = "RNA(cDNA)")) +
  #Sample type
  geom_text(data = data.bacteria.samples.bray.dendro.data$labels, 
    aes(x=x, y=y, label = sample_type.abbrv), 
    colour = "white", size =3, 
    position = position_nudge(y=-0.2), 
    fontface = "bold")+
  #Feedlot
  geom_text(data = data.bacteria.samples.bray.dendro.data$labels, 
            aes(x=x, y=y, label = factor(feedlot)), 
            colour = "white", size =4, 
            position = position_nudge(y=-0.5), fontface = "bold") +
  #Library Type
  geom_text(data = data.bacteria.samples.bray.dendro.data$labels, 
            aes(x=x, y=y, label = gen_material.abbrv), 
            colour = "white", size =4, 
            position = position_nudge(y=-0.8), fontface = "bold") +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  guides(color=guide_legend(title.position="top"))+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold", size = 24),
        legend.text = element_text(size = 25),
        plot.title = element_text(size = 36),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 21),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
dendro.bray.plot

##PHYLUM LEVEL####
dendro_bray_data.bacteria.samples_order <- data.bacteria.samples.bray.dendro.data$labels$label ##order the samples
data.bacteria.samples.phylum.filt <- merge_low_abundance_grouped_ra(data.bacteria.samples_phylum.ra, "sample_type", level = "Phylum", threshold = 0.5)
data.bacteria.samples.phylum.filt #17 phyla over 0.5% mean RA
data.bacteria.samples.phylum.filt.melt <- psmelt(data.bacteria.samples.phylum.filt)%>%
  mutate(Phylum = factor(Phylum, 
                        levels = c(setdiff(Phylum, 
                                           unique(grep("Others", Phylum, value = TRUE))), 
                                   unique(grep("Others", Phylum, value = TRUE)))))##Factoring the Phylum column so that "Others.." is the last category
levels(data.bacteria.samples.phylum.filt.melt$Phylum) ##ok

##Which are the top most abundant taxa by group? 
data.bacteria.samples.phylum.filt.melt %>%
  group_by(sample_type, Phylum) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(sample_type,  desc(mean_abun))%>%
  print(n=40)

#Per library type
most_abun_phyla <- data.bacteria.samples.phylum.filt.melt %>%
  group_by(sample_type, gen_material, Phylum) %>%
  summarise(
    mean_relative_abundance = mean(Abundance, na.rm = TRUE),
    sd = sd(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(sample_type, gen_material, desc(mean_relative_abundance))%>%
  mutate(`Sample Type` = str_replace(sample_type, "Water", "Catch Basins"), 
         `Library Type` = str_replace_all(gen_material, 
                                          c("cDNA" = "Metatranscriptomic (RNA (cDNA))", 
                                            "(?<!c)DNA" = "Metagenomic (DNA)") #"Match DNA only if it is NOT immediately preceded by c." 
                                               ), 
         `Mean Relative Abundance (%) ± SD` = paste0(
           round(mean_relative_abundance, 2),
           " ± ",
           round(sd, 3)))%>%
  select(-c("gen_material", "sample_type", "mean_relative_abundance", "sd"))

#####SUPPLEMENTARY TABLE 8.1####
stable8.1 <- most_abun_phyla%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)")%>%
  select(Dataset, `Library Type`, `Sample Type`, `Mean Relative Abundance (%) ± SD`, everything())
write_xlsx(stable8.1,
           "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryTable8_1.xlsx")

##Create color palette
set.seed(87)
phylum.filt.palette <- distinctColorPalette(length(unique(data.bacteria.samples.phylum.filt.melt$Phylum)))
phylum_filt_names <- unique(data.bacteria.samples.phylum.filt.melt$Phylum)# Create a named vector for the palette, where the names correspond to phlyum names
phylum_named_palette <- setNames((phylum.filt.palette)[1:length(phylum_filt_names)], phylum_filt_names)
phylum_named_palette$'Others <0.5% RA' <- "grey95"

#Plot
dendroRA.phylum.plot <- ggplot(data.bacteria.samples.phylum.filt.melt, aes(x=Sample, y= Abundance, fill = Phylum)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_bray_data.bacteria.samples_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =phylum_named_palette) +
  guides(fill=guide_legend(title.position="top", nrow = 3))+
  theme(legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 22),
    legend.text = element_text(size = 22),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
    axis.title.y = element_text(size = 25),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()) 
dendroRA.phylum.plot

#Put dendrogram together with relative abundance
dendroRA.phylum.plot.2 <- plot_grid(dendro.bray.plot, 
                                    dendroRA.phylum.plot, 
                                    align = "v", 
                                    ncol = 1, 
                                    rel_heights = c(0.3, 0.7))

dendroRA.phylum.plot.2

###FIGURE 4ABC - FECES VS WATER #####
#Want to add these on the side of the dendrogram
cDNA_BC_beta_div_spider_sampletype
DNA_BC_beta_div_spider_sampletype

#cDNA and DNA
cDNAandDNA_BC_beta_div_spider_sampletype <- ggarrange(
  DNA_BC_beta_div_spider_sampletype + 
    theme(
          # plot.title = element_text(size = 35),
          plot.title = element_blank(),
          axis.title = element_text(size = 20)),
  cDNA_BC_beta_div_spider_sampletype +
    theme(
      # plot.title = element_text(size = 35),
      plot.title = element_blank(),
      axis.title = element_text(size = 20)),
  align = "h",
  nrow = 1,
  common.legend = TRUE,
  labels = c("A", "B"),
  font.label = list(size = 30),
  legend = "none"
)
cDNAandDNA_BC_beta_div_spider_sampletype

##DNA and cDNA ordination plots, on top of dendrogram at the phylum level - Figure 5A&B
figure4ABC <- plot_grid(cDNAandDNA_BC_beta_div_spider_sampletype, 
                     dendroRA.phylum.plot.2, 
                     align = "v", 
                     ncol = 1,
                     labels = c(" ", "C"),
                     label_size = 30,
                     rel_heights = c(0.25, 0.75))
  #labs(title = "MICROBIOME")+
  #theme(plot.title = element_text(size = 30, face = "bold"))
figure4ABC
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure4ABC.svg", 
       plot = figure4ABC, 
       device = "svg", 
       #dpi = 600,
       width = 30, 
       height = 14, 
       bg = "white")


##ORDER LEVEL####
data.bacteria.samples.order.filt <- merge_low_abundance_grouped_ra(data.bacteria.samples_order.ra, "sample_type", level = "Order", 
                                                   threshold = 0.5)
data.bacteria.samples.order.filt #46 orders over 0.5% RA
data.bacteria.samples.order.filt.melt <- psmelt(data.bacteria.samples.order.filt)%>%
  mutate(Order = factor(Order, 
                        levels = c(
                          unlist(
                            lapply(unique(Class), function(c) {
                              setdiff(unique(Order[Class == c]), unique(grep("Others", Order, value = TRUE)))
                            })
                          ), 
                          unique(grep("Others", Order, value = TRUE))
                        ))) ##Factoring Order by Class and making sure "Others" is last 

levels(data.bacteria.samples.order.filt.melt$Order)

##Apply the function to obtain top orders (n=20)
top_orders <- top_taxa_legend(data.bacteria.samples.order.filt.melt, taxlevel = "Order", n = 25)
top_orders

##What are the abundances of those top orders
data.bacteria.samples.order.filt.melt %>%
  group_by(sample_type, Order) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(sample_type, desc(mean_abun))%>%
  tail(n = 200)%>%
  print(n=200)

#Per library type:
most_abun_orders <- data.bacteria.samples.order.filt.melt %>%
  group_by(sample_type, gen_material, Order) %>%
  summarise(
  mean_relative_abundance = mean(Abundance, na.rm = TRUE),
  sd = sd(Abundance, na.rm = TRUE),
  .groups = "drop") %>%
  arrange(gen_material, sample_type, desc(mean_relative_abundance))%>%
  mutate(`Sample Type` = str_replace(sample_type, "Water", "Catch Basins"), 
         `Library Type` = str_replace_all(gen_material, 
                                          c("cDNA" = "Metatranscriptomic (RNA (cDNA))", 
                                            "(?<!c)DNA" = "Metagenomic (DNA)") #"Match DNA only if it is NOT immediately preceded by c." 
         ), 
         `Mean Relative Abundance (%) ± SD` = paste0(
           round(mean_relative_abundance, 2),
           " ± ",
           round(sd, 3)))%>%
  select(-c("gen_material", "sample_type", "mean_relative_abundance", "sd"))
most_abun_orders

##Making color palette
# Function to generate a distinct color palette with a specified size, excluding colors in the phylum.filt.palette
generate_palette <- function(num_colors, exclude_colors) {
  # Generate a large enough distinct color palette (for example, 100 colors)
  candidate_palette <- distinctColorPalette(num_colors)
  
  # Exclude colors already present in phylum.filt.palette
  remaining_colors <- setdiff(candidate_palette, exclude_colors)
  
  # If there are fewer than num_colors, regenerate more distinct colors
  while (length(remaining_colors) < num_colors) {
    # Generate more colors and remove already used ones
    additional_colors <- distinctColorPalette(100)
    remaining_colors <- setdiff(c(remaining_colors, additional_colors), exclude_colors)
  }
  
  # Return exactly num_colors
  return(remaining_colors[1:num_colors])
}
##Make the color palette
order.filt.palette <- generate_palette(length(unique(data.bacteria.samples.order.filt.melt$Order)), phylum.filt.palette)# Generate the order palette with 48 colors, excluding phylum.filt.palette colors
order_filt_names <- unique(data.bacteria.samples.order.filt.melt$Order) ##extract order taxa names
order_named_palette <- setNames((order.filt.palette)[1:length(order_filt_names)], order_filt_names)# Create a named vector for the palette, where the names correspond to order names

##FINDING MATCHING PHYLUM-ORDER NAMES###
##Some taxa names include "unclassified", want to remove that so I can find if, for example, the Phylum is "Bacillota" and the order is "unclassified Bacillota"
standardize_names <- function(names_vector) {
  # Remove prefixes such as "unclassified" or similar, keeping the main taxonomy
  gsub("^(unclassified\\s+)?", "", names_vector, ignore.case = TRUE)
}
# Standardize the names in both palettes
standardized_order_names <- standardize_names(names(order_named_palette))
standardized_phylum_names <- standardize_names(names(phylum_named_palette))

# Match names based on the "standardized" versions
matching_names_phyla_order_st <- intersect(standardized_order_names, standardized_phylum_names)

# Match colors for orders that have the (exact) same name as their corresponding phylum
matching_names_phyla_order <- intersect(names(order_named_palette), names(phylum_named_palette)) ##unclassified in both order and species

##The rest that are matching are these ones, which are unclassified at the order level
phylum_named_palette[matching_names_phyla_order_st]##Orders that are unclassified but have a phylum 

# Assign colors from the phylum palette to the order palette for matching names
for (name in matching_names_phyla_order) {
  order_named_palette[name] <- phylum_named_palette[name]
}# Now order_named_palette will have colors matching those in phylum_named_palette for identical names

##Changing the other unclassified colors to make sure they match their respective phylum (check phylum_named_palette[matching_names_phyla_order_st])
order_named_palette$'Others <0.5% RA' <- "grey95" #give Others grey color
order_named_palette$'unclassified Pseudomonadota' <- "#E3BFBD"
order_named_palette$'unclassified Bacillota' <-  "#DB8350"
order_named_palette$'unclassified Bacteroidota' <- "#D3EAD5"
order_named_palette$'unknown Bacteria [uncultured bacterium]' <- "#E45189"


#Dendrogram plot
dendroRA.order.plot <- ggplot(data.bacteria.samples.order.filt.melt, aes(x=Sample, y= Abundance, fill = Order)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_bray_data.bacteria.samples_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = order_named_palette, breaks = top_orders) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 22),
        legend.text = element_text(size = 22),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  guides(
    fill = guide_legend(
      nrow = 5,
      title.position = "top"
      # Adjust the height of the legend key
    )
  )
dendroRA.order.plot

#Dendrogram together with RA plot at the order level
dendroRA.order.plot.2 <- plot_grid(dendro.bray.plot, dendroRA.order.plot, 
                                   align = "v", ncol = 1,
                                   rel_heights = c(0.3, 0.7))
dendroRA.order.plot.2


##FAMILY ####
data.bacteria.samples.family.filt <- merge_low_abundance_grouped_ra(data.bacteria.samples_family.ra, 
                                                    "sample_type", 
                                                    level = "Family", 
                                                    threshold = 0.5) ##merge the low abundance taxa into "Others"
data.bacteria.samples.family.filt #59 families over 0.3% RA 

data.bacteria.samples.family.filt.melt <- psmelt(data.bacteria.samples.family.filt)%>%
  mutate(Family = factor(Family, 
                        levels = c(
                          unlist(
                            lapply(unique(Order), function(c) {
                              setdiff(unique(Family[Order == c]), unique(grep("Others", Family, value = TRUE)))
                            })
                          ), 
                          unique(grep("Others", Family, value = TRUE))
                        ))) ##Factoring the Family column so that "Others.." is the last category, also factoring families by Order
levels(data.bacteria.samples.family.filt.melt$Family)

##Apply the function to obtain top orders (n=10)
top_families <- top_taxa_legend(data.bacteria.samples.family.filt.melt , taxlevel = "Family", n = 25)
top_families

##Making the color palette for family
family.filt.palette <- generate_palette(length(unique(data.bacteria.samples.family.filt.melt$Family)), order.filt.palette) # Create the distinct color palette for families
family_filt_names <- unique(data.bacteria.samples.family.filt.melt$Family) ##extract family taxa names
family_named_palette <- setNames((family.filt.palette)[1:length(family_filt_names)], family_filt_names)# Create a named vector for the palette, where the names correspond to family names

###FINDING MATCHING NAMES###
# Standardize the names 
standardized_family_names <- standardize_names(names(family_named_palette))

# Match names based for the "standardized" versions
matching_names_order_family_st <- intersect(standardized_order_names, standardized_family_names)
####

# Match colors for families that have the (exact) same name as their corresponding order
matching_names_order_family <- intersect(names(family_named_palette), names(order_named_palette)) ##unclassified in both order and family

##The rest that are matching are these 
order_named_palette[matching_names_order_family_st]##Families that are unclassified but have an order 

# Assign colors from the order palette to the family palette for matching names
for (name in matching_names_order_family) {
  family_named_palette[name] <- order_named_palette[name]
}# Now family_named_palette will have colors matching those in order_named_palette for identical names

##Changing the other unclassified colors to make sure they match their respective order (check order_named_palette[matching_names_order_family_st])
family_named_palette$'Others <0.3% RA' <- "grey95" #give Others grey color
family_named_palette$'unclassified Bacteroidales' <- "#A131E5"
family_named_palette$'unclassified Eubacteriales' <- "#E557AB"
family_named_palette$'unknown Bacteria [uncultured bacterium]' <- "#E45189"
family_named_palette$'unclassified Bacteroidales' <- "#D3EAD5"
family_named_palette$'Paenibacillaceae' <- "#E7E7CB"
family_named_palette$'Chromatiaceae' <- "dodgerblue"
family_named_palette$'Bacteroidaceae' <- "aquamarine"
family_named_palette$'Treponemataceae' <- "gold"
family_named_palette$'Clostridiaceae' <- "maroon"
family_named_palette$'Bacillaceae' <- "seagreen"
##Dendrogram
dendroRA.family.plot <- ggplot(data.bacteria.samples.family.filt.melt, aes(x=Sample, y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_bray_data.bacteria.samples_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =family_named_palette, breaks = top_families) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 22),
        legend.text = element_text(size = 22),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  guides(fill = guide_legend(nrow = 4, title.position = "top", override.aes = list(size =3)))
dendroRA.family.plot

#Dendrogram together with RA
dendroRA.family.plot.2 <- plot_grid(dendro.bray.plot, 
                                    dendroRA.family.plot, 
                                    align = "v", ncol = 1,
                                    rel_heights = c(0.3, 0.7))
dendroRA.family.plot.2

##What are the abundances of those top families
data.bacteria.samples.family.filt.melt %>%
  group_by(sample_type, Family) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(sample_type, desc(mean_abun))%>%
  tail(n = 200)%>%
  print(n=200)

#Per library type (40 most abundant families)
most_abun_families <- data.bacteria.samples.family.filt.melt %>%
  group_by(gen_material, sample_type, Family) %>%
  summarise(
    mean_relative_abundance = mean(Abundance, na.rm = TRUE),
    sd = sd(Abundance, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(gen_material, sample_type, desc(mean_relative_abundance))%>%
  group_by(gen_material, sample_type) %>%
  slice_max(order_by = mean_relative_abundance, n = 40) %>%
  ungroup()%>%
  mutate(`Sample Type` = str_replace(sample_type, "Water", "Catch Basins"), 
         `Library Type` = str_replace_all(gen_material, 
                                          c("cDNA" = "Metatranscriptomic (RNA (cDNA))", 
                                            "(?<!c)DNA" = "Metagenomic (DNA)") #"Match DNA only if it is NOT immediately preceded by c." 
         ), 
         `Mean Relative Abundance (%) ± SD` = paste0(
           round(mean_relative_abundance, 2),
           " ± ",
           round(sd, 3)))%>%
  select(-c("gen_material", "sample_type", "mean_relative_abundance", "sd"))
most_abun_families

#####SUPPLEMENTARY TABLE 9_1####
stable9.1 <- most_abun_families%>%
  mutate(Dataset = "Microbiome (Bacteria - Archaea)")%>%
  select(Dataset, `Library Type`, `Sample Type`, `Mean Relative Abundance (%) ± SD`, everything()) 
write_xlsx(stable9.1, 
           "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryTable9_1.xlsx")

###FIGURE 7ABC - EMPHASIS ON LIBRARY TYPE #######
water_genmat_BC_beta_div 
feces_genmat_BC_beta_div 

WaterandFeces_beta_div_spider_genmat <- ggarrange(
  feces_genmat_BC_beta_div + 
    theme(
      # plot.title = element_text(size = 35),
      plot.title = element_blank(),
      axis.title = element_text(size = 20)),
  water_genmat_BC_beta_div + 
    theme(
      # plot.title = element_text(size = 35),
      plot.title = element_blank(),
      axis.title = element_text(size = 20)),
  align = "h",
  nrow = 1,
  common.legend = TRUE,
  labels = c("A", "B"),
  font.label = list(size = 30),
  legend = "none"
)
WaterandFeces_beta_div_spider_genmat

##Add to dendrogram
figure7ABC<- plot_grid(WaterandFeces_beta_div_spider_genmat, 
                           dendroRA.family.plot.2,
                      align = "v", 
                      ncol = 1,
                      labels = c(" ", "C"),
                      label_size = 30,
                      rel_heights = c(0.25, 0.75))
  # labs(title = "MICROBIOME")+
  # theme(plot.title = element_text(size = 30, face = "bold"))
figure7ABC
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure7ABC.svg", 
       plot = figure7ABC, 
       device = "svg",
       dpi = 600,
       width = 30, 
       height =14, 
       bg = "white")


##GENUS LEVEL####
data.bacteria.samples.genus.filt <- merge_low_abundance_grouped_ra(data.bacteria.samples_genus.ra, "sample_type", 
                                                   level = "Genus", threshold = 0.15) ##merge the low abundance taxa into "Others"
data.bacteria.samples.genus.filt #155 genera over 0.15% RA 

#Looks like there are duplicates:
length(unique(phyloseq::tax_table(data.bacteria.samples.genus.filt)[, "Genus"])) #155
genus_vec <- as.character(phyloseq::tax_table(data.bacteria.samples.genus.filt)[, "Genus"])
duplicated_genera <- unique(genus_vec[duplicated(genus_vec)]) #Flintibacter

#Melt and grouping (so the duplicated genus will merge) (Summarize abundance per Sample + Genus)
summed_abundance_genus_noreps <- psmelt(data.bacteria.samples.genus.filt) %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")
nrow(summed_abundance_genus_noreps) #down to 18720 (Flintibacter is now present only once per sample)

#Get representative metadata per Sample + Genus
data.bacteria.samples.genus.metadata.noreps <- psmelt(data.bacteria.samples.genus.filt) %>%
  select(-Abundance)%>%
  group_by(Sample, Genus) %>%
  summarise(across(where(is.character), ~ first(.x)),
            across(where(is.numeric), ~ first(.x)),
            .groups = "drop")
nrow(data.bacteria.samples.genus.metadata.noreps) #down to 18600
#Join summed Abundance to genus metadata
data.bacteria.samples.genus.filt.melt <- summed_abundance_genus_noreps %>%
  left_join(data.bacteria.samples.genus.metadata.noreps, by = c("Sample", "Genus"))
nrow(data.bacteria.samples.genus.filt.melt) #18600

#Had to update since there are of duplicates
data.bacteria.samples.genus.filt.melt <- data.bacteria.samples.genus.filt.melt %>%
  mutate(Genus = factor(
    Genus,
    levels = unique(c(
      unlist(
        lapply(unique(Family), function(fam) {
          setdiff(unique(Genus[Family == fam]), grep("Others", Genus[Family == fam], value = TRUE))
        })
      ),
      grep("Others", Genus, value = TRUE)
    ))
  ))

##Apply the function to obtain top orders (n=10)
top_genera <- top_taxa_legend(data.bacteria.samples.genus.filt.melt , taxlevel = "Genus", n = 20)
top_genera
##What are the abundances of those top genera
data.bacteria.samples.genus.filt.melt %>%
  group_by(sample_type, Genus) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(sample_type, desc(mean_abun))%>%
  tail(n = 200)%>%
  print(n=200)

#Per library type
most_abun_genera <- data.bacteria.samples.genus.filt.melt %>%
  group_by(sample_type, gen_material, Genus) %>%
  summarise(
    mean_relative_abundance = mean(Abundance, na.rm = TRUE),
    sd = sd(Abundance, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(gen_material, sample_type, desc(mean_relative_abundance))%>%
  mutate(`Sample Type` = str_replace(sample_type, "Water", "Catch Basins"), 
         `Library Type` = str_replace_all(gen_material, 
                                          c("cDNA" = "Metatranscriptomic (RNA (cDNA))", 
                                            "(?<!c)DNA" = "Metagenomic (DNA)") #"Match DNA only if it is NOT immediately preceded by c." 
         ), 
         `Mean Relative Abundance (%) ± SD` = paste0(
           round(mean_relative_abundance, 2),
           " ± ",
           round(sd, 3)))%>%
  select(-c("gen_material", "sample_type", "mean_relative_abundance", "sd"))
most_abun_genera


##Making the color palette for genus
genus.filt.palette <- generate_palette(length(unique(data.bacteria.samples.genus.filt.melt$Genus)), order.filt.palette) # Create the distinct color palette for genera
genus_filt_names <- unique(data.bacteria.samples.genus.filt.melt$Genus) ##extract genus taxa names
genus_named_palette <- setNames((genus.filt.palette)[1:length(genus_filt_names)], genus_filt_names)# Create a named vector for the palette, where the names correspond to genus names

###FINDING MATCHING NAMES###
# Standardize the names 
standardized_genus_names <- standardize_names(names(genus_named_palette))

# Match names based for the "standardized" versions
matching_names_family_genus_st <- intersect(standardized_family_names, standardized_genus_names)
####

# Match colors for genera that have the (exact) same name as their corresponding family
matching_names_family_genus <- intersect(names(genus_named_palette), names(family_named_palette)) ##unclassified/unknown in both family and genus

##The rest that are matching are these 
family_named_palette[matching_names_family_genus_st]##genera that are unclassified but have a family

# Assign colors from the order palette to the genus palette for matching names
for (name in matching_names_family_genus) {
  genus_named_palette[name] <- family_named_palette[name]
}# Now genus_named_palette will have colors matching those in order_named_palette for identical names

##Changing the other unclassified colors to make sure they match their respective order (check family_named_palette[matching_names_family_genus_st])
genus_named_palette$'Others < 0.15 % RA' <- "grey95" #give zzzOther grey color
genus_named_palette$'unclassified Prevotellaceae' <- "#EAAE67"
genus_named_palette$'unclassified Chromatiaceae' <- "#83E9D2"
genus_named_palette$'unclassified Comamonadaceae' <- "#F26852"
genus_named_palette$'unclassified Oscillospiraceae' <- "#BA99F2"
genus_named_palette$'unclassified Enterobacteriaceae' <- "#59E8E4"
genus_named_palette$'unclassified Lachnospiraceae' <- "#638C8C"
genus_named_palette$'unclassified Bacteroidaceae' <- "#E2A6E8"
genus_named_palette$'unclassified Pseudomonadaceae' <- "#AEF234"
genus_named_palette$'unknown Bacteroidales [Bacteroidales bacterium]' <- "#6867D6"
genus_named_palette$'unknown Bacteroidales [Bacteroidales bacterium MB20-C3-3]' <- "#C6CD8F"
genus_named_palette$'unknown Eubacteriales [uncultured Eubacteriales bacterium]' <- "#ACB959"
genus_named_palette$'unknown Eubacteriales [Clostridiales bacterium]' <- "#B48E6D"
genus_named_palette$'unknown Clostridia [uncultured Clostridia bacterium]' <-  "#BCA49A"

##RA plot
dendroRA.genus.plot <- ggplot(data.bacteria.samples.genus.filt.melt, aes(x=Sample, y= Abundance, fill = Genus)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_bray_data.bacteria.samples_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =genus_named_palette, breaks = top_genera) +
  theme(legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 19),
    legend.text = element_text(size = 18),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()) +
  guides(fill = guide_legend(nrow = 7, title.position = "top", override.aes = list(size =3)))
dendroRA.genus.plot

#Dendrogram with RA plot
dendroRA.genus.plot.2 <- plot_grid(dendro.bray.plot, dendroRA.genus.plot,
                                   align = "v", 
                                   ncol = 1, 
                                   rel_heights = c(0.3, 0.7))
dendroRA.genus.plot.2


#DIFFERENTIAL ABUNDANCE####
##Feces (DNA vs cDNA) ##### 
###ANCOMBC#######
##Getting untransformed (raw) counts in fecal samples 
data.bacteria.samples.feces <- subset_samples(data.bacteria.samples, sample_type == "Feces") 
data.bacteria.samples.feces <- prune_taxa(taxa_sums(data.bacteria.samples.feces) > 0, data.bacteria.samples.feces) 
data.bacteria.samples.feces ## 36550 taxa and 96 samples
ancombc_feces.counts <-data.bacteria.samples.feces 

##GENUS
#Preprocessing, filtering out low (relative) abundance genera for feces samples
data.bacteria.samples.feces.ra <- subset_samples(data.bacteria.samples.tss, sample_type == "Feces") ##Only fecal samples (RA)
data.bacteria.samples.feces.ra <- prune_taxa(taxa_sums(data.bacteria.samples.feces.ra) > 0, data.bacteria.samples.feces.ra) 
data.bacteria.samples.feces.ra ##36550  taxa in fecal samples
data.bacteria.samples.feces.ra.genus <- tax_glom(data.bacteria.samples.feces.ra, taxrank = "Genus", NArm = F) ##Glom to the genus level
data.bacteria.samples.feces.ra.genus ##12453 genera in 96 fecal samples (RA)

##Filtering out the low relative abundance (less than 0.3 %) genera
data.bacteria.samples.feces_genus.ra.filt <- filter_taxa(data.bacteria.samples.feces.ra.genus, function(x) mean(x) > 0.3, TRUE) 
data.bacteria.samples.feces_genus.ra.filt ## 52 genera with mean RA > 0.3% across 96 samples (fecal samples) 
##Filtering those genera (> 0.3% RA) on the raw counts phyloseq object for feces
feces_genus.counts_filtered <- subset_taxa(ancombc_feces.counts, Genus %in% tax_table(data.bacteria.samples.feces_genus.ra.filt)[,"Genus"])
feces_genus.counts_filtered ##5294 taxa for those 50 genera (96 samples)
##To include feedlot as a random effect, making sure it is a factor (not continuous variable)
sample_data(feces_genus.counts_filtered)$feedlot <- factor(sample_data(feces_genus.counts_filtered)$feedlot, levels = c("1", "2", "3", "4", "5"))
##reorder gen_material as factor, DNA as "reference"
sample_data(feces_genus.counts_filtered)$gen_material <- factor(sample_data(feces_genus.counts_filtered)$gen_material, levels = c("DNA", "cDNA"))

##running ancombc on the variable of interest (gen_material)
ancombc_output_feces.genus <-ancombc2(data= feces_genus.counts_filtered, 
                                       assay_name = "counts", 
                                       tax_level = "Genus",
                                       fix_formula = "gen_material+feedlot",
                                       # fix_formula = "gen_material",
                                       # rand_formula =  "(1 | feedlot) + (1 | original_sample)", 
                                       rand_formula =  "(1 | original_sample)",
                                       prv_cut = 0.05, 
                                       lib_cut = 0, 
                                       group= "gen_material", 
                                       struc_zero = TRUE, 
                                       neg_lb = TRUE,
                                       alpha = 0.05, #default significance
                                       n_cl = 1, verbose = TRUE)
## extract results from comparisons 
res.feces.genus <- ancombc_output_feces.genus$res %>%
  select(-matches("feedlot"))

#Pivot longer the results
ancom_gen_material_feces.genus <- res.feces.genus %>%
  mutate(across(starts_with("lfc_gen_material"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancom_gen_material_feces.genus

##Getting rid of _rounded suffix using sub command
ancom_gen_material_feces.genus$group ##want to get rid of "_rounded"
ancom_gen_material_feces.genus$group<- sub("_rounded", "", ancom_gen_material_feces.genus$group) 
ancom_gen_material_feces.genus$group #names don't have "rounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_gen_material_feces.genus <- ancom_gen_material_feces.genus %>%
  mutate(group= case_when(
    group == "gen_materialcDNA" ~ "DNA vs cDNA",
    TRUE ~ group ##keeps original name for groups not specified (DNA)
  ))
ancom_gen_material_feces.genus$group ##Now the group names are shorter and more manageable
ancom_gen_material_feces.genus<- ancom_gen_material_feces.genus %>%
  rename(Genus = taxon) ##This ancombc was done at the genus level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancom_gen_material_feces.genus_2 <- ancom_gen_material_feces.genus %>%
  mutate(lower.ci = lfc_gen_materialcDNA - 1.96*se_gen_materialcDNA,
         upper.ci = lfc_gen_materialcDNA + 1.96*se_gen_materialcDNA)

##Final fix - up to make compatible with plotting
ancom_gen_material_feces.genus_3 <- ancom_gen_material_feces.genus_2 %>%
  filter (passed_ss_gen_materialcDNA == 1)%>%##Only want those that passed sensitivity testing
  #Will only include the classified taxa
  filter(!grepl("unknown|unclassified", Genus))%>%
  rename(stderr = se_gen_materialcDNA, ##Renaming
         pval = p_gen_materialcDNA, 
         qval = q_gen_materialcDNA)%>%
  mutate(direction = ifelse(value >0, "elevated", "depleted"),
         # Exponentiate the values first
         value_exp = exp(value),
         lower.ci_exp = exp(lower.ci),
         upper.ci_exp = exp(upper.ci)) %>%
  mutate(
    ##Getting log 2 so I can compare with maaslin
    coef =  log2(value_exp),
    lower.ci =  log2(lower.ci_exp),
    upper.ci =  log2(upper.ci_exp),
    plot = "Log2 Fold change with 95%CI", 
    test = "ANCOM-BC",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  # filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Genus, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
nrow(ancom_gen_material_feces.genus_3) ##26 DA genera between DNA and cDNA with ANCOM

###MaAsLin3#######
#Data (otu counts) and metadata for MaAslin
#Counts for data_maaslin_feces_ra are tss transformed
data.bacteria.samples.feces_genus.ra.filt #Will feed it this filtered ps object (filtered for those genera with mean RA > 0.3% across 96 samples (fecal samples))
data_maaslin_feces_ra  <- data.frame(t(data.bacteria.samples.feces_genus.ra.filt@otu_table)) 

##Sample metadata
metadata_maaslin_feces <- data.frame(data.bacteria.samples.feces_genus.ra.filt@sam_data) %>%
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")),
         gen_material = factor(gen_material, levels = c("DNA", "cDNA"))) #making feedlot a factor since I'll be adding it as a random effect on the MaAslin model

##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_feces_genus <-  maaslin3(
  input_data = data_maaslin_feces_ra, 
  input_metadata = metadata_maaslin_feces, 
  output = "MaAsLin3_feces",  
  fixed_effect = c("gen_material", 'feedlot'), 
  random_effects = c("original_sample"),
  # random_effects = c("feedlot", "original_sample"),
  min_prevalence=0.05,
  median_comparison_abundance = TRUE, #default
  median_comparison_prevalence = FALSE, #default  
  min_abundance = 0, ##input_data has already been filtered (>0.3% RA)
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_feces_genus$fit_data_abundance$results) #Abundance results from MaAslin - 260

##Taxonomy of feces OTUs 
input_taxonomy_feces <- data.frame(data.bacteria.samples.feces.ra.genus@tax_table) %>%
  rownames_to_column(var = "feature")

#Calculate confidence intervals, add taxonomy
maaslin_feces_genus_2 <- maaslin_feces_genus$fit_data_abundance$results %>%
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  left_join(input_taxonomy_feces, by = "feature")

##Final edits to put together for plot
maaslin_feces_genus_3 <- maaslin_feces_genus_2%>%
  #Will only include the classified taxa
  filter(!grepl("unknown|unclassified", Genus))%>%
  rename(group = name,
         pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(group= case_when(
    group == "gen_materialcDNA" ~ "DNA vs cDNA",
    TRUE ~ group ##keeps original name for groups not specified 
  ))%>%
  mutate(direction = ifelse(coef > 0, "elevated", "depleted"))%>%
  mutate(
    plot = "Log2 Fold change with 95%CI", 
    test = "MaAsLin3",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  # filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Genus, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancom_gen_material_feces.genus_3)
nrow(maaslin_feces_genus_3) ##21 DA genera between DNA and cDNA by MaAslin

##ANCOM and MaAslin together
DA_feces_plot_MaAslinANCOM.data <- rbind(ancom_gen_material_feces.genus_3, maaslin_feces_genus_3) %>%
  filter(Genus %in% intersect(maaslin_feces_genus_3$Genus,
                              ancom_gen_material_feces.genus_3$Genus)) ##Only going to plot those taxa DA by both tests

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_feces <- maaslin_feces_genus$transformed_data %>% #transformed data is TSS transformed otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "feature", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy_feces, by = "feature")%>%
  filter(Genus %in% intersect(maaslin_feces_genus_3$Genus,
                              ancom_gen_material_feces.genus_3$Genus))%>%
  left_join(metadata_maaslin_feces%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(Genus, Sample, logvalue, gen_material, plot)

#Bias-corrected abundances (ANCOM)
feces.genus_log_corr_abn <- ancombc_output_feces.genus$bias_correct_log_table %>%
  data.frame()%>% ##make into data frame
  rownames_to_column("Genus")%>%
  filter(Genus %in% intersect(maaslin_feces_genus_3$Genus,
                              ancom_gen_material_feces.genus_3$Genus))%>% #keep only those genera in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         gen_material = ifelse(grepl("c", Sample), "cDNA", "DNA"))

#put together objects to plot DA
DA_feces_plot_together <- bind_rows(DA_feces_plot_MaAslinANCOM.data, 
                                    feces.genus_log_corr_abn, 
                                    RA_MaaslinAncom_feces) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_feces_plot_together$plot <- factor(DA_feces_plot_together$plot, 
                                      levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_feces_plot_together$gen_material <- factor(DA_feces_plot_together$gen_material, levels = c("DNA", "cDNA"))


###PLOTTING DA#
##Ordering how I want the "Genus" taxlevel to show up on the plot - according to family
input_taxonomy_feces ##dataframe object for Taxonomy of feces OTUs 

# Create the taxonomy plot data and modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_fecal <- DA_feces_plot_together %>%
  left_join(input_taxonomy_feces, by= "Genus")%>%
  distinct(Phylum, Order, Family, Genus) %>%
  arrange(Phylum, Order, Family) %>%
  mutate(Genus = factor(Genus, levels = rev(Genus)))%>% ##Since I arranged by family, this is the order I want the genera to show up
  dplyr::group_by(Family) %>%
  dplyr::mutate(label_Family = ifelse(row_number() == 1, Family, "")) %>%  # Create 'label_Family' with only the first occurrence of each Family
  ungroup()
levels(taxonomy_plot_data_fecal$Genus)

##Factor "Genus" level by the order I want (taxonomy_plot_data_fecal$Genus)
DA_feces_plot_together$Genus <- factor(DA_feces_plot_together$Genus, levels = rev(taxonomy_plot_data_fecal$Genus))


# Create the updated taxonomy plot
taxonomy_plot_feces <- ggplot(taxonomy_plot_data_fecal) +
  geom_text(aes(x =0, y = Genus, label = label_Family), hjust = 0, size = 8, family = "sans") +  # Move text left by adjusting x
  labs(title = "Family") +
  theme_void() +
  #scale_y_discrete(limits =  rev(taxonomy_plot_data_DNA$Genus)) +  # Ensure y-axis matches the taxon order
  theme(plot.title = element_text(hjust = 0, size = 22, vjust = -0.5, face = "bold", family = "sans"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )+
  coord_cartesian(xlim = c(0.05, 1)) #had to add this to move the geom_text more to the left
taxonomy_plot_feces

#Plotting 
DA_feces_plot_MaAslinANCOM <-
  ggplot(data=DA_feces_plot_together%>%filter(grepl("abundances", plot)),
         aes(x=Genus, y=logvalue, fill = gen_material, color = gen_material)) +
  geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order = 1,title="Library Type", title.position="top"))+
  scale_color_manual(values = gen.material.palette,
                     labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
  scale_fill_manual(values=gen.material.palette,
                    labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
  new_scale_color()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_feces_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=Genus, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_feces_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=Genus, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_feces_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed'
             , alpha=0.5) +
  facet_nested(. ~ plot, scales='free_x',
               space='free_y',
               switch='y',
               strip=strip_nested(text_y=list(element_text(angle=0))),
               labeller=labeller(group=label_wrap_gen(width=10),
                                 sub_group=label_wrap_gen(width=10))) +
  scale_x_discrete(position='bottom') +
  scale_y_continuous(position='right') +  # Default to break_labels otherwise 
  coord_flip() +
  scale_color_manual(values=c("red", "blue")) +
  scale_linetype_manual(values=c("11", "solid")) +
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Library Type", title.position="top"),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 2.5)),
         linetype = guide_legend(title = "Fold change source", title.position = "top",
                                 override.aes = list(linewidth = 1),
                                 theme = theme(legend.key.width = unit(1.5, "cm"))),
         pch=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  theme_bw()+
  labs(title = "FECES")+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=26), legend.text=element_text(size=26),
        plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title.x=element_blank(), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size= 22, angle=0, vjust= 1.03, face = "bold"), 
        axis.text.y=element_text(size=20, vjust = 0.5),
        strip.text=element_text(size=16, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_feces_plot_MaAslinANCOM

##Adding the q values
DA_feces_plot_MaAslinANCOM_q <- DA_feces_plot_MaAslinANCOM +
  geom_text(inherit.aes=FALSE,
            data = DA_feces_plot_together %>% filter(plot == 'Log2 Fold change with 95%CI'),
            aes(x=Genus, y = coef, label = DA),
            position = position_dodge2(width = 0.75),
            vjust = -0.6, size = 3)
DA_feces_plot_MaAslinANCOM_q

##Putting together DA plot with taxonomy (at the family level) plot
combined_plot_feces <- plot_grid(
  taxonomy_plot_feces+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_feces_plot_MaAslinANCOM  + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                    plot.title = element_blank(),
                                    strip.text=element_text(size=16, color = "white", face = "bold")),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "FECES")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_feces

##Adding q values 
combined_plot_feces_q <- plot_grid(
  taxonomy_plot_feces+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_feces_plot_MaAslinANCOM_q + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                     plot.title = element_blank(),
                                     strip.text=element_text(size=16, color = "white", face = "bold")),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "FECES")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_feces_q

##Catch Basins (DNA vs cDNA) #######
###ANCOMBC#######
##Getting untransformed (raw) counts in water samples 
data.bacteria.samples.water <- subset_samples(data.bacteria.samples, sample_type == "Water") 
data.bacteria.samples.water <- prune_taxa(taxa_sums(data.bacteria.samples.water) > 0, data.bacteria.samples.water) 
data.bacteria.samples.water ## 32601 taxa and 24 samples
ancombc_water.counts <-data.bacteria.samples.water 
ancombc_water.counts@sam_data$gen_material <- factor(ancombc_water.counts@sam_data$gen_material, levels = c("DNA", "cDNA"))##reorder gen_material as factor, DNA as "reference"

##GENUS
#Preprocessing, filtering out low (relative) abundance genera for water samples
data.bacteria.samples.water.ra <- subset_samples(data.bacteria.samples.tss, sample_type == "Water") ##Only water samples (RA)
data.bacteria.samples.water.ra <- prune_taxa(taxa_sums(data.bacteria.samples.water.ra) > 0, data.bacteria.samples.water.ra) 
data.bacteria.samples.water.ra ##32601 taxa in water samples (24)
data.bacteria.samples.water.ra.genus <- tax_glom(data.bacteria.samples.water.ra, taxrank = "Genus", NArm = F) ##Glom to the genus level
data.bacteria.samples.water.ra.genus ##11590 genera in 24 water samples


##Filtering out the low relative abundance (less than 0.3 %) genera
data.bacteria.samples.water_genus.ra.filt <- filter_taxa(data.bacteria.samples.water.ra.genus, function(x) mean(x) > 0.3, TRUE) 
data.bacteria.samples.water_genus.ra.filt ## 51 genera with mean RA > 0.3% across 24 samples (water samples) 
##Filtering those genera (> 0.3% RA) on the raw counts phyloseq object for water
water_genus.counts_filtered <- subset_taxa(ancombc_water.counts, Genus %in% tax_table(data.bacteria.samples.water_genus.ra.filt)[,"Genus"])
water_genus.counts_filtered ##4894 taxa for those 51 genera (24 samples)
##To include feedlot as a random effect, making sure it is a factor (not continuous variable)
sample_data(water_genus.counts_filtered)$feedlot <- factor(sample_data(water_genus.counts_filtered)$feedlot, 
                                                           levels = c("1","2","3","4","5"))

##running ancombc on the variable of interest (gen_material)
ancombc_output_water.genus <-ancombc2(data= water_genus.counts_filtered, 
                                      assay_name = "counts", 
                                      tax_level = "Genus",
                                      fix_formula = "gen_material+feedlot",
                                      # fix_formula = "gen_material",
                                      # rand_formula =  "(1 | feedlot) + (1 | original_sample)", 
                                      rand_formula =  "(1 | original_sample)",
                                      prv_cut = 0.05, 
                                      lib_cut = 0, 
                                      group= "gen_material", 
                                      struc_zero = TRUE, 
                                      neg_lb = TRUE,
                                      alpha = 0.05, #default significance
                                      n_cl = 1, verbose = TRUE)

## extract results from comparisons 
res.water.genus <- ancombc_output_water.genus$res %>%
  select(-matches("feedlot")) #not doing feedlot comparisons


#Pivot longer the results
ancom_gen_material_water.genus <- res.water.genus %>%
  mutate(across(starts_with("lfc_gen_material"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancom_gen_material_water.genus

##Getting rid of _rounded suffix using sub command
ancom_gen_material_water.genus$group ##want to get rid of "_rounded"
ancom_gen_material_water.genus$group<- sub("_rounded", "", ancom_gen_material_water.genus$group) 
ancom_gen_material_water.genus$group #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_gen_material_water.genus <- ancom_gen_material_water.genus %>%
  mutate(group= case_when(
    group == "gen_materialcDNA" ~ "DNA vs cDNA",
    TRUE ~ group ##keeps original name for groups not specified (DNA)
  ))
ancom_gen_material_water.genus$group ##Now the group names are shorter and more manageable
ancom_gen_material_water.genus<- ancom_gen_material_water.genus %>%
  rename(Genus = taxon) ##This ancombc was done at the genus level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancom_gen_material_water.genus_2 <- ancom_gen_material_water.genus %>%
  mutate(lower.ci = lfc_gen_materialcDNA - 1.96*se_gen_materialcDNA,
         upper.ci = lfc_gen_materialcDNA + 1.96*se_gen_materialcDNA)

##Final fix - up to make compatible with plotting
ancom_gen_material_water.genus_3 <- ancom_gen_material_water.genus_2 %>%
  #filter (q_gen_materialcDNA < 0.05 & passed_ss_gen_materialcDNA == 1) %>% ##Only want significant differences that passed sensitivity testing
  filter (passed_ss_gen_materialcDNA == 1)%>% #only want those that passed sensitivity testing
  #Will only include the classified taxa
  filter(!grepl("unknown|unclassified", Genus))%>%
  rename(stderr = se_gen_materialcDNA, ##Renaming
         pval = p_gen_materialcDNA, 
         qval = q_gen_materialcDNA)%>%
  mutate(direction = ifelse(value >0, "elevated", "depleted"),
         # Exponentiate the values first
         value_exp = exp(value),
         lower.ci_exp = exp(lower.ci),
         upper.ci_exp = exp(upper.ci)) %>%
  mutate(
    ##Getting log 2 so I can compare with maaslin
    coef =  log2(value_exp),
    lower.ci =  log2(lower.ci_exp),
    upper.ci =  log2(upper.ci_exp),
    
    plot = "Log2 Fold change with 95%CI", 
    test = "ANCOM-BC",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  # filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Genus, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
nrow(ancom_gen_material_water.genus_3) ##7 DA genera between DNA and cDNA with ANCOM

###MaAsLin3#######
#Data (otu counts) and metadata for MaAslin
#Counts for data_maaslin_water_ra are tss transformed
data.bacteria.samples.water_genus.ra.filt #Will feed it this filtered ps object (filtered for those genera with mean RA > 0.3% across 96 samples (water samples))
data_maaslin_water_ra  <- data.frame(t(data.bacteria.samples.water_genus.ra.filt@otu_table)) 

##Sample metadata
metadata_maaslin_water <- data.frame(data.bacteria.samples.water_genus.ra.filt@sam_data) %>%
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")),
         gen_material = factor(gen_material, levels = c("DNA", "cDNA"))) #making feedlot a factor since I'll be adding it as a random effect on the MaAslin model

##Will also try out with raw counts, glommed at genus level but not filtered (so when MaAslin does tss normalization, it takes into account all OTUs)
data.bacteria.samples.water.genus <- tax_glom(data.bacteria.samples.water, taxrank = "Genus", NArm = F) #data.bacteria.samples.water has raw counts
data_maaslin_water  <- data.frame(t(data.bacteria.samples.water.genus@otu_table))


##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_water_genus <-  maaslin3(
  input_data = data_maaslin_water_ra, 
  input_metadata = metadata_maaslin_water, 
  output = "MaAsLin3_water",  
  fixed_effect = c("gen_material", 'feedlot'), 
  random_effects = c("original_sample"),
  # fixed_effect = "gen_material",
  # random_effects = c("feedlot", "original_sample"),
  min_prevalence=0.05,
  median_comparison_abundance = T, #default
  median_comparison_prevalence = FALSE, #default  
  min_abundance = 0, ##input_data has already been filtered (>0.3% RA)
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_water_genus$fit_data_abundance$results) #Abundance results from MaAslin - 255


##Taxonomy of water OTUs 
input_taxonomy_water <- data.frame(data.bacteria.samples.water.ra.genus@tax_table) %>%
  rownames_to_column(var = "feature")

#Calculate confidence intervals, add taxonomy
maaslin_water_genus_2 <- maaslin_water_genus$fit_data_abundance$results %>%
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) %>%
  filter(metadata != "feedlot")%>%
  left_join(input_taxonomy_water, by = "feature")

##Final edits to put together for plot
maaslin_water_genus_3 <- maaslin_water_genus_2%>%
  #Will only include the classified taxa
  filter(!grepl("unknown|unclassified", Genus))%>%
  rename(group = name,
         pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(group= case_when(
    group == "gen_materialcDNA" ~ "DNA vs cDNA",
    TRUE ~ group ##keeps original name for groups not specified 
  ))%>%
  mutate(direction = ifelse(coef > 0, "elevated", "depleted"))%>%
  mutate(
    plot = "Log2 Fold change with 95%CI", 
    test = "MaAsLin3",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  # filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Genus, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns (same as ancom_gen_material_feces.genus_3)
nrow(maaslin_water_genus_3) ##8 DA genera between DNA and cDNA by MaAslin

##ANCOM and MaAslin together
DA_water_plot_MaAslinANCOM.data <- rbind(ancom_gen_material_water.genus_3, maaslin_water_genus_3) %>%
  filter(Genus %in% intersect(maaslin_water_genus_3$Genus,
                              ancom_gen_material_water.genus_3$Genus)) ##Only going to plot those taxa DA by both tests

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_water <- maaslin_water_genus$transformed_data %>% #transformed data is TSS transformed otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "feature", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy_water, by = "feature")%>%
  filter(Genus %in% intersect(maaslin_water_genus_3$Genus,
                              ancom_gen_material_water.genus_3$Genus))%>%
  left_join(metadata_maaslin_water%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(Genus, Sample, logvalue, gen_material, plot)


#Bias-corrected abundances (ANCOM)
water.genus_log_corr_abn <- ancombc_output_water.genus$bias_correct_log_table %>%
  data.frame()%>% ##make into data frame
  rownames_to_column("Genus")%>%
  filter(Genus %in% intersect(maaslin_water_genus_3$Genus,
                              ancom_gen_material_water.genus_3$Genus))%>% #keep only those genera in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         gen_material = ifelse(grepl("c", Sample), "cDNA", "DNA"))

#put together objects to plot DA
DA_water_plot_together <- bind_rows(DA_water_plot_MaAslinANCOM.data, water.genus_log_corr_abn, RA_MaaslinAncom_water) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_water_plot_together$plot <- factor(DA_water_plot_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_water_plot_together$gen_material <- factor(DA_water_plot_together$gen_material, levels = c("DNA", "cDNA"))


###PLOTTING DA#
##Ordering how I want the "Genus" taxlevel to show up on the plot 
input_taxonomy_water ##dataframe object for Taxonomy of water OTUs 

##IF DOING FAMILY:
# Create the taxonomy plot data and modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_water <- DA_water_plot_together %>%
  left_join(input_taxonomy_water, by= "Genus")%>%
  distinct(Order, Family, Genus) %>%
  arrange(Order, Family) %>%
  mutate(Genus = factor(Genus, levels = rev(Genus)))%>% ##Since I arranged by family, this is the order I want the genera to show up
  dplyr::group_by(Family) %>%
  dplyr::mutate(label_Family = ifelse(row_number() == 1, Family, "")) %>%  # Create 'label_Family' with only the first occurrence of each family
  ungroup()
levels(taxonomy_plot_data_water$Genus)


##Factor "Genus" level by the order I want (taxonomy_plot_data_water$Genus)
DA_water_plot_together$Genus <- factor(DA_water_plot_together$Genus, levels = rev(taxonomy_plot_data_water$Genus))

# Create the updated taxonomy plot
taxonomy_plot_water <- ggplot(taxonomy_plot_data_water) +
  geom_text(aes(x =0, y = Genus, label = label_Family), hjust = 0, size = 8, family = "sans") +  # Move text left by adjusting x
  labs(title = "Family") +
  theme_void() +
  #scale_y_discrete(limits =  rev(taxonomy_plot_data_DNA$Genus)) +  # Ensure y-axis matches the taxon order
  theme(plot.title = element_text(hjust = 0, size = 22, vjust = -0.5, face = "bold", family = "sans"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )+
  coord_cartesian(xlim = c(0.05, 1)) #had to add this to move the geom_text more to the left


#Plotting 
DA_water_plot_MaAslinANCOM <-
  ggplot(data=DA_water_plot_together%>%filter(grepl("abundances", plot)),
         aes(x=Genus, y=logvalue, fill = gen_material, color = gen_material)) +
  geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order = 1,title="Library Type", title.position="top"))+
  scale_color_manual(values = gen.material.palette,
                     labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
  scale_fill_manual(values=gen.material.palette,
                    labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
  new_scale_color()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_water_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=Genus, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_water_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=Genus, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_water_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed'
             , alpha=0.5) +
  facet_nested(. ~ plot, scales='free_x',
               space='free_y',
               switch='y',
               strip=strip_nested(text_y=list(element_text(angle=0))),
               labeller=labeller(group=label_wrap_gen(width=10),
                                 sub_group=label_wrap_gen(width=10))) +
  scale_x_discrete(position='bottom') +
  scale_y_continuous(position='right') +  # Default to break_labels otherwise 
  coord_flip() +
  scale_color_manual(values=c("red", "blue")) +
  scale_linetype_manual(values=c("11", "solid")) +
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Library Type", title.position="top"),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 2.5)),
         linetype = guide_legend(title = "Fold change source", title.position = "top",
                                 override.aes = list(linewidth = 1),
                                 theme = theme(legend.key.width = unit(1.5, "cm"))),
         pch=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  labs(title = "CATCH BASINS")+
  theme_bw()+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=26), legend.text=element_text(size=26),
        plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y=element_text(size= 22, angle=0, vjust= 1.06, face = "bold"), 
        axis.text.y=element_text(size=20, vjust = 0.5),
        strip.text=element_text(size=16, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_water_plot_MaAslinANCOM

##Adding the q values
DA_water_plot_MaAslinANCOM_q <- DA_water_plot_MaAslinANCOM +
  geom_text(inherit.aes=FALSE,
            data = DA_water_plot_together %>% filter(plot == 'Log2 Fold change with 95%CI'),
            aes(x = Genus, y = coef, label = DA),
            position = position_dodge2(width = 0.75),
            vjust = -0.6, size = 3)

##Putting together DA plot with taxonomy (family level) plot
combined_plot_water <- plot_grid(
  taxonomy_plot_water+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_water_plot_MaAslinANCOM  + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                    plot.title = element_blank(),
                                    # strip.text = element_blank(),
                                    # strip.background = element_rect(fill = "white"),
                                    legend.position = "none"
                                    ),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "CATCH BASINS")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_water

##Adding q values 
combined_plot_water_q <- plot_grid(
  taxonomy_plot_water+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_water_plot_MaAslinANCOM_q + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                     plot.title = element_blank(),
                                     strip.text=element_text(size=16, color = "white", face = "bold")),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "CATCH BASINS")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_water_q

####FIGURE 8 A & B - DA EMPHASIS ON DNA vs cDNA#####
#Put together these 2 plots:
combined_plot_feces
combined_plot_water
figure8AB <-plot_grid(combined_plot_feces+
                         theme(plot.title = element_blank()), 
                    combined_plot_water+
                      theme(plot.title = element_blank()),
                    align = "v",
                    labels = c("A", "B"),
                    label_size = 32,
                    ncol = 1,
                    rel_heights = c(0.75, 0.25))
  # labs(title = "MICROBIOME")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
figure8AB
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure8AB.png", 
       plot = figure8AB, 
       dpi = 600,
       device = "png",
       width = 19, 
       height = 15, 
       bg = "white")

##DNA (CB vs Feces)#####
###ANCOMBC#######
##Getting untransformed (raw) counts in DNA samples 
data.bacteria.samples.DNA <- subset_samples(data.bacteria.samples, gen_material == "DNA") 
data.bacteria.samples.DNA <- prune_taxa(taxa_sums(data.bacteria.samples.DNA) > 0, data.bacteria.samples.DNA) 
data.bacteria.samples.DNA ## 39081 taxa and 60 samples
ancombc_DNA.counts <-data.bacteria.samples.DNA 
ancombc_DNA.counts@sam_data$sample_type <- factor(ancombc_DNA.counts@sam_data$sample_type, levels = c("Feces", "Water"))##reorder sample_type as factor, Feces as "reference"

##GENUS
#Preprocessing, filtering out low (relative) abundance genera for DNA samples
data.bacteria.samples.DNA.ra <- subset_samples(data.bacteria.samples.tss, gen_material == "DNA") ##Only DNA samples (RA)
data.bacteria.samples.DNA.ra <- prune_taxa(taxa_sums(data.bacteria.samples.DNA.ra) > 0, data.bacteria.samples.DNA.ra) 
data.bacteria.samples.DNA.ra ##39081 taxa and 60 DNA samples (RA)
data.bacteria.samples.DNA.ra.genus <- tax_glom(data.bacteria.samples.DNA.ra, taxrank = "Genus", NArm = F) ##Glom to the genus level
data.bacteria.samples.DNA.ra.genus ##12416 genera in 60 DNA samples


##Filtering out the low relative abundance (less than 0.3 %) genera
data.bacteria.samples.DNA_genus.ra.filt <- filter_taxa(data.bacteria.samples.DNA.ra.genus, function(x) mean(x) > 0.3, TRUE) 
data.bacteria.samples.DNA_genus.ra.filt ## 53 genera with mean RA > 0.3% across 60 samples (DNA samples) 
##Filtering those genera (> 0.3% RA) on the raw counts phyloseq object for DNA
DNA_genus.counts_filtered <- subset_taxa(ancombc_DNA.counts, Genus %in% tax_table(data.bacteria.samples.DNA_genus.ra.filt)[,"Genus"])
DNA_genus.counts_filtered ##5014 taxa for those 53 genera (60 samples)
##To include feedlot as a random effect, making sure it is a factor (not continuous variable)
sample_data(DNA_genus.counts_filtered)$feedlot <- factor(sample_data(DNA_genus.counts_filtered)$feedlot, levels = c("1", "2", "3", "4", "5"))

##running ancombc on the variable of interest (sample_type)
ancombc_output_DNA.genus <-ancombc2(data= DNA_genus.counts_filtered, 
                                      assay_name = "counts", 
                                      tax_level = "Genus",
                                      fix_formula = "sample_type + feedlot", 
                                      # rand_formula =  "(1 | feedlot)", 
                                      prv_cut = 0.05, 
                                      lib_cut = 0, 
                                      group= "sample_type", 
                                      struc_zero = TRUE, 
                                      neg_lb = TRUE,
                                      alpha = 0.05, #default significance
                                      n_cl = 1, verbose = TRUE)

## extract results from comparisons 
res.DNA.genus <- ancombc_output_DNA.genus$res %>%
  select(-matches("feedlot")) #not doing feedlot comparisons

#Pivot longer the results
ancom_sample_type_DNA.genus <- res.DNA.genus %>%
  mutate(across(starts_with("lfc_sample_type"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancom_sample_type_DNA.genus

##Getting rid of _rounded suffix using sub command
ancom_sample_type_DNA.genus$group ##want to get rid of "_rounded"
ancom_sample_type_DNA.genus$group<- sub("_rounded", "", ancom_sample_type_DNA.genus$group) 
ancom_sample_type_DNA.genus$group #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_sample_type_DNA.genus <- ancom_sample_type_DNA.genus %>%
  mutate(group= case_when(
    group == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ group ##keeps original name for groups not specified (DNA)
  ))
ancom_sample_type_DNA.genus$group ##Now the group names are shorter and more manageable
ancom_sample_type_DNA.genus<- ancom_sample_type_DNA.genus %>%
  rename(Genus = taxon) ##This ancombc was done at the genus level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancom_sample_type_DNA.genus_2 <- ancom_sample_type_DNA.genus %>%
  mutate(lower.ci = lfc_sample_typeWater - 1.96*se_sample_typeWater,
         upper.ci = lfc_sample_typeWater + 1.96*se_sample_typeWater)

##Final fix - up to make compatible with plotting
ancom_sample_type_DNA.genus_3 <- ancom_sample_type_DNA.genus_2 %>%
  filter (passed_ss_sample_typeWater == 1) %>% ##Only want those that passed sensitivity testing
  #Will only include the classified taxa
  filter(!grepl("unknown|unclassified", Genus))%>%
  rename(stderr = se_sample_typeWater, ##Renaming
         pval = p_sample_typeWater, 
         qval = q_sample_typeWater)%>%
  mutate(direction = ifelse(value >0, "elevated", "depleted"),
    # Exponentiate the values first
    value_exp = exp(value),
    lower.ci_exp = exp(lower.ci),
    upper.ci_exp = exp(upper.ci)) %>%
  mutate(
    ##Getting log 2 so I can compare with maaslin
    coef =  log2(value_exp),
    lower.ci =  log2(lower.ci_exp),
    upper.ci =  log2(upper.ci_exp),
    plot = "Log2 Fold change with 95%CI", 
    test = "ANCOM-BC",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  # filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Genus, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns
nrow(ancom_sample_type_DNA.genus_3) ##25 DA genera between Feces and Water with ANCOM

###MaAsLin3#######
#Data (otu counts) and metadata for MaAslin
#Counts for data_maaslin_DNA_ra are tss transformed
data.bacteria.samples.DNA_genus.ra.filt #Will feed it this filtered ps object (filtered for those genera with mean RA > 0.3% across 60 samples (DNA samples))
data_maaslin_DNA_ra  <- data.frame(t(data.bacteria.samples.DNA_genus.ra.filt@otu_table)) 

##Sample metadata
metadata_maaslin_DNA <- data.frame(data.bacteria.samples.DNA_genus.ra.filt@sam_data) %>%
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")),
         sample_type = factor(sample_type, levels = c("Feces", "Water"))) #making feedlot a factor since I'll be adding it as a random effect on the MaAslin model

##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_DNA_genus <-  maaslin3(
  input_data = data_maaslin_DNA_ra, 
  input_metadata = metadata_maaslin_DNA, 
  output = "MaAsLin3_DNA",  
  fixed_effect = c("sample_type", "feedlot"),
  # fixed_effect = "sample_type",
  # random_effects = "feedlot",
  min_prevalence=0.05,
  median_comparison_abundance = T, #default
  median_comparison_prevalence = FALSE, #default  
  min_abundance = 0, ##input_data has already been filtered (>0.3% RA)
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_DNA_genus$fit_data_abundance$results) #Abundance results from MaAslin - 265

##Taxonomy of DNA OTUs 
input_taxonomy_DNA <- data.frame(data.bacteria.samples.DNA.ra.genus@tax_table) %>%
  rownames_to_column(var = "feature")

#Calculate confidence intervals, add taxonomy
maaslin_DNA_genus_2 <- maaslin_DNA_genus$fit_data_abundance$results %>%
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  left_join(input_taxonomy_DNA, by = "feature")

##Final edits to put together for plot
maaslin_DNA_genus_3 <- maaslin_DNA_genus_2%>%
  #Will only include the classified taxa
  filter(!grepl("unknown|unclassified", Genus))%>%
  rename(group = name,
         pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(group= case_when(
    group == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ group ##keeps original name for groups not specified 
  ))%>%
  mutate(direction = ifelse(coef > 0, "elevated", "depleted"))%>%
  mutate(
    
    plot = "Log2 Fold change with 95%CI", 
    test = "MaAsLin3", 
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  # filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Genus, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancom_sample_type_DNA.genus_3)
nrow(maaslin_DNA_genus_3) ##23 DA genera between Feces and Water by MaAslin

##ANCOM and MaAslin together
DA_DNA_plot_MaAslinANCOM.data <- rbind(ancom_sample_type_DNA.genus_3, maaslin_DNA_genus_3) %>%
  ##Only going to plot those taxa DA by both tests
  filter(Genus %in% intersect(maaslin_DNA_genus_3$Genus,
                              ancom_sample_type_DNA.genus_3$Genus))

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_DNA <- maaslin_DNA_genus$transformed_data %>% #transformed data is TSS (ra) log2 transformed otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "feature", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy_DNA, by = "feature")%>%
  filter(Genus %in% intersect(maaslin_DNA_genus_3$Genus,
                              ancom_sample_type_DNA.genus_3$Genus))%>%
  left_join(metadata_maaslin_DNA%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(Genus, Sample, logvalue, sample_type, plot)


#Bias-corrected abundances (ANCOM)
DNA.genus_log_corr_abn <- ancombc_output_DNA.genus$bias_correct_log_table %>%
  data.frame()%>% ##make into data frame
  rownames_to_column("Genus")%>%
  filter(Genus %in% intersect(maaslin_DNA_genus_3$Genus,
                              ancom_sample_type_DNA.genus_3$Genus))%>% #keep only those genera in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         sample_type = ifelse(grepl("W", Sample), "Water", "Feces"))

#put together objects to plot DA
DA_DNA_plot_together <- bind_rows(DA_DNA_plot_MaAslinANCOM.data, 
                                  DNA.genus_log_corr_abn, 
                                  RA_MaaslinAncom_DNA) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_DNA_plot_together$plot <- factor(DA_DNA_plot_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_DNA_plot_together$sample_type <- factor(DA_DNA_plot_together$sample_type, levels = c("Feces", "Water"))


###PLOTTING DA#
##Ordering how I want the "Genus" taxlevel to show up on the plot 
input_taxonomy_DNA ##dataframe object for Taxonomy of DNA OTUs 

##PHYLUM LEVEL 
# Create the taxonomy plot data and modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_DNA <- DA_DNA_plot_together %>%
  left_join(input_taxonomy_DNA, by= "Genus")%>%
  distinct(Phylum, Order, Family, Genus) %>%
  arrange(Phylum, Order, Family) %>%
  mutate(Genus = factor(Genus, levels = rev(Genus)))%>% ##Since I arranged by family, this is the order I want the genera to show up
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Family' with only the first occurrence of each class
  ungroup() 
levels(taxonomy_plot_data_DNA$Genus)

##Factor "Genus" level by the order I want (taxonomy_plot_data_DNA$Genus)
DA_DNA_plot_together$Genus <- factor(DA_DNA_plot_together$Genus, levels = rev(taxonomy_plot_data_DNA$Genus))

# Create the updated taxonomy plot
taxonomy_plot_DNA <- ggplot(taxonomy_plot_data_DNA) +
  geom_text(aes(x =0, y = Genus, label = label_Phylum), hjust = 0, size = 8, family = "sans") +  # Move text left by adjusting x
  labs(title = "Phylum") +
  theme_void() +
  #scale_y_discrete(limits =  rev(taxonomy_plot_data_DNA$Genus)) +  # Ensure y-axis matches the taxon order
  theme(plot.title = element_text(hjust = 0, size = 22, vjust = -0.5, face = "bold", family = "sans"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )+
  coord_cartesian(xlim = c(0.05, 1)) #had to add this to move the geom_text more to the left
taxonomy_plot_DNA

#Plotting 
DA_DNA_plot_MaAslinANCOM <-
  ggplot(data=DA_DNA_plot_together%>%filter(grepl("abundances", plot)),
         aes(x=Genus, y=logvalue, fill = sample_type, color = sample_type)) +
  geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order = 1,title="Sample Type", title.position="top"),
         fill=guide_legend(order = 1,title="Sample Type", title.position="top"))+
  scale_fill_manual(values=sample.type.palette, 
                    labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  scale_color_manual(values = sample.type.palette, 
                     labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  new_scale_color()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_DNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=Genus, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_DNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=Genus, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_DNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed'
             , alpha=0.5) +
  facet_nested(. ~ plot, scales='free_x',
               space='free_y',
               switch='y',
               strip=strip_nested(text_y=list(element_text(angle=0))),
               labeller=labeller(group=label_wrap_gen(width=10),
                                 sub_group=label_wrap_gen(width=10))) +
  scale_x_discrete(position='bottom') +
  scale_y_continuous(position='right') +  # Default to break_labels otherwise 
  coord_flip() +
  scale_color_manual(values=c("red", "blue")) +
  scale_linetype_manual(values=c("11", "solid")) +
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Sample Type", title.position="top"),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 2.5)),
         linetype = guide_legend(title = "Fold change source", title.position = "top",
                                 override.aes = list(linewidth = 1),
                                 theme = theme(legend.key.width = unit(1.5, "cm"))),
         shape=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  labs(title = "Metagenomic libraries (DNA)")+
  theme_bw()+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=26), legend.text=element_text(size=26),
        plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size= 22, angle=0, vjust= 1.045, face = "bold"), 
        axis.text.y=element_text(size=20, vjust = 0.5),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text=element_text(size=16, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_DNA_plot_MaAslinANCOM

##Adding the q values
DA_DNA_plot_MaAslinANCOM_q <- DA_DNA_plot_MaAslinANCOM +
  geom_text(inherit.aes=FALSE,
            data = DA_DNA_plot_together %>% filter(plot == 'Log2 Fold change with 95%CI'),
            aes(x = Genus, y = coef, label = DA),
            position = position_dodge2(width = 0.75),
            vjust = -0.6, size = 3)

##Putting together DA plot with taxonomy (family level) plot
combined_plot_DNA <- plot_grid(
  taxonomy_plot_DNA+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_DNA_plot_MaAslinANCOM  + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                    plot.title = element_blank(),
                                    strip.text=element_text(size=16, color = "white", face = "bold")),
  ncol = 2, 
  rel_widths = c(0.3,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "Metagenomic libraries (DNA)")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_DNA

##Adding q values 
combined_plot_DNA_q <- plot_grid(
  taxonomy_plot_DNA+ 
    theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_DNA_plot_MaAslinANCOM_q + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                    plot.title = element_blank(),
                                    strip.text=element_text(size=16, color = "white", face = "bold")),
  ncol = 2, 
  rel_widths = c(0.3,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "Metagenomic libraries (DNA)")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_DNA_q

##cDNA (CB vs Feces)#####
###ANCOMBC#######
##Getting untransformed (raw) counts in fecal samples 
data.bacteria.samples.cDNA <- subset_samples(data.bacteria.samples, gen_material == "cDNA") 
data.bacteria.samples.cDNA <- prune_taxa(taxa_sums(data.bacteria.samples.cDNA) > 0, data.bacteria.samples.cDNA) 
data.bacteria.samples.cDNA ## 33954 taxa and 60 samples
ancombc_cDNA.counts <-data.bacteria.samples.cDNA 
ancombc_cDNA.counts@sam_data$sample_type <- factor(ancombc_cDNA.counts@sam_data$sample_type, levels = c("Feces", "Water"))##reorder sample_type as factor, Feces as "reference"

##GENUS
#Preprocessing, filtering out low (relative) abundance genera for cDNA samples
data.bacteria.samples.cDNA.ra <- subset_samples(data.bacteria.samples.tss, gen_material == "cDNA") ##Only cDNA samples (RA)
data.bacteria.samples.cDNA.ra <- prune_taxa(taxa_sums(data.bacteria.samples.cDNA.ra) > 0, data.bacteria.samples.cDNA.ra) 
data.bacteria.samples.cDNA.ra ##33954 taxa in 60 cDNA samples (RA)
data.bacteria.samples.cDNA.ra.genus <- tax_glom(data.bacteria.samples.cDNA.ra, taxrank = "Genus", NArm = F) ##Glom to the genus level
data.bacteria.samples.cDNA.ra.genus ##11262  genera in 60 samples


##Filtering out the low relative abundance (less than 0.3 %) genera
data.bacteria.samples.cDNA_genus.ra.filt <- filter_taxa(data.bacteria.samples.cDNA.ra.genus, function(x) mean(x) > 0.3, TRUE) 
data.bacteria.samples.cDNA_genus.ra.filt ## 53 genera with mean RA > 0.3% across 60 samples (cDNA samples) 
##Filtering those genera (> 0.3% RA) on the raw counts phyloseq object for cDNA
cDNA_genus.counts_filtered <- subset_taxa(ancombc_cDNA.counts, Genus %in% tax_table(data.bacteria.samples.cDNA_genus.ra.filt)[,"Genus"])
cDNA_genus.counts_filtered ##3886 taxa for those 51 genera (60 cDNA samples)
##To include feedlot as a random effect, making sure it is a factor (not continuous variable)
sample_data(cDNA_genus.counts_filtered)$feedlot <- factor(sample_data(cDNA_genus.counts_filtered)$feedlot, levels = c("1", "2", "3", "4", "5"))


##running ancombc on the variable of interest (sample_type)
ancombc_output_cDNA.genus <-ancombc2(data= cDNA_genus.counts_filtered, 
                                    assay_name = "counts", 
                                    tax_level = "Genus",
                                    fix_formula = "sample_type + feedlot", 
                                    # rand_formula =  "(1 | feedlot)", 
                                    prv_cut = 0.05, 
                                    lib_cut = 0, 
                                    group= "sample_type", 
                                    struc_zero = TRUE, 
                                    neg_lb = TRUE,
                                    alpha = 0.05, #default significance
                                    n_cl = 1, verbose = TRUE)

## extract results from comparisons 
res.cDNA.genus <- ancombc_output_cDNA.genus$res %>%
  select(-matches("feedlot")) #not doing feedlot comparisons

#Pivot longer the results
ancom_sample_type_cDNA.genus <- res.cDNA.genus %>%
  mutate(across(starts_with("lfc_sample_type"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancom_sample_type_cDNA.genus

##Getting rid of _rounded suffix using sub command
ancom_sample_type_cDNA.genus$group ##want to get rid of "_rounded"
ancom_sample_type_cDNA.genus$group<- sub("_rounded", "", ancom_sample_type_cDNA.genus$group) 
ancom_sample_type_cDNA.genus$group #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_sample_type_cDNA.genus <- ancom_sample_type_cDNA.genus %>%
  mutate(group= case_when(
    group == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ group ##keeps original name for groups not specified (cDNA)
  ))
ancom_sample_type_cDNA.genus$group ##Now the group names are shorter and more manageable
ancom_sample_type_cDNA.genus<- ancom_sample_type_cDNA.genus %>%
  rename(Genus = taxon) ##This ancombc was done at the genus level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancom_sample_type_cDNA.genus_2 <- ancom_sample_type_cDNA.genus %>%
  mutate(lower.ci = lfc_sample_typeWater - 1.96*se_sample_typeWater,
         upper.ci = lfc_sample_typeWater + 1.96*se_sample_typeWater)

##Final fix - up to make compatible with plotting
ancom_sample_type_cDNA.genus_3 <- ancom_sample_type_cDNA.genus_2 %>%
  #Only want those that passed sensitivity testing
  filter (passed_ss_sample_typeWater == 1) %>%
  #Will only include the classified taxa
  filter(!grepl("unknown|unclassified", Genus))%>%
   rename(stderr = se_sample_typeWater, ##Renaming
         pval = p_sample_typeWater, 
         qval = q_sample_typeWater)%>%
  mutate(direction = ifelse(value >0, "elevated", "depleted"),
         # Exponentiate the values first
         value_exp = exp(value),
         lower.ci_exp = exp(lower.ci),
         upper.ci_exp = exp(upper.ci)) %>%
  mutate(
    ##Getting log 2 so I can compare with maaslin
    coef =  log2(value_exp),
    lower.ci =  log2(lower.ci_exp),
    upper.ci =  log2(upper.ci_exp),
    plot = "Log2 Fold change with 95%CI", 
    test = "ANCOM-BC", 
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  # filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Genus, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns
nrow(ancom_sample_type_cDNA.genus_3) #26 DA genera between Feces and Water with ANCOM

###MaAsLin3#######
#Data (otu counts) and metadata for MaAslin
#Counts for data_maaslin_cDNA_ra are tss transformed
data.bacteria.samples.cDNA_genus.ra.filt #Will feed it this filtered ps object (filtered for those genera with mean RA > 0.3% across 60 samples (cDNA samples))
data_maaslin_cDNA_ra  <- data.frame(t(data.bacteria.samples.cDNA_genus.ra.filt@otu_table)) 

##Sample metadata
metadata_maaslin_cDNA <- data.frame(data.bacteria.samples.cDNA_genus.ra.filt@sam_data) %>%
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")),
         sample_type = factor(sample_type, levels = c("Feces", "Water"))) #making feedlot a factor since I'll be adding it as a random effect on the MaAslin model

##Will also try out with raw counts, glommed at genus level but not filtered (so when MaAslin does tss normalization, it takes into account all OTUs)
data.bacteria.samples.cDNA.genus <- tax_glom(data.bacteria.samples.cDNA, taxrank = "Genus", NArm = F) #data.bacteria.samples.cDNA has raw counts
data_maaslin_cDNA  <- data.frame(t(data.bacteria.samples.cDNA.genus@otu_table))


##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_cDNA_genus <-  maaslin3(
  input_data = data_maaslin_cDNA_ra, 
  input_metadata = metadata_maaslin_cDNA, 
  output = "MaAsLin3_cDNA",  
  fixed_effect = c("sample_type", "feedlot"),
  # fixed_effect = "sample_type",
  # random_effects = "feedlot",
  min_prevalence=0.05,
  median_comparison_abundance = T, #default
  median_comparison_prevalence = FALSE, #default 
  min_abundance = 0, ##input_data has already been filtered (>0.3% RA)
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_cDNA_genus$fit_data_abundance$results) #Abundance results from MaAslin - 53

##Taxonomy of cDNA OTUs 
input_taxonomy_cDNA <- data.frame(data.bacteria.samples.cDNA.ra.genus@tax_table) %>%
  rownames_to_column(var = "feature")

#Calculate confidence intervals, add taxonomy
maaslin_cDNA_genus_2 <- maaslin_cDNA_genus$fit_data_abundance$results %>%
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  left_join(input_taxonomy_cDNA, by = "feature")

##Final edits to put together for plot (keeping log2 values)
maaslin_cDNA_genus_3 <- maaslin_cDNA_genus_2%>%
  #Will only include the classified taxa
  filter(!grepl("unknown|unclassified", Genus))%>%
  rename(group = name,
         pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(group= case_when(
    group == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ group ##keeps original name for groups not specified 
  ))%>%
  mutate(direction = ifelse(coef > 0, "elevated", "depleted"))%>%
  mutate(
    plot = "Log2 Fold change with 95%CI", 
    test = "MaAsLin3",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  # filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Genus, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancom_sample_type_DNA.genus_3)
nrow(maaslin_cDNA_genus_3) ##20 DA genera between Feces and Water by MaAslin

##ANCOM and MaAslin together
DA_cDNA_plot_MaAslinANCOM.data <- rbind(ancom_sample_type_cDNA.genus_3, maaslin_cDNA_genus_3) %>%
  filter(Genus %in% intersect(maaslin_cDNA_genus_3$Genus,
                              ancom_sample_type_cDNA.genus_3$Genus)) ##Only going to plot those taxa DA by both tests

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_cDNA <- maaslin_cDNA_genus$transformed_data %>% #transformed data is TSS (ra) log2 transformed otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "feature", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy_cDNA, by = "feature")%>%
  filter(Genus %in% intersect(maaslin_cDNA_genus_3$Genus,
                              ancom_sample_type_cDNA.genus_3$Genus))%>%
  left_join(metadata_maaslin_cDNA%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(Genus, Sample, logvalue, sample_type, plot)

#Bias-corrected abundances (ANCOM)
cDNA.genus_log_corr_abn <- ancombc_output_cDNA.genus$bias_correct_log_table %>%
  data.frame()%>% ##make into data frame
  rownames_to_column("Genus")%>%
  filter(Genus %in% intersect(maaslin_cDNA_genus_3$Genus,
                              ancom_sample_type_cDNA.genus_3$Genus))%>% #keep only those genera in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         sample_type = ifelse(grepl("W", Sample), "Water", "Feces"))

#put together objects to plot DA
DA_cDNA_plot_together <- bind_rows(DA_cDNA_plot_MaAslinANCOM.data, cDNA.genus_log_corr_abn, RA_MaaslinAncom_cDNA) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_cDNA_plot_together$plot <- factor(DA_cDNA_plot_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_cDNA_plot_together$sample_type <- factor(DA_cDNA_plot_together$sample_type, levels = c("Feces", "Water"))


###PLOTTING DA#
##Ordering how I want the "Genus" taxlevel to show up on the plot 
input_taxonomy_cDNA ##dataframe object for Taxonomy of cDNA OTUs 


#PHYLUM LEVEL TAXONOMY
# Create the taxonomy plot data and modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_cDNA <- DA_cDNA_plot_together %>%
  left_join(input_taxonomy_cDNA, by= "Genus")%>%
  distinct(Phylum, Order, Family, Genus) %>%
  arrange(Phylum, Order, Family) %>%
  mutate(Genus = factor(Genus, levels = rev(Genus)))%>% ##Since I arranged by family, this is the order I want the genera to show up
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Family' with only the first occurrence of each class
  ungroup()

levels(taxonomy_plot_data_cDNA$Genus)

##Factor "Genus" level by the order I want (taxonomy_plot_data_cDNA$Genus)
DA_cDNA_plot_together$Genus <- factor(DA_cDNA_plot_together$Genus, levels = rev(taxonomy_plot_data_cDNA$Genus))

# Create the updated taxonomy plot
taxonomy_plot_cDNA <- ggplot(taxonomy_plot_data_cDNA) +
  geom_text(aes(x =0, y = Genus, label = label_Phylum), hjust = 0, size = 8, family = "sans") +  # Move text left by adjusting x
  labs(title = "Phylum") +
  theme_void() +
  #scale_y_discrete(limits =  rev(taxonomy_plot_data_cDNA$Genus)) +  # Ensure y-axis matches the taxon order
  theme(plot.title = element_text(hjust = 0, size = 22, vjust = -0.5, face = "bold", family = "sans"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )+
  coord_cartesian(xlim = c(0.05, 1)) #had to add this to move the geom_text more to the left
taxonomy_plot_cDNA

#Plotting 
DA_cDNA_plot_MaAslinANCOM <-
  ggplot(data=DA_cDNA_plot_together%>%filter(grepl("abundances", plot)),
         aes(x=Genus, y=logvalue, fill = sample_type, color = sample_type)) +
  geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order = 1,title="Sample Type", title.position="top"))+
  scale_fill_manual(values=sample.type.palette, 
                    labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  scale_color_manual(values = sample.type.palette, 
                     labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  new_scale_color()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_cDNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=Genus, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_cDNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=Genus, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_cDNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed'
             , alpha=0.5) +
  facet_nested(. ~ plot, scales='free_x',
               space='free_y',
               switch='y',
               strip=strip_nested(text_y=list(element_text(angle=0))),
               labeller=labeller(group=label_wrap_gen(width=10),
                                 sub_group=label_wrap_gen(width=10))) +
  scale_x_discrete(position='bottom') +
  scale_y_continuous(position='right') +  # Default to break_labels otherwise 
  coord_flip() +
  scale_color_manual(values=c("red", "blue")) +
  scale_linetype_manual(values=c("11", "solid")) +
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Sample Type", title.position="top"),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 2.5)),
         linetype = guide_legend(title = "Fold change source", title.position = "top",
                                 override.aes = list(linewidth = 1),
                                 theme = theme(legend.key.width = unit(1.5, "cm"))),
         shape=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  labs(title = "RNA (cDNA)")+
  theme_bw()+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=26), legend.text=element_text(size=26),
        plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y=element_text(size= 22, angle=0, vjust= 1.045, face = "bold"),         
        axis.text.y=element_text(size=20, vjust = 0.5),
        strip.text=element_text(size=16, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_cDNA_plot_MaAslinANCOM

##Adding the q values
DA_cDNA_plot_MaAslinANCOM_q <- DA_cDNA_plot_MaAslinANCOM +
  geom_text(inherit.aes=FALSE,
            data = DA_cDNA_plot_together %>% filter(plot == 'Log2 Fold change with 95%CI'),
            aes(x = Genus, y = coef, label = DA),
            position = position_dodge2(width = 0.75),
            vjust = -0.6, size = 3)
DA_cDNA_plot_MaAslinANCOM_q

##Putting together DA plot with taxonomy (family) plot
combined_plot_cDNA <- plot_grid(
  taxonomy_plot_cDNA+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_cDNA_plot_MaAslinANCOM  + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                    plot.title = element_blank(),
                                    # strip.text = element_blank(),
                                    # strip.background = element_rect(fill = "white"),
                                    legend.position = "none",
                                    ),
  ncol = 2, 
  rel_widths = c(0.3,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  theme(
    legend.position = "none"
  )
  # labs(title = "Metatranscriptomic libraries (RNA (cDNA))")+
  # theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_cDNA

##Adding q values 
combined_plot_cDNA_q <- plot_grid(
  taxonomy_plot_cDNA+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_cDNA_plot_MaAslinANCOM_q + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                     plot.title = element_blank(),
                                     strip.text=element_text(size=16, color = "white", face = "bold"),
                                     legend.position = "none"),
  ncol = 2, 
  rel_widths = c(0.3,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "Metatranscriptomic libraries (RNA (cDNA))")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_cDNA_q

####FIGURE 5AandB#####
#Put together these 2 plots:
combined_plot_DNA
combined_plot_cDNA

figure5AandB <-plot_grid(combined_plot_DNA +
                           theme(plot.title = element_blank()), 
          combined_plot_cDNA +
            theme(plot.title = element_blank()), 
          align = "v",
          labels = c("A", "B"),
          label_size = 32,
          rel_heights = c(1.2,0.8),
          ncol = 1) 
  # labs(title = "MICROBIOME")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
figure5AandB
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure5AB.png", 
       plot = figure5AandB, 
       device = "png",
       width = 19, 
       height = 16, 
       dpi = 600,
       bg = "white")

#UPSET PLOTS#####
##DNA (CB vs. Feces)-Figure 4A#####
##Making the dataset for upset of UpSetR (binary matrix, present or absent)
upset.data.bacteria.samples.DNA <- MicrobiotaProcess::get_upset(data.bacteria.samples.DNA, factorNames="sample_type") 
upset.data.bacteria.samples.DNA
upset.data.bacteria.samples.DNA <- upset.data.bacteria.samples.DNA%>%
  rename("CB" = "Water")

#Plot
upset_plot_DNA <-upset(upset.data.bacteria.samples.DNA,
                       sets.bar.color = c("brown",
                                          "#4C72B0"),
                       order.by = "freq", text.scale = c(3, 2.5, 2.5, 1.5, 3, 2),
                       point.size = 5, line.size = 2, mainbar.y.label= "OTU count",
                       sets.x.label = "OTU count", 
                       set_size.show = F)
upset_plot_DNA
svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure4A.svg",
        width=7, 
        height=6)
upset_plot_DNA
dev.off()

##cDNA (CB vs. Feces) -Figure 4B#####
##making the dataset for upset of UpSetR (binary matrix, present or absent)
upset.data.bacteria.samples.cDNA <- MicrobiotaProcess::get_upset(data.bacteria.samples.cDNA, factorNames="sample_type") 
upset.data.bacteria.samples.cDNA
upset.data.bacteria.samples.cDNA <- upset.data.bacteria.samples.cDNA%>%
  rename("CB" = "Water")

#Plot
upset_plot_cDNA <-upset(upset.data.bacteria.samples.cDNA, 
                        sets.bar.color = c("brown",
                                           "#4C72B0"),
                        order.by = "freq", text.scale = c(3, 2.5, 2.5, 1.5, 3, 2),
                        point.size = 5, line.size = 2, mainbar.y.label= "OTU count",
                        sets.x.label = "OTU count", 
                        set_size.show = F)
upset_plot_cDNA

##saving the upset plot
svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure4B.svg",
        width=7, height=6)
upset_plot_cDNA
dev.off()


##Feces (cDNA vs. DNA) - Figure 8A#####
##making the dataset for upset of UpSetR (binary matrix, present or absent)
upset.data.bacteria.samples.feces <- MicrobiotaProcess::get_upset(data.bacteria.samples.feces, factorNames="gen_material") 
upset.data.bacteria.samples.feces
upset.data.bacteria.samples.feces <- upset.data.bacteria.samples.feces%>%
  rename("RNA (cDNA)" = "cDNA")

##saving the upset plot
upset_plot_feces <-UpSetR::upset(upset.data.bacteria.samples.feces, 
                                 sets.bar.color = c("#CC79A7", "#009E73"), 
                                 order.by = "freq", text.scale = c(3, 2.5, 2.5, 1.5, 3, 2),
                                 point.size = 5, line.size = 2, mainbar.y.label= "OTU count",
                                 sets.x.label = "OTU count", 
                                 set_size.show = F) 
upset_plot_feces

##saving the upset plot
svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure8A.svg",
        width=7, height=6)
upset_plot_feces
dev.off()

##CB (cDNA vs. DNA)-Figure 8B#####
##making the dataset for upset of UpSetR (binary matrix, present or absent)
upset.data.bacteria.samples.water <- MicrobiotaProcess::get_upset(data.bacteria.samples.water, factorNames="gen_material") 
upset.data.bacteria.samples.water
upset.data.bacteria.samples.water <- upset.data.bacteria.samples.water%>%
  rename("RNA (cDNA)" = "cDNA")

##saving the upset plot
upset_plot_water <-UpSetR::upset(upset.data.bacteria.samples.water, 
                                 sets.bar.color = c("#CC79A7", "#009E73"), 
                                 order.by = "freq", text.scale = c(3, 2.5, 2.5, 1.5, 3, 2),
                                 point.size = 5, line.size = 2, mainbar.y.label= "OTU count",
                                 sets.x.label = "OTU count", 
                                 set_size.show = F)
upset_plot_water
svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure8B.svg",
        width=7, height=6)
upset_plot_water
dev.off()



#SALMONELLA POSITIVES (CULTURE) AND RA####
data.bacteria.samples.tss.df ##Metadata
##Are there any samples where PCR and culture results do not match?
check_column_mismatch <- function(df, column1, column2) {
  # Filter rows where values in column1 and column2 do not match
  mismatches <- df[df[[column1]] != df[[column2]], ]
  
  # Return the "SampleID" column for the mismatched rows
  return(mismatches$SampleID)
}
check_column_mismatch(data.bacteria.samples.tss.df,
                      "salmonella_culture_status", 
                      "salmonella_PCR_results") ##F2F11 was culture negative but PCR positive

#Salmonella
salmonella_ra<- subset_taxa(data.bacteria.samples.tss, Genus =="Salmonella") ##Only need Salmonella, taking it from the RA object
salmonella_ra<- prune_samples(sample_sums(salmonella_ra) > 0, salmonella_ra) ##Getting rid of samples without any Salmonella (though they all had it)
salmonella_ra <- tax_glom(salmonella_ra, taxrank = "Genus", NArm = F) ##Tax_glomming, to get total RA of the Salmonella Genus in each sample
salmonella_ra.cDNA <- subset_samples(salmonella_ra, gen_material=="cDNA") ##Getting only cDNA samples
salmonella_ra.cDNA <- prune_taxa(taxa_sums(salmonella_ra.cDNA) > 0, salmonella_ra.cDNA)  ##Pruning taxa with 0 total counts in cDNA samples
salmonella_ra.DNA <- subset_samples(salmonella_ra, gen_material=="DNA") ##Getting only cDNA samples
salmonella_ra.DNA <- prune_taxa(taxa_sums(salmonella_ra.DNA) > 0, salmonella_ra.DNA)##Pruning taxa with 0 total counts in DNA samples
##cDNA 
#Feces
salmonella_ra.cDNA.feces<- subset_samples(salmonella_ra.cDNA, sample_type == "Feces") ##Getting only fecal cDNA samples
salmonella_ra.cDNA.feces <- prune_taxa(taxa_sums(salmonella_ra.cDNA.feces) > 0, salmonella_ra.cDNA.feces)  ##Pruning taxa with 0 total counts in fecal cDNA samples
#Catch Basins
salmonella_ra.cDNA.water<- subset_samples(salmonella_ra.cDNA, sample_type == "Water") ##Getting only water cDNA samples
salmonella_ra.cDNA.water <- prune_taxa(taxa_sums(salmonella_ra.cDNA.water) > 0, salmonella_ra.cDNA.water)  ##Pruning taxa with 0 total counts in water cDNA samples

##DNA 
#Feces
salmonella_ra.DNA.feces<- subset_samples(salmonella_ra.DNA, sample_type == "Feces") ##Getting only fecal DNA samples
salmonella_ra.DNA.feces <- prune_taxa(taxa_sums(salmonella_ra.DNA.feces) > 0, salmonella_ra.DNA.feces)  ##Pruning taxa with 0 total counts in fecal DNA samples
#Catch basins
salmonella_ra.DNA.water<- subset_samples(salmonella_ra.DNA, sample_type == "Water") ##Getting only water DNA samples
salmonella_ra.DNA.water <- prune_taxa(taxa_sums(salmonella_ra.DNA.water) > 0, salmonella_ra.DNA.water)  ##Pruning taxa with 0 total counts in water DNA samples


##Melting objects
salmonella_ra.melt <- psmelt(salmonella_ra) ##melting to long format (both cDNA and DNA samples included)
salmonella_ra.cDNA.melt <- psmelt(salmonella_ra.cDNA) ##melting to long format (only cDNA)
salmonella_ra.DNA.melt <- psmelt(salmonella_ra.DNA) ##melting to long format (only DNA)
salmonella_ra.cDNA.feces.melt<- psmelt(salmonella_ra.cDNA.feces) ##melting to long format (only feces cDNA)
salmonella_ra.DNA.feces.melt<- psmelt(salmonella_ra.DNA.feces) ##melting to long format (only feces DNA)
salmonella_ra.cDNA.water.melt<- psmelt(salmonella_ra.cDNA.water) ##melting to long format (only water cDNA)
salmonella_ra.DNA.water.melt<- psmelt(salmonella_ra.DNA.water) ##melting to long format (only water DNA)

##CULTURE STATUS######
###Feces DNA #########
salmonella_status_RA_DNA_feces <- ggplot(salmonella_ra.DNA.feces.melt, 
                                         aes(x= salmonella_culture_status, y =  Abundance, color = salmonella_culture_status, fill =salmonella_culture_status))+
  geom_point(size = 3, shape = 18) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"), 
                     labels =  c("positive"= "Positive", "negative"= "Negative"))+
  scale_fill_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"),
                    labels =  c("positive"= "Positive", "negative"= "Negative"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"), 
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 35),
        panel.border = element_rect(colour = "black", linewidth= 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "FECES (DNA)",
       y = expression(italic("Salmonella") ~ "RA (%)"),
       color = "Culture Status",
       fill = "Culture Status"
  ) +
  geom_pwc(method = "wilcox_test",label = "p = {p}",
           hide.ns = TRUE,
           label.size = 7,
           step.increase = 0.08,
           tip.length = 0.02)
salmonella_status_RA_DNA_feces
##Stat test
wilcox_test(salmonella_ra.DNA.feces.melt, Abundance ~ salmonella_culture_status)#n.s. p =  0.09

###Feces cDNA########
salmonella_status_RA_cDNA_feces <- ggplot(salmonella_ra.cDNA.feces.melt, 
                                          aes(x= salmonella_culture_status, y =  Abundance, color = salmonella_culture_status, fill =salmonella_culture_status))+
  geom_point(size = 3, shape = 18) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"), 
                     labels =  c("positive"= "Positive", "negative"= "Negative"))+
  scale_fill_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"),
                    labels =  c("positive"= "Positive", "negative"= "Negative"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"), 
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 35),
        panel.border = element_rect(colour = "black", linewidth= 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "FECES (RNA (cDNA))",
       y = expression(italic("Salmonella") ~ "RA (%)"),
       color = "Culture Status",
       fill = "Culture Status"
  ) +
  geom_pwc(method = "wilcox_test",label = "p = {p}",
           hide.ns = TRUE,
           label.size = 7,
           step.increase = 0.08,
           tip.length = 0.02)
salmonella_status_RA_cDNA_feces
##Stat test
wilcox_test(salmonella_ra.cDNA.feces.melt, Abundance ~ salmonella_culture_status)#n.s. p = 0.3

###Catch Basin DNA#######
salmonella_status_RA_DNA_water <- ggplot(salmonella_ra.DNA.water.melt, 
                                         aes(x= salmonella_culture_status, y =  Abundance, color = salmonella_culture_status, fill =salmonella_culture_status))+
  geom_point(size = 3, shape = 18) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"), 
                     labels =  c("positive"= "Positive", "negative"= "Negative"))+
  scale_fill_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"),
                    labels =  c("positive"= "Positive", "negative"= "Negative"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"), 
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 35),
        panel.border = element_rect(colour = "black", linewidth= 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "PONDS (DNA)",
       y = expression(italic("Salmonella") ~ "RA (%)"),
       color = "Culture Status",
       fill = "Culture Status"
  ) +
  geom_pwc(method = "wilcox_test",label = "p = {p}",
           hide.ns = TRUE,
           label.size = 7,
           step.increase = 0.08,
           tip.length = 0.02)
salmonella_status_RA_DNA_water
##Stat test
wilcox_test(salmonella_ra.DNA.water.melt, Abundance ~ salmonella_culture_status)#n.s. p =  0.53

###Catch basin cDNA########
salmonella_status_RA_cDNA_water <- ggplot(salmonella_ra.cDNA.water.melt, 
                                          aes(x= salmonella_culture_status, y =  Abundance, color = salmonella_culture_status, fill =salmonella_culture_status))+
  geom_point(size = 3, shape = 18, position = position_dodge(width= 0.75)) +
  geom_boxplot(alpha = 0.1, position = position_dodge2()) +
  scale_color_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"), 
                     labels =  c("positive"= "Positive", "negative"= "Negative"))+
  scale_fill_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"),
                    labels =  c("positive"= "Positive", "negative"= "Negative"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"), 
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 35),
        panel.border = element_rect(colour = "black", linewidth= 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "PONDS (RNA (cDNA))",
       y = expression(italic("Salmonella") ~ "RA (%)"),
       color = "Culture Status",
       fill = "Culture Status"
  ) +
  geom_pwc(method = "wilcox_test",label = "p = {p}",
           hide.ns = TRUE,
           label.size = 7,
           step.increase = 0.08,
           tip.length = 0.02)
salmonella_status_RA_cDNA_water
##Stat test
wilcox_test(salmonella_ra.cDNA.water.melt, Abundance ~ salmonella_culture_status)#n.s. p = 0.343


#####SUPPLEMENTARY FIGURE 5#####
salmonella_ra.melt$sample_type <- factor(salmonella_ra.melt$sample_type, levels = c("Feces", "Water"))
salmonella_ra.melt$gen_material <- factor(salmonella_ra.melt$gen_material, levels = c("DNA", "cDNA"))

sfigure5 <- ggplot(salmonella_ra.melt, 
                   aes(x= salmonella_culture_status, 
                       y =  Abundance, 
                       color = salmonella_culture_status, 
                       fill =salmonella_culture_status))+
  facet_grid(sample_type~gen_material, scales = "free",
             labeller = as_labeller(c("Feces" = "FECES", "Water" = "CATCH BASINS",
                                      "cDNA" = "RNA (cDNA)", "DNA" = "DNA")))+
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"), 
                     labels =  c("positive"= "Positive", "negative"= "Negative"))+
  scale_fill_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"),
                    labels =  c("positive"= "Positive", "negative"= "Negative"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"), 
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 35),strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(y = expression(italic("Salmonella") ~ "RA (%)"),
       color = "Culture Status",
       fill = "Culture Status"
  ) +
  geom_pwc(method = "wilcox_test",label = "p = {p}",
           hide.ns = TRUE,
           label.size = 7,
           step.increase = 0.08,
           tip.length = 0.02)
sfigure5 
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryFigure5.png", 
       plot = sfigure5, 
       device = "png", 
       dpi = 600, 
       width = 11, height =10)

##TABLE 3 - Average RA of salmonella in sample types #######
salmonella_ra.melt %>%
  group_by(gen_material, sample_type)%>%
  summarise(mean_salmonella_abundance = mean(Abundance),
            sd_salmonella_abundance = sd(Abundance))

#How many samples per group? 
salmonella_ra.melt %>%
  group_by(gen_material, sample_type, salmonella_culture_status)%>%
  count()
# gen_material sample_type salmonella_culture_status     n
# DNA          Feces       negative                     29
# DNA          Feces       positive                     19
# DNA          Water       negative                      7
# DNA          Water       positive                      5
# cDNA         Feces       negative                     29
# cDNA         Feces       positive                     19
# cDNA         Water       negative                      7
# cDNA         Water       positive                      5

##PCR STATUS#######
###Feces DNA #######
salmonella_statusPCR_RA_DNA_feces <- ggplot(salmonella_ra.DNA.feces.melt, 
                                            aes(x= salmonella_PCR_results, y =  Abundance, color = salmonella_PCR_results, fill =salmonella_PCR_results))+
  geom_point(size = 3, shape = 18) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"), 
                     labels =  c("positive"= "Positive", "negative"= "Negative"))+
  scale_fill_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"),
                    labels =  c("positive"= "Positive", "negative"= "Negative"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"), 
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 35),
        panel.border = element_rect(colour = "black", linewidth= 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "FECES (DNA)",
       y = expression(italic("Salmonella") ~ "RA (%)"),
       color = "PCR Status",
       fill = "PCR Status"
  ) +
  geom_pwc(method = "wilcox_test",label = "p = {p}",
           hide.ns = TRUE,
           label.size = 7,
           step.increase = 0.08,
           tip.length = 0.02)
salmonella_statusPCR_RA_DNA_feces
##Stat test
wilcox_test(salmonella_ra.DNA.feces.melt, Abundance ~ salmonella_PCR_results)#n.s. p =  0.19

###Feces cDNA samples plot#######
salmonella_statusPCR_RA_cDNA_feces <- ggplot(salmonella_ra.cDNA.feces.melt, 
                                             aes(x= salmonella_PCR_results, y =  Abundance, color = salmonella_PCR_results, fill =salmonella_PCR_results))+
  geom_point(size = 3, shape = 18) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"), 
                     labels =  c("positive"= "Positive", "negative"= "Negative"))+
  scale_fill_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"),
                    labels =  c("positive"= "Positive", "negative"= "Negative"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"), 
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 35),
        panel.border = element_rect(colour = "black", linewidth= 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "FECES RNA (cDNA)",
       y = expression(italic("Salmonella") ~ "RA (%)"),
       color = "PCR Status",
       fill = "PCR Status"
  ) +
  geom_pwc(method = "wilcox_test",label = "p = {p}",
           hide.ns = TRUE,
           label.size = 7,
           step.increase = 0.08,
           tip.length = 0.02)
salmonella_statusPCR_RA_cDNA_feces
##Stat test
wilcox_test(salmonella_ra.cDNA.feces.melt, Abundance ~ salmonella_PCR_results)#n.s. p = 0.461


###Catch basin DNA #######
salmonella_statusPCR_RA_DNA_water <- ggplot(salmonella_ra.DNA.water.melt, 
                                            aes(x= salmonella_PCR_results, y =  Abundance, color = salmonella_PCR_results, fill =salmonella_PCR_results))+
  geom_point(size = 3, shape = 18) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"), 
                     labels =  c("positive"= "Positive", "negative"= "Negative"))+
  scale_fill_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"),
                    labels =  c("positive"= "Positive", "negative"= "Negative"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"), 
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 35),
        panel.border = element_rect(colour = "black", linewidth= 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "PONDS (DNA)",
       y = expression(italic("Salmonella") ~ "RA (%)"),
       color = "PCR Status",
       fill = "PCR Status"
  ) +
  geom_pwc(method = "wilcox_test",label = "p = {p}",
           hide.ns = TRUE,
           label.size = 7,
           step.increase = 0.08,
           tip.length = 0.02)
salmonella_statusPCR_RA_DNA_water
##Stat test
wilcox_test(salmonella_ra.DNA.water.melt, Abundance ~ salmonella_PCR_results)#n.s. p =  0.53


###Catch basin cDNA #######
salmonella_statusPCR_RA_cDNA_water <- ggplot(salmonella_ra.cDNA.water.melt, 
                                             aes(x= salmonella_PCR_results, y =  Abundance, color = salmonella_PCR_results, fill =salmonella_PCR_results))+
  geom_point(size = 3, shape = 18) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"), 
                     labels =  c("positive"= "Positive", "negative"= "Negative"))+
  scale_fill_manual(values = c("positive"= "#fc8d62", "negative"= "#8da0cb"),
                    labels =  c("positive"= "Positive", "negative"= "Negative"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"), 
        legend.text = element_text(size = 24),
        plot.title = element_text(size = 35),
        panel.border = element_rect(colour = "black", linewidth= 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 28),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "PONDS (RNA (cDNA))",
       y = expression(italic("Salmonella") ~ "RA (%)"),
       color = "PCR Status",
       fill = "PCR Status"
  ) +
  geom_pwc(method = "wilcox_test",label = "p = {p}",
           hide.ns = TRUE,
           label.size = 7,
           step.increase = 0.08,
           tip.length = 0.02)
salmonella_statusPCR_RA_cDNA_water
##Stat test
wilcox_test(salmonella_ra.cDNA.water.melt, Abundance ~ salmonella_PCR_results)#n.s. p = 0.2


##SALMONELLA SPECIES DENDROGRAM#####
data.bacteria.samples.Salmonella<- subset_taxa(data.bacteria.samples, Genus == "Salmonella") ##Choosing only those OTUs that are part of the Salmonella Genus
data.bacteria.samples.Salmonella <- prune_samples(sample_sums(data.bacteria.samples.Salmonella) > 0, data.bacteria.samples.Salmonella) ## Keep only samples with nonzero Salmonella reads
data.bacteria.samples.Salmonella ##24 OTUs are Salmonella (120 samples)
setdiff(sample_names(data.bacteria.samples), sample_names(data.bacteria.samples.Salmonella))##All samples had Salmonella 

#Relative abundance 
data.bacteria.samples.Salmonella.tss <- transform_sample_counts(data.bacteria.samples.Salmonella, function(x) x/sum(x)*100) ##normalizing counts

##Creating distance matrix
data.bacteria.samples.Salmonella.tss.bray <- vegdist(t(data.bacteria.samples.Salmonella.tss@otu_table), method = "bray") 

##Clustering (Wards distance) at the species level 
data.bacteria.samples.Salmonella.bray.hclust <- hclust(data.bacteria.samples.Salmonella.tss.bray, method = "ward.D2")
plot(data.bacteria.samples.Salmonella.bray.hclust, hang = -1)

##Building dendrogram
data.bacteria.samples.Salmonella.bray.dendro <- as.dendrogram(data.bacteria.samples.Salmonella.bray.hclust) # Build dendrogram object from hclust results
data.bacteria.samples.Salmonella.bray.dendro.data <- dendro_data(data.bacteria.samples.Salmonella.bray.dendro, type = "rectangle") # Extract the dendrogram plot data
data.bacteria.samples.Salmonella.tss@sam_data$sampleID <- rownames(data.bacteria.samples.Salmonella.tss@sam_data) ##need a column with sample_ID
data.bacteria.samples.Salmonella.bray.dendro.metadata <- 
  data.frame(data.bacteria.samples.Salmonella.tss@sam_data) %>% ##making a metadata dataframe
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"),
         gen_material.abbrv = dplyr::recode(gen_material, "cDNA"= "T", "DNA"= "G"))##adding a shorter name for Water and Feces
  
##Ordering sample labels
data.bacteria.samples.Salmonella.bray.dendro.data$labels <- data.bacteria.samples.Salmonella.bray.dendro.data$labels %>%
  left_join(data.bacteria.samples.Salmonella.bray.dendro.metadata, by = c("label" = "sampleID")) 

##Plotting dendrogram
dendro.bray.plot.Salmonella <- ggplot(data.bacteria.samples.Salmonella.bray.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = data.bacteria.samples.Salmonella.bray.dendro.data$labels, aes(x=x, y=y, colour = sample_type),
             size = 4, shape = 15, position = position_nudge(y = -0.03)) +
  scale_color_manual(name = "Sample Type", values = sample.type.palette,
                     labels = c("Water" = "CB", "Feces" = "Feces")) +
  guides(color = guide_legend (title.position = "top", override.aes = list(size = 8)))+
  new_scale_color()+ 
  geom_point(data = data.bacteria.samples.Salmonella.bray.dendro.data$labels,  aes(x = x, y = y, color = factor(feedlot)), 
             size = 4, shape = 15, position = position_nudge(y = -0.08)) +
  scale_color_manual(name = "Feedlot", values = feedlot_palette) +
  guides(color = guide_legend (title.position = "top", override.aes = list(size = 8)))+
  new_scale_color()+ 
  geom_point(data = data.bacteria.samples.Salmonella.bray.dendro.data$labels,  aes(x = x, y = y, color = gen_material), 
             size = 4, shape = 15, position = position_nudge(y = -0.13)) +
  scale_color_manual(name = "Library Type", values = gen.material.palette,
                     labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)")) +
  geom_text(data = data.bacteria.samples.Salmonella.bray.dendro.data$labels, aes(x=x, y=y, label = sample_type.abbrv), 
            colour = "white", size =2, fontface = "bold", position = position_nudge(y = -0.03)) +
  geom_text(data = data.bacteria.samples.Salmonella.bray.dendro.data$labels, 
            aes(x=x, y=y, label = factor(feedlot)), colour = "white", 
            size =3, position = position_nudge(y=-0.08), fontface = "bold")+
  geom_text(data = data.bacteria.samples.Salmonella.bray.dendro.data$labels, aes(x=x, y=y, label = gen_material.abbrv),
            colour = "white", size =3, position = position_nudge(y=-0.13), fontface = "bold") +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 18),
        legend.title = element_text(face = "bold", size = 18),
        plot.title = element_text(size = 30),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  guides(color = guide_legend (title.position = "top", 
                               override.aes = list(size = 5)))
dendro.bray.plot.Salmonella 


### RA plot
dendro_bray_data.bacteria.samples.Salmonella_order <- data.bacteria.samples.Salmonella.bray.dendro.data$labels$label ##order the samples
data.bacteria.samples.Salmonella.species.ra<- tax_glom(data.bacteria.samples.Salmonella.tss, taxrank = "Species") #Glom at the species level

#Merge low abundance salmonella species
data.bacteria.samples.Salmonella.species.ra.merged <- merge_low_abundance_ra(data.bacteria.samples.Salmonella.species.ra, 
                                                                     threshold = 0.5)

#Melt to long format
data.bacteria.samples.Salmonella.species.ra.merged.melt <- psmelt(data.bacteria.samples.Salmonella.species.ra.merged)
unique(data.bacteria.samples.Salmonella.species.ra.merged.melt$Species) #Not sure why I am getting NA as a name, going to manually substitute it
data.bacteria.samples.Salmonella.species.ra.merged.melt[is.na(data.bacteria.samples.Salmonella.species.ra.merged.melt$Species), "Species"] <- "Others < 0.5% RA"
#Now, factor
data.bacteria.samples.Salmonella.species.ra.merged.melt <- data.bacteria.samples.Salmonella.species.ra.merged.melt%>%
  mutate(Species = factor(Species, 
                        levels = c(setdiff(Species, 
                                           unique(grep("Others", Species, value = TRUE))), 
                                   unique(grep("Others", Species, value = TRUE)))))##Factoring the Species column so that "Others.." is the last category

####Relative abundances of each of the top Salmonella species#####
#Salmonella enterica
data.bacteria.samples.Salmonella.species.ra.merged.melt%>%
  filter(Species == "Salmonella enterica")%>%
  summarise(mean_abundance = mean(Abundance),
            sd_abundance = sd(Abundance)) #76.7%
#Salmonella bongori
data.bacteria.samples.Salmonella.species.ra.merged.melt%>%
  filter(Species == "Salmonella bongori")%>%
  summarise(mean_abundance = mean(Abundance),
            sd_abundance = sd(Abundance)) #6.2%

#Unclassified Salmonella
data.bacteria.samples.Salmonella.species.ra.merged.melt%>%
  filter(Species == "unclassified Salmonella")%>%
  summarise(mean_abundance = mean(Abundance),
            sd_abundance = sd(Abundance)) #15.97%

#Salmonella sp.
data.bacteria.samples.Salmonella.species.ra.merged.melt%>%
  filter(Species == "Salmonella sp.")%>%
  summarise(mean_abundance = mean(Abundance)) #0.57%


##Color palette for Salmonella
Salmonella_palette <- distinctColorPalette(10)
Salmonella_palette_names <- unique(data.bacteria.samples.Salmonella.species.ra.merged.melt$Species)# Create a named vector for the palette, where the names correspond to Species names
Salmonella_palette_names <- setNames((Salmonella_palette )[1:length(Salmonella_palette_names)], 
                                     Salmonella_palette_names) ##Assign a color to the species names
Salmonella_palette_names$'unclassified Salmonella' <- "grey30"
Salmonella_palette_names$'Salmonella bongori' <- "#8fd7d7"
Salmonella_palette_names$'Salmonella sp.' <- "purple"
Salmonella_palette_names$'Salmonella enterica' <- "#e25759"
Salmonella_palette_names$ "Others < 0.5% RA" <- "grey95"


##PLOT
dendroRA.Salmonella.species.plot <- ggplot(data.bacteria.samples.Salmonella.species.ra.merged.melt, aes(x=Sample, y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_bray_data.bacteria.samples.Salmonella_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =Salmonella_palette_names) +
  theme(legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 18),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())+
  guides(fill = guide_legend(title.position = "top", nrow = 2))
dendroRA.Salmonella.species.plot

#RA with Dendrogram
dendroRA.Salmonella.species.plot.2 <- plot_grid(dendro.bray.plot.Salmonella , 
                                                dendroRA.Salmonella.species.plot, 
                                                align = "v", 
                                                ncol = 1, 
                                                rel_heights = c(0.4,0.6))
dendroRA.Salmonella.species.plot.2

####FIGURE 2#####
figure2 <- plot_grid(dendroRA.Salmonella.species.plot.2, 
                     WaterandFeces.DNA_salmonella_only_BC_beta_div, 
                     align = "v", 
                     nrow = 1,
                     labels = c("A", ""),
                     label_size = 22,
                     rel_widths  = c(0.8, 0.2))
figure2
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure2.png", 
       plot = figure2, 
       dpi = 600,
       device = "png", 
       width = 20, 
       height =10,
       bg = "white")
