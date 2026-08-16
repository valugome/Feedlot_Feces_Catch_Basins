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
  "Polychrome", "colorspace", "devtools", "remotes", "patchwork"
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


####Importing the count tables
##KO groups raw table, to get reads mapped 
unmapped_reads <- read.csv("Data/HUMAnN_GeneFamilies_Unmapped_Feedlot_CatchBasins.tsv",
                         sep = "\t", 
                         check.names = F,
                         row.names = "# Gene Family HUMAnN v4.0.0.alpha.1 Adjusted CPMs")

#KO group RPK abundances
KO_groups <- read.csv('Data/HUMAnN_KO_GroupCounts_Feedlot_CatchBasins.tsv',
                      sep = "\t", 
                      check.names = F,
                      row.names = "# Gene Family HUMAnN v4.0.0.alpha.1 Adjusted CPMs")

#Pathway RPK abundances
pathways <- read.csv('Data/HUMAnN_Pathways_Counts_Feedlot_CatchBasins.tsv',
                     sep = "\t", 
                     check.names = F,
                     row.names = "# Pathway HUMAnN v4.0.0.alpha.1")

##Cleaning up names 
##KO groups
##Changing col names (sample IDs), to match them as they are in the metadata file
KO_groups_newcolnames <- sapply(strsplit(colnames(KO_groups), "_"), `[`, 1) #Splitting col names by "_", then extracting the first part of each split column name 
colnames(KO_groups) <- KO_groups_newcolnames #replacing col names for new ones
colnames(KO_groups) ##good, now sample names are OK!

##Pathway abundances
##Changing col names (sample IDs), to match them as they are in the metadata file
pathways_newcolnames <- sapply(strsplit(colnames(pathways), "_"), `[`, 1)#Splitting col names by "_", then extracting the first part of each split column name 
colnames(pathways) <- pathways_newcolnames  #replacing col names for new ones
colnames(pathways) ##good, now sample names are OK

#Change zymo- to zymo. as they are in metadata
KO_groups <- KO_groups %>%
  rename("Zymo.1a" = "Zymo-1a",
         "Zymo.1b" = "Zymo-1b")
pathways <- pathways %>%
  rename("Zymo.1a" = "Zymo-1a",
         "Zymo.1b" = "Zymo-1b")

#METADATA #####
metadata <- read.csv('Data/Metadata_Feedlot_CatchBasins.csv', 
                     check.names = F,
                     row.names = "sampleID")

metadata$SampleID<- rownames(metadata) #Making a SampleID column 
metadata$original_sample <- sub("c$", "", metadata$SampleID) #Making a column for the sample that both the metagenome and metatranscriptomes come from (match 2)


##Host free reads####
hostrem <- read.csv('Data/HostRem_Reads_Feedlot_CatchBasins.csv')

#Change zymo- to zymo. as they are in metadata
hostrem[which(hostrem$SampleID == "Zymo-1a"), "SampleID"] <- "Zymo.1a"
hostrem[which(hostrem$SampleID == "Zymo-1b"), "SampleID"] <- "Zymo.1b"

#Merge with metadata
hostrem <- hostrem %>%
  select(SampleID, Hostrem_output_total_num_seqs)%>%
  left_join(metadata, by = "SampleID")%>%
  rename(HostFree_Reads=Hostrem_output_total_num_seqs)#Hostrem_output_total_num_seqs is the total number of host free reads
rownames(hostrem) <- hostrem$SampleID

##Make into phyloseq object
##OTU tables for KO and pathways
OTU_table_KO <-phyloseq::otu_table(KO_groups%>%as.matrix(), taxa_are_rows = TRUE) ##otu table from KO_groups
OTU_table_pathways <- phyloseq::otu_table(pathways%>%as.matrix(), taxa_are_rows = T) ##otu table from 

#Sample data
sampledata <-sample_data(metadata) ##use phyloseq function sample_data() to make metadata into phyloseq sample data object

#COLOR PALETTES#####
feedlot_palette <- c("1" = "#fcca46", 
                     "2" = "#fe7f2d", 
                     "3" = "#233d4d", 
                     "4"= "#3b9ab2", 
                     "5"= "#e1b6ff")
###for sample type
sample.type.palette <- c("Water" = "#4C72B0",
                         "Feces" = "brown") 
# Library Type
gen.material.palette <- c("cDNA" = "#009E73",  
                          "DNA"  = "#CC79A7" )  

#sequencing batches
batch_palette <- c("no" = "#d19bac", 
                   "yes" = "#6a9c55") #sequencing batches

#PHYLOSEQ OBJECT####
##KO####
KO_ps <- phyloseq(OTU_table_KO,sampledata)
KO_ps ##8209 taxa (KO gene families) and 126 samples (120 samples plus 3 no template controls, 1 EB, 2 zymo mock communities )

##PATHWAYS####
pathways_ps <- phyloseq(OTU_table_pathways,sampledata)
pathways_ps ##765 taxa (pathways) and 126 samples (120 samples plus 3 no template controls, 1 EB, 2 zymo mock communities )


###SUBSETTING ONLY SAMPLES####
####KO######
KO_ps <- subset_samples(KO_ps, sample_type=="Water" | sample_type=="Feces")
KO_ps <- prune_taxa(taxa_sums(KO_ps) > 0, KO_ps) 
KO_ps ##7205 taxa (KO gene families) in 120 samples 

####PATHWAYS######
pathways_ps <- subset_samples(pathways_ps, sample_type=="Water" | sample_type=="Feces")
pathways_ps <- prune_taxa(taxa_sums(pathways_ps) > 0, pathways_ps) 
pathways_ps ##741 taxa (metabolic pathways) in 120 samples


#PERCENTAGE OF READS ALIGNED BY HUMAnN (NUCLEOTIDE AND TRANSLATED)####
#Total unmapped reads
total_unmapped_reads <- unmapped_reads
#Make into df 
total_unmapped_reads_df <- data.frame(
  SampleID = names(total_unmapped_reads),
  Unmapped_reads = as.numeric(total_unmapped_reads))

#Merge with hostrem data 
reads_hostrem_humann_mapping <- left_join(
  hostrem,
  total_unmapped_reads_df,
  by = "SampleID"
)%>%
  mutate(mapped_reads = (HostFree_Reads- Unmapped_reads))%>% #Now, mapped reads = hostfreereads - unmapped reads
  mutate(Percentage_reads_mapped = (mapped_reads / HostFree_Reads) * 100) #And for percentage calculation of mapped reads

#Stats
summary(reads_hostrem_humann_mapping$Percentage_reads_mapped)
sort(reads_hostrem_humann_mapping$Percentage_reads_mapped) #Ok, not dropping any

#Descriptive stats per group
reads_hostrem_humann_mapping %>%
  filter(sample_type %in% c("Water", "Feces"))%>%
  group_by(sample_type, gen_material)%>%
  summarise(mean_percentage_mapped_reads = mean(Percentage_reads_mapped),
            sd_percentage_mapped_reads = sd(Percentage_reads_mapped),
            min_percentage_mapped_reads = min(Percentage_reads_mapped),
            max_percentage_mapped_reads = max(Percentage_reads_mapped))
# sample_type gen_material mean_percentage_mapped_reads sd_percentage_mapped_reads min_percentage_mapped_reads max_percentage_mapped_reads
# Feces       DNA                                  65.7                       2.40                        60.4                        70.8
# Feces       cDNA                                 66.7                       1.96                        63.5                        70.7
# Water       DNA                                  62.9                       2.11                        58.2                        65.1
# Water       cDNA                                 62.4                       3.66                        55.8                        66.9

##Keep only samples for plots 
reads_hostrem_humann_mapping_samples <- reads_hostrem_humann_mapping %>%
  filter(sample_type %in% c("Water", "Feces"))

#Factor variables 
reads_hostrem_humann_mapping_samples$sample_type <- factor(
  reads_hostrem_humann_mapping_samples$sample_type,
  levels = c("Feces", "Water")
)
reads_hostrem_humann_mapping_samples$gen_material <- factor(
  reads_hostrem_humann_mapping_samples$gen_material,
  levels = c("DNA", "cDNA")
)

##PLOTS####
###cDNA vs DNA faceted by CB and Feces ####
percent_alignedreads.cDNAvsDNA.WaF<- ggplot(reads_hostrem_humann_mapping_samples%>%arrange(SampleID),
                                                aes(x = gen_material, y= Percentage_reads_mapped, color = gen_material, fill = gen_material)) +
  theme_bw() +
  labs(y= "Percentage (%) of\nAligned Reads", color = "Library Type", fill = "Library Type") +
  facet_grid(~sample_type, scales = "free",  labeller = as_labeller(c("Feces" = "FECES", "Water" = "CATCH BASINS"))) +
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
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))  +
  geom_pwc(method = "wilcox_test",
           label = "Wilcoxon, p = {p}",
           hide.ns = TRUE,
           method.args = list(paired = T),
           step.increase = 0.08,
           label.size = 6,
           tip.length = 0.02)
percent_alignedreads.cDNAvsDNA.WaF

#Stats
water.reads_hostrem_humann_mapping <- subset(reads_hostrem_humann_mapping, sample_type == "Water")
wilcox_test(water.reads_hostrem_humann_mapping%>%arrange(SampleID), 
            Percentage_reads_mapped ~ gen_material,
            paired = T) #n.s., p = 0.91
feces.reads_hostrem_humann_mapping <- subset(reads_hostrem_humann_mapping, sample_type == "Feces")
wilcox_test(feces.reads_hostrem_humann_mapping%>%arrange(SampleID), 
            Percentage_reads_mapped ~ gen_material,
            paired = T) #s., p = 0.00254

###CB vs Feces faceted by cDNA and DNA#####
percent_alignedreads_WvF.cDNAandDNA<- ggplot(reads_hostrem_humann_mapping_samples, 
                                             aes(x = sample_type, 
                                                                               y= Percentage_reads_mapped, color = sample_type, fill = sample_type)) +
  theme_bw() +
  labs(y= "Percentage (%) of\nAligned Reads", color = "Sample Type", fill = "Sample Type") +
  facet_grid(~gen_material,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("DNA" = "DNA", 
                                      "cDNA" = "RNA (cDNA)"))) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  scale_color_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  scale_fill_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 24, colour = "black"),
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
percent_alignedreads_WvF.cDNAandDNA

#Stats
cDNA.reads_hostrem_humann_mapping <- subset(reads_hostrem_humann_mapping, gen_material == "cDNA")
wilcox_test(cDNA.reads_hostrem_humann_mapping, 
            Percentage_reads_mapped ~ sample_type) #s., p = 0.0000675
DNA.reads_hostrem_humann_mapping <- subset(reads_hostrem_humann_mapping, gen_material == "DNA")
wilcox_test(DNA.reads_hostrem_humann_mapping, 
            Percentage_reads_mapped~sample_type) #s., p = 0.000385

###Supplementary Figure? Percentage aligned reads####
percent_alignedreads_WvF.cDNAandDNA
percent_alignedreads.cDNAvsDNA.WaF
sfigureX_CD <- plot_grid(percent_alignedreads_WvF.cDNAandDNA + 
                        theme(axis.title.y = element_text(size = 19)),
                        percent_alignedreads.cDNAvsDNA.WaF+ 
                        theme(axis.title.y = element_text(size = 19)), 
                      align = "v",
                      ncol = 1,
                      labels = c("C", "D"),
                      label_size = 22)
  #labs(title = "NON-HOST READ ALIGNMENT TO FUNCTIONAL GENES\n (HUMAnN)")+
  #theme(plot.title = element_text(size = 30, face = "bold"))
sfigureX_CD


#UMAPPED AND UNGROUPED#########
##KO GROUPS############
#There were no unmapped groups
KO_ungrouped <- prune_taxa(taxa_names(KO_ps) == "UNGROUPED", KO_ps)
KO_ungrouped
summary(sample_sums(KO_ungrouped))

##Percentages of ungrouped in total abundance 
ungrouped_KO_abundance <- KO_ps %>%
  transform_sample_counts(., function(x) x/sum(x)*100) %>%#convert to relative abundance
  prune_taxa(taxa_names(.) == "UNGROUPED", .)%>% #get only ungrouped 
  psmelt()%>%
  group_by(OTU, sample_type, gen_material) %>%  #group by OTU
  summarize(`mean group ra abundance` = mean(Abundance),
            `min group ra abundance` = min(Abundance),
            `max group ra abundance` = max(Abundance))%>% # Mean abundance per OTU
  rename("Sample type" = sample_type,
         "Library Type" = gen_material,
         "KO group" = OTU,
         )%>%
  ungroup()
ungrouped_KO_abundance

#Mean abundance across all samples of ungrouped gene families:
ungrouped_KO_abundance %>%
  summarise(mean_group_tss_abundance = mean(`mean group ra abundance`),
            sd_group_tss_abundance = sd(`mean group ra abundance`))

#How many actually 'grouped'?
KO_ps %>%
  transform_sample_counts(., function(x) x/sum(x)*100) %>%#convert to relative abundance
  prune_taxa(!taxa_names(.) == "UNGROUPED", .)%>%
  psmelt()%>%
  group_by(sample_type, gen_material) %>%  #group by OTU
  summarize(mean_grouped_KO_tss_abundance = mean(Abundance),
            sd_grouped_KO_tss_abundance = sd(Abundance))%>%# Mean abundance per OTU
  rename("Sample type" = sample_type,
         "Library Type" = gen_material)

##PATHWAYS#########
pathways_unmapped <- prune_taxa(taxa_names(pathways_ps) == "UNMAPPED", pathways_ps)
pathways_unmapped
summary(sample_sums(pathways_unmapped))

pathways_unintegrated <- prune_taxa(taxa_names(pathways_ps) == "UNINTEGRATED", pathways_ps)
pathways_unintegrated
summary(sample_sums(pathways_unintegrated))

##Percentages of unmapped and unintegrated in total abundance 
unintegrated_pathways_abundance <- pathways_ps %>%
  transform_sample_counts(., function(x) x/sum(x)*100) %>%#convert to relative abundance
  prune_taxa(taxa_names(.) %in% c("UNINTEGRATED"), .)%>%
  psmelt()%>%
  group_by(OTU, sample_type, gen_material) %>%  #group by OTU
  summarize(`Mean Pathway Abundance (%)` = mean(Abundance),
            `SD Pathway Abundance (%)` = sd(Abundance),
            `Min Pathway Abundance (%)` = min(Abundance),
            `Max Pathway Abundance (%)` = max(Abundance))%>% # Mean abundance per OTU
  rename("Sample type" = sample_type,
         "Library Type" = gen_material,
         "Pathway" = OTU,
  )%>%
  ungroup()
unintegrated_pathways_abundance

#Mean abundance across all samples of unintegrated and unmapped pathways:
unintegrated_pathways_abundance %>%
  group_by(Pathway)%>%
  summarise(mean_group_tss_abundance = mean(`Mean Pathway Abundance (%)`),
            sd_group_tss_abundance = sd(`Mean Pathway Abundance (%)`))

#Dropping Unmapped and Ungrouped####
KO_ps #7205 taxa 
KO_ps_filt <- prune_taxa(!(taxa_names(KO_ps) %in% c("UNMAPPED", "UNGROUPED")), KO_ps)
KO.tss_ps <- transform_sample_counts(KO_ps_filt, function(x) x/sum(x))
KO.tss_ps #7204 KO gene families (excluding "UNMAPPED" and "UNGROUPED")

pathways_ps #741 taxa 
pathways_ps_filt <- prune_taxa(!(taxa_names(pathways_ps) %in% c("UNMAPPED", "UNINTEGRATED")), pathways_ps)
pathways.tss_ps <- transform_sample_counts(pathways_ps_filt, function(x) x/sum(x))
pathways.tss_ps ##739 pathways (excluding "UNMAPPED" and "UNINTEGRATED")


#BETA_DIV ####
##BRAY CURTIS####
###GENE FAMILIES#####
#### FECES SAMPLES (DNA vs cDNA)#######
##Subsetting only feces samples
KO.feces.tss <- subset_samples(KO.tss_ps, sample_type=="Feces")
KO.feces.tss <- prune_taxa(taxa_sums(KO.feces.tss) > 0, KO.feces.tss) 
KO.feces.tss #96 samples
##Distance matrix
KO.feces.tss.bray <- vegdist(t(KO.feces.tss@otu_table), method = "bray") 
KO.feces.tss.df <- as(KO.feces.tss@sam_data,"data.frame") # make DF from metadata

##Adding a column concatenating "gen_material" and "feedlot", then making feedlot a factor
KO.feces.tss.df<- KO.feces.tss.df %>%
  mutate (gen_material_spec_2 = paste(gen_material, feedlot, sep = '_'),
          feedlot = factor(feedlot))

#####PERMANOVA#####
#Is there an interaction? 
set.seed(87)
feces_KO_BC_adonis_interaction <- adonis2(KO.feces.tss.bray ~ gen_material * feedlot,
                                          by = "margin",
                                          KO.feces.tss.df, 
                                          p.adjust.methods = "BH", permutations = 9999)
feces_KO_BC_adonis_interaction #No

##Model Library Type and feedlot, stratified by original sample from which DNA and RNA were extracted 
set.seed(87)
feces_KO_BC_adonis <- adonis2(KO.feces.tss.bray ~ gen_material + feedlot, 
                              strata = KO.feces.tss.df$original_sample,
                              by = "margin",
                              KO.feces.tss.df, 
                              p.adjust.methods = "BH", permutations = 9999)
feces_KO_BC_adonis #18.6% of the variation is due to Library Type, p = 0.0001
#12.1% of the variation is due to feedlot, p = 0.0001


#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.feces.feedlot.disp.KO <- betadisper(KO.feces.tss.bray, KO.feces.tss.df$feedlot)
bray.feces.feedlot.disp.KO
##Then test by permuting
set.seed(87)
bray.feces.feedlot.permdisp.KO <- permutest(bray.feces.feedlot.disp.KO, permutations = 9999, pairwise = 1)
bray.feces.feedlot.permdisp.KO 
##Different dispersions of variance ( p = 0.0114)

#PERMDISP- Gen_material
# Run the betadisper function, average distance to centroid
bray.feces.genmat.disp.KO <- betadisper(KO.feces.tss.bray, KO.feces.tss.df$gen_material)
bray.feces.genmat.disp.KO
##Then test by permuting
set.seed(87)
bray.feces.genmat.permdisp.KO <- permutest(bray.feces.genmat.disp.KO, permutations = 9999)
bray.feces.genmat.permdisp.KO 
#No difference in dispersions of variance between DNA and cDNA (p = 0.108)

#### ORDINATION
set.seed(87)
KO.feces.tss.bray.ord <- metaMDS(KO.feces.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

#Extract ordination coordinates, add metadata
KO.feces.tss.bray.plot <- ordiplot(KO.feces.tss.bray.ord$points)
#Extract the scores (coordinates) of the points from the ordination plot object 
KO.feces.tss.bray.scrs <- scores(KO.feces.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
KO.feces.tss.bray.scrs <- cbind(as.data.frame(KO.feces.tss.bray.scrs), 
                                   gen_material = KO.feces.tss.df$gen_material, 
                                   feedlot = factor(KO.feces.tss.df$feedlot), 
                                   gen_material_spec_2 = KO.feces.tss.df$gen_material_spec_2,
                                   sampleID = rownames(KO.feces.tss.df))

#####FEEDLOT EFFECT#####
## BC
KO.feces.tss.bray.scrs #have coodrinates and metadata
##Centroids according to feedlot
KO.feces.feedlot.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = KO.feces.tss.bray.scrs, FUN = mean) 
KO.feces.feedlot.tss.bray.segs <- merge(KO.feces.tss.bray.scrs, 
                                           setNames(KO.feces.feedlot.tss.bray.cent, c("feedlot","cMDS1","cMDS2")),
                                           by = 'feedlot', sort = F)

# Extract R2 and p-values
##Model gen_material and feedlot
feces_KO_BC_adonis #have the model
R2_feces_feedlot_BC_adonis_KO <- feces_KO_BC_adonis$R2[2] 
pvalue_feces_feedlot_BC_adonis_KO<-feces_KO_BC_adonis$`Pr(>F)`[2]

#### PLOT
fecesKO_feedlot_BC_beta_div <- ggplot(KO.feces.feedlot.tss.bray.segs) + theme_bw() +
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
  guides(color = guide_legend(override.aes = list(size = 7)))+
  annotate("text", x = -0.4, y = 0.1, ##change coordinates as needed
           label = "Feedlot", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (feedlot)
  annotate("text", x = -0.4, y = 0.1, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces_feedlot_BC_adonis_KO * 100, 1), "%",
                         "\np = ", round(pvalue_feces_feedlot_BC_adonis_KO, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
fecesKO_feedlot_BC_beta_div 

#####LIBRARY TYPE EFFECT#####
## BC
KO.feces.tss.bray.scrs #have coordinates and metadata
##Centroids according to gen_material (cDNA/DNA)
KO.feces.genmat.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ gen_material, 
                                              data = KO.feces.tss.bray.scrs, FUN = mean) 
KO.feces.genmat.tss.bray.segs <- merge(KO.feces.tss.bray.scrs, 
                                          setNames(KO.feces.genmat.tss.bray.cent, c("gen_material","cMDS1","cMDS2")), by = 'gen_material', sort = F)

# Extract R2 and p-values for genmat
feces_KO_BC_adonis #have the model
R2_feces_KO_BC_adonis_genmat <- feces_KO_BC_adonis$R2[1] 
pvalue_feces_KO_BC_adonis_genmat<-feces_KO_BC_adonis$`Pr(>F)`[1]

#### PLOT
feces_genmat_KO_BC_beta_div <- ggplot(KO.feces.genmat.tss.bray.segs) + theme_bw() +
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
  geom_text(aes (x= cMDS1, y = cMDS2,label= gen_material), colour= "white", size = 3, fontface = "bold") +
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
  annotate("text", x = -0.3, y = 0.1, ##change coordinates as needed
           label = "Library Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Library Type)
  annotate("text", x = -0.3, y = 0.1, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces_KO_BC_adonis_genmat * 100, 1), "%",
                         "\np = ", round(pvalue_feces_KO_BC_adonis_genmat, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
feces_genmat_KO_BC_beta_div 

####CATCH BASIN SAMPLES (DNA vs cDNA)####
##Subsetting only water samples
KO.water.tss <- subset_samples(KO.tss_ps, sample_type=="Water")
KO.water.tss <- prune_taxa(taxa_sums(KO.water.tss) > 0, KO.water.tss) 
KO.water.tss #24 samples
##Distance matrix
KO.water.tss.bray <- vegdist(t(KO.water.tss@otu_table), method = "bray") 
KO.water.tss.df <- as(KO.water.tss@sam_data,"data.frame") # make DF from metadata

##Adding a column concatenating "gen_material" and "feedlot", then making feedlot a factor
KO.water.tss.df<- KO.water.tss.df %>%
  mutate (gen_material_spec_2 = paste(gen_material, feedlot, sep = '_'),
          feedlot = factor(feedlot))

#####PERMANOVA#####
#Is there an interaction? 
set.seed(87)
water_KO_BC_adonis_interaction <- adonis2(KO.water.tss.bray ~ gen_material * feedlot,
                                       by = "margin",
                                       KO.water.tss.df, 
                                       p.adjust.methods = "BH", permutations = 9999)
water_KO_BC_adonis_interaction #No

##Model Library Type and feedlot, stratified by original sample from which DNA and RNA were extracted 
set.seed(87)
water_KO_BC_adonis <- adonis2(KO.water.tss.bray ~ gen_material + feedlot, 
                           strata = KO.water.tss.df$original_sample,
                           by = "margin",KO.water.tss.df, 
                           p.adjust.methods = "BH", permutations = 9999)
water_KO_BC_adonis #1.5% of the variation is due to Library Type, p = 0.0005
#70.7% of the variation is due to feedlot, p = 0.0005

#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.water.feedlot.disp.KO <- betadisper(KO.water.tss.bray, KO.water.tss.df$feedlot)
bray.water.feedlot.disp.KO
##Then test by permuting
set.seed(87)
bray.water.feedlot.permdisp.KO <- permutest(bray.water.feedlot.disp.KO, permutations = 9999, pairwise = 1)
bray.water.feedlot.permdisp.KO 
##Different dispersions of variance ( p = 0.0043)

#PERMDISP- Gen_material
# Run the betadisper function, average distance to centroid
bray.water.genmat.disp.KO <- betadisper(KO.water.tss.bray, KO.water.tss.df$gen_material)
bray.water.genmat.disp.KO
##Then test by permuting
set.seed(87)
bray.water.genmat.permdisp.KO <- permutest(bray.water.genmat.disp.KO, permutations = 9999)
bray.water.genmat.permdisp.KO 
#No difference in dispersions of variance between DNA and cDNA (p = 0.49)

#### ORDINATION
set.seed(87)
KO.water.tss.bray.ord <- metaMDS(KO.water.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

#Extract ordination coordinates, add metadata
KO.water.tss.bray.plot <- ordiplot(KO.water.tss.bray.ord$points)
#Extract the scores (coordinates) of the points from the ordination plot object 
KO.water.tss.bray.scrs <- scores(KO.water.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
KO.water.tss.bray.scrs <- cbind(as.data.frame(KO.water.tss.bray.scrs), 
                               gen_material = KO.water.tss.df$gen_material, 
                               feedlot = factor(KO.water.tss.df$feedlot), 
                               gen_material_spec_2 = KO.water.tss.df$gen_material_spec_2,
                               sampleID = rownames(KO.water.tss.df))

#####FEEDLOT EFFECT#####
## BC
KO.water.tss.bray.scrs #have coodrinates and metadata
##Centroids according to feedlot
KO.water.feedlot.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = KO.water.tss.bray.scrs, FUN = mean) 
KO.water.feedlot.tss.bray.segs <- merge(KO.water.tss.bray.scrs, 
                                       setNames(KO.water.feedlot.tss.bray.cent, c("feedlot","cMDS1","cMDS2")),
                                       by = 'feedlot', sort = F)

# Extract R2 and p-values
##Model gen_material and feedlot
water_KO_BC_adonis #have the model
R2_water_feedlot_KO_BC_adonis <- water_KO_BC_adonis$R2[2] 
pvalue_water_feedlot_KO_BC_adonis<-water_KO_BC_adonis$`Pr(>F)`[2]

#### PLOT
water_feedlot_KO_BC_beta_div <- ggplot(KO.water.feedlot.tss.bray.segs) + theme_bw() +
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
  guides(color = guide_legend(override.aes = list(size = 7)))+
  annotate("text", x = 0.4, y = 0.6, ##change coordinates as needed
           label = "Feedlot", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (feedlot)
  annotate("text", x = 0.4, y = 0.6, ##change coordinates as needed
           label = paste("R² = ", round(R2_water_feedlot_KO_BC_adonis * 100, 1), "%",
                         "\np = ", round(pvalue_water_feedlot_KO_BC_adonis, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
water_feedlot_KO_BC_beta_div 

#####LIBRARY TYPE EFFECT#####
## BC
KO.water.tss.bray.scrs #have coordinates and metadata
##Centroids according to gen_material (cDNA/DNA)
KO.water.genmat.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ gen_material, 
                                          data = KO.water.tss.bray.scrs, FUN = mean) 
KO.water.genmat.tss.bray.segs <- merge(KO.water.tss.bray.scrs, 
                                      setNames(KO.water.genmat.tss.bray.cent, c("gen_material","cMDS1","cMDS2")), by = 'gen_material', sort = F)

# Extract R2 and p-values for genmat
water_KO_BC_adonis #have the model
R2_water_KO_BC_adonis_genmat <- water_KO_BC_adonis$R2[1] 
pvalue_water_KO_BC_adonis_genmat<-water_KO_BC_adonis$`Pr(>F)`[1]

#### PLOT
water_genmat_KO_BC_beta_div <- ggplot(KO.water.genmat.tss.bray.segs) + theme_bw() +
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
  geom_text(aes (x= cMDS1, y = cMDS2,label= gen_material), colour= "white", size = 3, fontface = "bold") +
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
  annotate("text", x = 0.5, y = 0.6, ##change coordinates as needed
           label = "Library Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Library Type)
  annotate("text", x = 0.5, y = 0.6, ##change coordinates as needed
           label = paste("R² = ", round(R2_water_KO_BC_adonis_genmat * 100, 1), "%",
                         "\np = ", round(pvalue_water_KO_BC_adonis_genmat, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
water_genmat_KO_BC_beta_div 


####ONLY DNA (Feces vs CB)####
##Subsetting only DNA samples
KO.DNA.tss <- subset_samples(KO.tss_ps, gen_material=="DNA")
KO.DNA.tss <- prune_taxa(taxa_sums(KO.DNA.tss) > 0, KO.DNA.tss) 
KO.DNA.tss 

##Distance matrix
KO.DNA.tss.bray <- vegdist(t(KO.DNA.tss@otu_table), method = "bray") 
KO.DNA.tss.df <- as(KO.DNA.tss@sam_data,"data.frame") # make DF from metadata

#### ORDINATION
set.seed(87)
KO.DNA.tss.bray.ord <- metaMDS(KO.DNA.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

##Get ordination coordinates
KO.DNA.tss.bray.plot <- ordiplot(KO.DNA.tss.bray.ord$points)
#Extracts the scores (coordinates) of the points from the ordination plot object:
KO.DNA.tss.bray.scrs <- scores(KO.DNA.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata:
KO.DNA.tss.bray.scrs <- cbind(as.data.frame(KO.DNA.tss.bray.scrs), 
                              sample_type = KO.DNA.tss.df$sample_type, 
                              feedlot = factor(KO.DNA.tss.df$feedlot)) 
#####PERMANOVA########
#Is there an interaction with feedlot?
set.seed(87)
DNA_KO_BC_adonis_interaction  <- adonis2(KO.DNA.tss.bray ~ sample_type*feedlot,
                                         by = "margin",
                                         KO.DNA.tss.df, permutations = 9999)
DNA_KO_BC_adonis_interaction  #Yes

#Modelling sample type and feedlot
set.seed(87)
DNA_KO_BC_adonis_sampletype_feedlot  <- adonis2(KO.DNA.tss.bray ~ sample_type + feedlot, 
                                     #strata = KO.DNA.tss.df$feedlot, 
                                     by = "margin",
                                     KO.DNA.tss.df, 
                                     permutations = 9999)
DNA_KO_BC_adonis_sampletype_feedlot #49.5% of variation is due to sample type (Water vs Feces) p = 1e-04
#2.9% of the variation due to feedlot, p = 0.0160



#PERMDISP - Sample type
# Run the betadisper function, average distance to centroid
bray.DNA.sampletype.disp.KO <- betadisper(KO.DNA.tss.bray, KO.DNA.tss.df$sample_type)
bray.DNA.sampletype.disp.KO
##Then test by permuting
set.seed(87)
bray.DNA.sampletype.permdisp.KO <- permutest(bray.DNA.sampletype.disp.KO, permutations = 9999)
bray.DNA.sampletype.permdisp.KO 
##Feces vs Water p-value 1e-04

#PERMDISP - Feedlot
# Run the betadisper function, average distance to centroid
bray.DNA.feedlot.disp.KO <- betadisper(KO.DNA.tss.bray, KO.DNA.tss.df$feedlot)
bray.DNA.feedlot.disp.KO
##Then test by permuting
set.seed(87)
bray.DNA.feedlot.permdisp.KO <- permutest(bray.DNA.feedlot.disp.KO, permutations = 9999)
bray.DNA.feedlot.permdisp.KO 
##Not significant for feedlots (p = 0.9)

#####ADDING CENTROIDS FOR PLOTTING - SAMPLE TYPE#######
## BC
KO.DNA.tss.bray.cent.sample_type <- aggregate(cbind(MDS1,MDS2) ~ sample_type, data = KO.DNA.tss.bray.scrs, FUN = mean) ##Centroids according to sample type (water and feces)
KO.DNA.tss.bray.segs.sample_type <- merge(KO.DNA.tss.bray.scrs, 
                                          setNames(KO.DNA.tss.bray.cent.sample_type, 
                                                   c("sample_type", "cMDS1","cMDS2")), by = 'sample_type', sort = F) ##add centroids to main scrs dataframe

##Abbreviated version of sample_type for easier plotting
KO.DNA.tss.bray.segs.sample_type  <- KO.DNA.tss.bray.segs.sample_type  %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"))

# Extract R2 and p-values
R2_DNA_adonis_sample_type <- DNA_KO_BC_adonis_sampletype_feedlot$R2[1] 
pvalue_DNA_adonis_sample_type<-  DNA_KO_BC_adonis_sampletype_feedlot$`Pr(>F)`[1]

######SAMPLE TYPE EFFECT#####
DNA_KO_BC_beta_div_spider_sampletype <- ggplot(KO.DNA.tss.bray.segs.sample_type) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "Metagenomic libraries (DNA)", color = "Sample type") +
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
  guides(
    color= guide_legend(override.aes = list(size = 7)))+
  annotate("text", x = 1.2, y = 0.3, ##change coordinates as needed
           label = "Sample type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = 1.2, y = 0.3, ##change coordinates as needed
           label = paste0("R² = ", round(R2_DNA_adonis_sample_type * 100, 1), "%",
                          "\np = ", round(pvalue_DNA_adonis_sample_type, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values 
DNA_KO_BC_beta_div_spider_sampletype

######FEEDLOT EFFECT#######
## BC
KO.DNA.tss.bray.scrs #Already have the ordination coordinates with metadata
KO.DNA.tss.bray.cent.feedlot <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = KO.DNA.tss.bray.scrs, 
                                          FUN = mean) ##Centroids according to feedlot
KO.DNA.tss.bray.segs.feedlot <- merge(KO.DNA.tss.bray.scrs, 
                                      setNames(KO.DNA.tss.bray.cent.feedlot, c("feedlot", "cMDS1","cMDS2")), 
                                      by = 'feedlot', 
                                      sort = F) ##add centroids to main scrs dataframe

# Extract R2 and p-values
R2_DNA_adonis_feedlot <- DNA_KO_BC_adonis_sampletype_feedlot$R2[2] 
pvalue_DNA_adonis_feedlot<-  DNA_KO_BC_adonis_sampletype_feedlot$`Pr(>F)`[2]

#Plot
DNA_KO_BC_beta_div_spider_feedlot <- ggplot(KO.DNA.tss.bray.segs.feedlot) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "Metagenomic libraries (DNA)", 
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
        legend.text = element_text(colour = "black", size = 22, face = "bold"),
        legend.title = element_text(colour = "black", size = 22, face = "bold"),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  guides(
    color= guide_legend(override.aes = list(size = 7))
  )+
  annotate("text", x = 1, y = 0.4, ##change coordinates as needed
           label = "Feedlot",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = 1, y = 0.4, ##change coordinates as needed
           label = paste0("R² = ", round(R2_DNA_adonis_feedlot * 100, 1), "%",
                          "\np = ", round(pvalue_DNA_adonis_feedlot, 4)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values

DNA_KO_BC_beta_div_spider_feedlot

####ONLY cDNA (CB vs feces)#############
##Subsetting only cDNA samples
KO.cDNA.tss <- subset_samples(KO.tss_ps, gen_material=="cDNA")
KO.cDNA.tss <- prune_taxa(taxa_sums(KO.cDNA.tss) > 0, KO.cDNA.tss) 
KO.cDNA.tss 

##Distance matrix
KO.cDNA.tss.bray <- vegdist(t(KO.cDNA.tss@otu_table), method = "bray") 
KO.cDNA.tss.df <- as(KO.cDNA.tss@sam_data,"data.frame") # make DF from metadata

#### ORDINATION
set.seed(87)
KO.cDNA.tss.bray.ord <- metaMDS(KO.cDNA.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

##Get ordination coordinates
KO.cDNA.tss.bray.plot <- ordiplot(KO.cDNA.tss.bray.ord$points)
#Extracts the scores (coordinates) of the points from the ordination plot object:
KO.cDNA.tss.bray.scrs <- scores(KO.cDNA.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata:
KO.cDNA.tss.bray.scrs <- cbind(as.data.frame(KO.cDNA.tss.bray.scrs), 
                                  sample_type = KO.cDNA.tss.df$sample_type, 
                                  feedlot = factor(KO.cDNA.tss.df$feedlot)) 
#####PERMANOVA########
#Is there an interaction with feedlot?
set.seed(87)
cDNA_KO_BC_adonis_interaction  <- adonis2(KO.cDNA.tss.bray ~ sample_type*feedlot,
                                         by = "margin",
                                         KO.cDNA.tss.df, permutations = 9999)
cDNA_KO_BC_adonis_interaction  #Yes

#Sample type and feedlot in model
set.seed(87)
cDNA_KO_BC_adonis_sampletype_feedlot  <- adonis2(KO.cDNA.tss.bray ~ sample_type + feedlot, 
                                        by = "margin",
                                        KO.cDNA.tss.df, 
                                        permutations = 9999)
cDNA_KO_BC_adonis_sampletype_feedlot #47.2% of the varaition is due to sample type, p = 0.001
#3% of the variation due to feedlot


#PERMDISP - Sample type
# Run the betadisper function, average distance to centroid
bray.cDNA.sampletype.disp.KO <- betadisper(KO.cDNA.tss.bray, KO.cDNA.tss.df$sample_type)
bray.cDNA.sampletype.disp.KO
##Then test by permuting
set.seed(87)
bray.cDNA.sampletype.permdisp.KO <- permutest(bray.cDNA.sampletype.disp.KO, permutations = 9999)
bray.cDNA.sampletype.permdisp.KO 
##Feces vs Water p-value 1e-04

#PERMDISP - Feedlot
# Run the betadisper function, average distance to centroid
bray.cDNA.feedlot.disp.KO <- betadisper(KO.cDNA.tss.bray, KO.cDNA.tss.df$feedlot)
bray.cDNA.feedlot.disp.KO
##Then test by permuting
set.seed(87)
bray.cDNA.feedlot.permdisp.KO <- permutest(bray.cDNA.feedlot.disp.KO, permutations = 9999)
bray.cDNA.feedlot.permdisp.KO 
##Not significant for feedlots (p = 0.9)

#####SAMPLE TYPE EFFECT#######
## BC
KO.cDNA.tss.bray.cent.sample_type <- aggregate(cbind(MDS1,MDS2) ~ sample_type, data = KO.cDNA.tss.bray.scrs, FUN = mean) ##Centroids according to sample type (water and feces)
KO.cDNA.tss.bray.segs.sample_type <- merge(KO.cDNA.tss.bray.scrs, 
                                              setNames(KO.cDNA.tss.bray.cent.sample_type, 
                                                       c("sample_type", "cMDS1","cMDS2")), by = 'sample_type', sort = F) ##add centroids to main scrs dataframe

##Abbreviated version of sample_type for easier plotting
KO.cDNA.tss.bray.segs.sample_type  <- KO.cDNA.tss.bray.segs.sample_type  %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"))

# Extract R2 and p-values
R2_cDNA_adonis_sample_type_KO <- cDNA_KO_BC_adonis_sampletype_feedlot$R2[1] 
pvalue_cDNA_adonis_sample_type_KO<-cDNA_KO_BC_adonis_sampletype_feedlot$`Pr(>F)`[1]

#Plot
cDNA_BC_KO_beta_div_spider_sampletype <- ggplot(KO.cDNA.tss.bray.segs.sample_type) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "Metatranscriptomic libraries (RNA (cDNA))", color = "Sample type") +
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
  guides(
    color= guide_legend(override.aes = list(size = 7)))+
  annotate("text", x = 1.2, y = 0.3, ##change coordinates as needed
           label = "Sample type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = 1.2, y = 0.3, ##change coordinates as needed
           label = paste0("R² = ", round(R2_cDNA_adonis_sample_type_KO * 100, 1), "%",
                          "\np = ", round(pvalue_cDNA_adonis_sample_type_KO, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values 
cDNA_BC_KO_beta_div_spider_sampletype

#####FEEDLOT EFFECT#######
## BC
KO.cDNA.tss.bray.scrs #Already have the ordination coordinates with metadata
KO.cDNA.tss.bray.cent.feedlot <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = KO.cDNA.tss.bray.scrs, 
                                              FUN = mean) ##Centroids according to feedlot
KO.cDNA.tss.bray.segs.feedlot <- merge(KO.cDNA.tss.bray.scrs, 
                                          setNames(KO.cDNA.tss.bray.cent.feedlot, c("feedlot", "cMDS1","cMDS2")), 
                                          by = 'feedlot', 
                                          sort = F) ##add centroids to main scrs dataframe

# Extract R2 and p-values
R2_cDNA_adonis_feedlot_KO <- cDNA_KO_BC_adonis_sampletype_feedlot$R2[2] 
pvalue_cDNA_adonis_feedlot_KO<-  cDNA_KO_BC_adonis_sampletype_feedlot$`Pr(>F)`[2]

#Plot
cDNA_BC_KO_beta_div_spider_feedlot <- ggplot(KO.cDNA.tss.bray.segs.feedlot) + theme_bw() +
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
        legend.text = element_text(colour = "black", size = 22, face = "bold"),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  annotate("text", x = 1.65, y = 0.9, ##change coordinates as needed
           label = "Feedlot",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = 1.65, y = 0.9, ##change coordinates as needed
           label = paste0("R² = ", round(R2_cDNA_adonis_feedlot_KO * 100, 1), "%",
                          "\np = ", round(pvalue_cDNA_adonis_feedlot_KO, 4)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") +# Annotate R² and p-values
  guides(
    color= guide_legend(override.aes = list(size = 7))
  )
cDNA_BC_KO_beta_div_spider_feedlot

###PATHWAYS#####
####ONLY FECES#######
##Subsetting only feces samples
pathways.feces.tss <- subset_samples(pathways.tss_ps, sample_type=="Feces")
pathways.feces.tss <- prune_taxa(taxa_sums(pathways.feces.tss) > 0, pathways.feces.tss) 
pathways.feces.tss #655 pathways and 96 samples
##Distance matrix
pathways.feces.tss.bray <- vegdist(t(pathways.feces.tss@otu_table), method = "bray") 
pathways.feces.tss.df <- as(pathways.feces.tss@sam_data,"data.frame") # make DF from metadata

##Adding a column concatenating "gen_material" and "feedlot", then making feedlot a factor
pathways.feces.tss.df<- pathways.feces.tss.df %>%
  mutate (gen_material_spec_2 = paste(gen_material, feedlot, sep = '_'),
          feedlot = factor(feedlot))
#####PERMANOVA#####
#Is there an interaction? 
set.seed(87)
feces_pathways_BC_adonis_interaction <- adonis2(pathways.feces.tss.bray ~ gen_material * feedlot,
                                          by = "margin",
                                          pathways.feces.tss.df, 
                                          p.adjust.methods = "BH", permutations = 9999)
feces_pathways_BC_adonis_interaction #No

##Model Library Type and feedlot, stratified by original sample from which DNA and RNA were extracted 
set.seed(87)
feces_pathways_BC_adonis <- adonis2(pathways.feces.tss.bray ~ gen_material + feedlot, 
                              strata = pathways.feces.tss.df$original_sample,
                              by = "margin",pathways.feces.tss.df, 
                              p.adjust.methods = "BH", permutations = 9999)
feces_pathways_BC_adonis #29.9% of the variation is due to Library Type (p = 1e-04)
#12.8% of the variation is due to feedlot (p = 1e-04)

######SUPPLEMENTARY TABLE 5.15#######
stable5.15 <- data.frame(feces_pathways_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Functional Pathways",
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
stable5.15

#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.feces.feedlot.disp.pathways <- betadisper(pathways.feces.tss.bray, pathways.feces.tss.df$feedlot)
bray.feces.feedlot.disp.pathways
##Then test by permuting
set.seed(87)
bray.feces.feedlot.permdisp.pathways <- permutest(bray.feces.feedlot.disp.pathways, permutations = 9999, pairwise = 1)
bray.feces.feedlot.permdisp.pathways 
##Different variances between feedlots p = 0.0125

#PERMDISP- Gen_material
# Run the betadisper function, average distance to centroid
bray.feces.genmat.disp.pathways <- betadisper(pathways.feces.tss.bray, pathways.feces.tss.df$gen_material)
bray.feces.genmat.disp.pathways
##Then test by permuting
set.seed(87)
bray.feces.genmat.permdisp.pathways <- permutest(bray.feces.genmat.disp.pathways, permutations = 9999)
bray.feces.genmat.permdisp.pathways 
#Different variances between DNA and cDNA p = 0.023

#### ORDINATION
set.seed(87)
pathways.feces.tss.bray.ord <- metaMDS(pathways.feces.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

#Extract ordination coordinates, add metadata
pathways.feces.tss.bray.plot <- ordiplot(pathways.feces.tss.bray.ord$points)
#Extract the scores (coordinates) of the points from the ordination plot object 
pathways.feces.tss.bray.scrs <- scores(pathways.feces.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
pathways.feces.tss.bray.scrs <- cbind(as.data.frame(pathways.feces.tss.bray.scrs), 
                                gen_material = pathways.feces.tss.df$gen_material, 
                                feedlot = factor(pathways.feces.tss.df$feedlot), 
                                gen_material_spec_2 = pathways.feces.tss.df$gen_material_spec_2,
                                sampleID = rownames(pathways.feces.tss.df))
#####FEEDLOT EFFECT#####
## BC
pathways.feces.tss.bray.scrs #have coodrinates and metadata
##Centroids according to feedlot
pathways.feces.feedlot.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = pathways.feces.tss.bray.scrs, FUN = mean) 
pathways.feces.feedlot.tss.bray.segs <- merge(pathways.feces.tss.bray.scrs, 
                                        setNames(pathways.feces.feedlot.tss.bray.cent, c("feedlot","cMDS1","cMDS2")),
                                        by = 'feedlot', sort = F)

# Extract R2 and p-values
##Model gen_material and feedlot
feces_pathways_BC_adonis #have the model
R2_feces_feedlot_BC_adonis_pathways <- feces_pathways_BC_adonis$R2[2] 
pvalue_feces_feedlot_BC_adonis_pathways<-feces_pathways_BC_adonis$`Pr(>F)`[2]

#### PLOT
fecespathways_feedlot_BC_beta_div <- ggplot(pathways.feces.feedlot.tss.bray.segs) + theme_bw() +
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
  guides(color = guide_legend(override.aes = list(size = 7)))+
  annotate("text", x = 0.15, y = 0.1, ##change coordinates as needed
           label = "Feedlot", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (feedlot)
  annotate("text", x = 0.15, y = 0.1, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces_feedlot_BC_adonis_pathways * 100, 1), "%",
                         "\np = ", round(pvalue_feces_feedlot_BC_adonis_pathways, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
fecespathways_feedlot_BC_beta_div 

#####LIBRARY TYPE EFFECT#####
## BC
pathways.feces.tss.bray.scrs #have coordinates and metadata
##Centroids according to gen_material (cDNA/DNA)
pathways.feces.genmat.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ gen_material, 
                                           data = pathways.feces.tss.bray.scrs, FUN = mean) 
pathways.feces.genmat.tss.bray.segs <- merge(pathways.feces.tss.bray.scrs, 
                                       setNames(pathways.feces.genmat.tss.bray.cent, c("gen_material","cMDS1","cMDS2")), by = 'gen_material', sort = F)%>%
  mutate(gen_material_plot = dplyr::recode(gen_material, 
                                           "DNA" = "DNA", 
                                           "cDNA" = "RNA"))

# Extract R2 and p-values for genmat
feces_pathways_BC_adonis #have the model
R2_feces_pathways_BC_adonis_genmat <- feces_pathways_BC_adonis$R2[1] 
pvalue_feces_pathways_BC_adonis_genmat<-feces_pathways_BC_adonis$`Pr(>F)`[1]

#### PLOT
feces_genmat_pathways_BC_beta_div <- ggplot(pathways.feces.genmat.tss.bray.segs) + theme_bw() +
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
  geom_text(aes (x= cMDS1, y = cMDS2,label= gen_material_plot), 
            colour= "white", size = 3, fontface = "bold") +
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
  annotate("text", x = 0.15, y = 0.1, ##change coordinates as needed
           label = "Library Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Library Type)
  annotate("text", x = 0.15, y = 0.1, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces_pathways_BC_adonis_genmat * 100, 1), "%",
                         "\np = ", round(pvalue_feces_pathways_BC_adonis_genmat, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
feces_genmat_pathways_BC_beta_div 

####ONLY CATCH BASINS (DNA vs cDNA)####
##Subsetting only water samples
pathways.water.tss <- subset_samples(pathways.tss_ps, sample_type=="Water")
pathways.water.tss <- prune_taxa(taxa_sums(pathways.water.tss) > 0, pathways.water.tss) 
pathways.water.tss #736 pathways in 24 samples
##Distance matrix
pathways.water.tss.bray <- vegdist(t(pathways.water.tss@otu_table), method = "bray") 
pathways.water.tss.df <- as(pathways.water.tss@sam_data,"data.frame") # make DF from metadata

##Adding a column concatenating "gen_material" and "feedlot", then making feedlot a factor
pathways.water.tss.df<- pathways.water.tss.df %>%
  mutate (gen_material_spec_2 = paste(gen_material, feedlot, sep = '_'),
          feedlot = factor(feedlot))

#####PERMANOVA#####
#Is there an interaction? 
set.seed(87)
water_pathways_BC_adonis_interaction <- adonis2(pathways.water.tss.bray ~ gen_material * feedlot,
                                          by = "margin",
                                          pathways.water.tss.df, 
                                          p.adjust.methods = "BH", permutations = 9999)
water_pathways_BC_adonis_interaction #No

##Model Library Type and feedlot, stratified by original sample from which DNA and RNA were extracted 
set.seed(87)
water_pathways_BC_adonis <- adonis2(pathways.water.tss.bray ~ gen_material + feedlot, 
                              strata = pathways.water.tss.df$original_sample,
                              by = "margin",pathways.water.tss.df, 
                              p.adjust.methods = "BH", permutations = 9999)
water_pathways_BC_adonis #2.5% of the variation is due to Library Type (p=0.0004883 NS)
#63.7% of the variation is due to feedlot (P = 0.0004883)

######SUPPLEMENTARY TABLE 5.16#######
stable5.16 <- data.frame(water_pathways_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Functional Pathways",
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
stable5.16

#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.water.feedlot.disp.pathways <- betadisper(pathways.water.tss.bray, pathways.water.tss.df$feedlot)
bray.water.feedlot.disp.pathways
##Then test by permuting
set.seed(87)
bray.water.feedlot.permdisp.pathways <- permutest(bray.water.feedlot.disp.pathways, permutations = 9999, pairwise = 1)
bray.water.feedlot.permdisp.pathways 
##Different variances between feedlots (p = 0.0389)

#PERMDISP- Library Type
# Run the betadisper function, average distance to centroid
bray.water.genmat.disp.pathways <- betadisper(pathways.water.tss.bray, pathways.water.tss.df$gen_material)
bray.water.genmat.disp.pathways
##Then test by permuting
set.seed(87)
bray.water.genmat.permdisp.pathways <- permutest(bray.water.genmat.disp.pathways, permutations = 9999)
bray.water.genmat.permdisp.pathways 
#No difference in dispersions of variance between DNA and cDNA (p = 0.5579)

#### ORDINATION
set.seed(87)
pathways.water.tss.bray.ord <- metaMDS(pathways.water.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

#Extract ordination coordinates, add metadata
pathways.water.tss.bray.plot <- ordiplot(pathways.water.tss.bray.ord$points)
#Extract the scores (coordinates) of the points from the ordination plot object 
pathways.water.tss.bray.scrs <- scores(pathways.water.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
pathways.water.tss.bray.scrs <- cbind(as.data.frame(pathways.water.tss.bray.scrs), 
                                gen_material = pathways.water.tss.df$gen_material, 
                                feedlot = factor(pathways.water.tss.df$feedlot), 
                                gen_material_spec_2 = pathways.water.tss.df$gen_material_spec_2,
                                sampleID = rownames(pathways.water.tss.df))

#####FEEDLOT EFFECT#####
## BC
pathways.water.tss.bray.scrs #have coodrinates and metadata
##Centroids according to feedlot
pathways.water.feedlot.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = pathways.water.tss.bray.scrs, FUN = mean) 
pathways.water.feedlot.tss.bray.segs <- merge(pathways.water.tss.bray.scrs, 
                                        setNames(pathways.water.feedlot.tss.bray.cent, c("feedlot","cMDS1","cMDS2")),
                                        by = 'feedlot', sort = F)

# Extract R2 and p-values
##Model gen_material and feedlot
water_pathways_BC_adonis #have the model
R2_water_feedlot_pathways_BC_adonis <- water_pathways_BC_adonis$R2[2] 
pvalue_water_feedlot_pathways_BC_adonis<-water_pathways_BC_adonis$`Pr(>F)`[2]

#### PLOT
water_feedlot_pathways_BC_beta_div <- ggplot(pathways.water.feedlot.tss.bray.segs) + theme_bw() +
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
  guides(color = guide_legend(override.aes = list(size = 7)))+
  annotate("text", x = 0.25, y = 0.2, ##change coordinates as needed
           label = "Feedlot", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (feedlot)
  annotate("text", x = 0.25, y = 0.2, ##change coordinates as needed
           label = paste("R² = ", round(R2_water_feedlot_pathways_BC_adonis * 100, 1), "%",
                         "\np = ", round(pvalue_water_feedlot_pathways_BC_adonis, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
water_feedlot_pathways_BC_beta_div 

#####LIBRARY TYPE EFFECT#####
## BC
pathways.water.tss.bray.scrs #have coordinates and metadata
##Centroids according to gen_material (cDNA/DNA)
pathways.water.genmat.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ gen_material, 
                                           data = pathways.water.tss.bray.scrs, FUN = mean) 
pathways.water.genmat.tss.bray.segs <- merge(pathways.water.tss.bray.scrs, 
                                       setNames(pathways.water.genmat.tss.bray.cent, c("gen_material","cMDS1","cMDS2")), by = 'gen_material', sort = F)%>%
  mutate(gen_material_plot = dplyr::recode(gen_material, 
                                           "DNA" = "DNA", 
                                           "cDNA" = "RNA"))

# Extract R2 and p-values for genmat
water_pathways_BC_adonis #have the model
R2_water_pathways_BC_adonis_genmat <- water_pathways_BC_adonis$R2[1] 
pvalue_water_pathways_BC_adonis_genmat<-water_pathways_BC_adonis$`Pr(>F)`[1]

#### PLOT
water_genmat_pathways_BC_beta_div <- ggplot(pathways.water.genmat.tss.bray.segs) + theme_bw() +
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
  geom_text(aes (x= cMDS1, y = cMDS2,label= gen_material_plot), colour= "white", size = 3, fontface = "bold") +
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
  annotate("text", x = 0.25, y = 0.2, ##change coordinates as needed
           label = "Library Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Library Type)
  annotate("text", x = 0.25, y = 0.2, ##change coordinates as needed
           label = paste("R² = ", round(R2_water_pathways_BC_adonis_genmat * 100, 1), "%",
                         "\np = ", round(pvalue_water_pathways_BC_adonis_genmat, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
water_genmat_pathways_BC_beta_div 


####ONLY DNA (Feces vs CB)####
##Subsetting only DNA samples
pathways.DNA.tss <- subset_samples(pathways.tss_ps, gen_material=="DNA")
pathways.DNA.tss <- prune_taxa(taxa_sums(pathways.DNA.tss) > 0, pathways.DNA.tss) 
pathways.DNA.tss 

##Distance matrix
pathways.DNA.tss.bray <- vegdist(t(pathways.DNA.tss@otu_table), method = "bray") 
pathways.DNA.tss.df <- as(pathways.DNA.tss@sam_data,"data.frame")%>% # make DF from metadata
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")))

#### ORDINATION
set.seed(87)
pathways.DNA.tss.bray.ord <- metaMDS(pathways.DNA.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

##Get ordination coordinates
pathways.DNA.tss.bray.plot <- ordiplot(pathways.DNA.tss.bray.ord$points)
#Extracts the scores (coordinates) of the points from the ordination plot object:
pathways.DNA.tss.bray.scrs <- scores(pathways.DNA.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata:
pathways.DNA.tss.bray.scrs <- cbind(as.data.frame(pathways.DNA.tss.bray.scrs), 
                              sample_type = pathways.DNA.tss.df$sample_type, 
                              feedlot = factor(pathways.DNA.tss.df$feedlot)) 
#####PERMANOVA########
#Modelling sample type and feedlot
set.seed(87)
DNA_pathways_BC_adonis_sampletype_feedlot  <- adonis2(pathways.DNA.tss.bray ~ sample_type+feedlot,
                                         by = "margin",
                                         pathways.DNA.tss.df, permutations = 9999)
DNA_pathways_BC_adonis_sampletype_feedlot  #Significant effect of sample type and feedlot


######SUPPLEMENTARY TABLE 5.13 #######
stable5.13 <- data.frame(DNA_pathways_BC_adonis_sampletype_feedlot, check.names = F)%>%
  mutate(Dataset = "Functional Pathways",
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
stable5.13


#PERMDISP - Sample type
# Run the betadisper function, average distance to centroid
bray.DNA.sampletype.disp.pathways <- betadisper(pathways.DNA.tss.bray, pathways.DNA.tss.df$sample_type)
bray.DNA.sampletype.disp.pathways
##Then test by permuting
set.seed(87)
bray.DNA.sampletype.permdisp.pathways <- permutest(bray.DNA.sampletype.disp.pathways, permutations = 9999)
bray.DNA.sampletype.permdisp.pathways 
##Feces vs Water dispersions p-value 1e-04

#PERMDISP - Feedlot
# Run the betadisper function, average distance to centroid
bray.DNA.feedlot.disp.pathways <- betadisper(pathways.DNA.tss.bray, pathways.DNA.tss.df$feedlot)
bray.DNA.feedlot.disp.pathways
##Then test by permuting
set.seed(87)
bray.DNA.feedlot.permdisp.pathways <- permutest(bray.DNA.feedlot.disp.pathways, permutations = 9999)
bray.DNA.feedlot.permdisp.pathways 
##Not significant for feedlots (p = 0.9)

#####SAMPLE TYPE EFFECT#######
## BC
pathways.DNA.tss.bray.cent.sample_type <- aggregate(cbind(MDS1,MDS2) ~ sample_type, data = pathways.DNA.tss.bray.scrs, FUN = mean) ##Centroids according to sample type (water and feces)
pathways.DNA.tss.bray.segs.sample_type <- merge(pathways.DNA.tss.bray.scrs, 
                                          setNames(pathways.DNA.tss.bray.cent.sample_type, 
                                                   c("sample_type", "cMDS1","cMDS2")), by = 'sample_type', sort = F) ##add centroids to main scrs dataframe

##Abbreviated version of sample_type for easier plotting
pathways.DNA.tss.bray.segs.sample_type  <- pathways.DNA.tss.bray.segs.sample_type  %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"))

# Extract R2 and p-values
R2_DNA_adonis_sample_type <- DNA_pathways_BC_adonis_sampletype_feedlot$R2[1] 
pvalue_DNA_adonis_sample_type<-  DNA_pathways_BC_adonis_sampletype_feedlot$`Pr(>F)`[1]

###PLOT
DNA_pathways_BC_beta_div_spider_sampletype <- ggplot(pathways.DNA.tss.bray.segs.sample_type) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "Metagenomic libraries (DNA)", color = "Sample type") +
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
  guides(
    color= guide_legend(override.aes = list(size = 7)))+
  annotate("text", x = 0.15, y = 0.12, ##change coordinates as needed
           label = "Sample Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = 0.15, y = 0.12, ##change coordinates as needed
           label = paste0("R² = ", round(R2_DNA_adonis_sample_type * 100, 1), "%",
                          "\np = ", round(pvalue_DNA_adonis_sample_type, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values 
DNA_pathways_BC_beta_div_spider_sampletype


#####FEEDLOT EFFECT#######
## BC
pathways.DNA.tss.bray.scrs #Already have the ordination coordinates with metadata
pathways.DNA.tss.bray.cent.feedlot <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = pathways.DNA.tss.bray.scrs, 
                                          FUN = mean) ##Centroids according to feedlot
pathways.DNA.tss.bray.segs.feedlot <- merge(pathways.DNA.tss.bray.scrs, 
                                      setNames(pathways.DNA.tss.bray.cent.feedlot, c("feedlot", "cMDS1","cMDS2")), 
                                      by = 'feedlot', 
                                      sort = F) ##add centroids to main scrs dataframe

# Extract R2 and p-values
R2_DNA_adonis_feedlot <- DNA_pathways_BC_adonis_sampletype_feedlot$R2[2] 
pvalue_DNA_adonis_feedlot<-  DNA_pathways_BC_adonis_sampletype_feedlot$`Pr(>F)`[2]

###PLOT
DNA_pathways_BC_beta_div_spider_feedlot <- ggplot(pathways.DNA.tss.bray.segs.feedlot) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "Metagenomic libraries (DNA)", 
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
        legend.text = element_text(colour = "black", size = 22, face = "bold"),
        legend.title = element_text(colour = "black", size = 22, face = "bold"),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  guides(
    color= guide_legend(override.aes = list(size = 7))
  )+
  annotate("text",  x = -0.05, y = 0.12, ##change coordinates as needed
           label = "Feedlot",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = -0.05, y = 0.12, ##change coordinates as needed
           label = paste0("R² = ", round(R2_DNA_adonis_feedlot * 100, 1), "%",
                          "\np = ", round(pvalue_DNA_adonis_feedlot, 3)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
DNA_pathways_BC_beta_div_spider_feedlot

####ONLY cDNA samples (Feces vs CB)#############
##Subsetting only cDNA samples
pathways.cDNA.tss <- subset_samples(pathways.tss_ps, gen_material=="cDNA")
pathways.cDNA.tss <- prune_taxa(taxa_sums(pathways.cDNA.tss) > 0, pathways.cDNA.tss) 
pathways.cDNA.tss 

##Distance matrix
pathways.cDNA.tss.bray <- vegdist(t(pathways.cDNA.tss@otu_table), method = "bray") 
pathways.cDNA.tss.df <- as(pathways.cDNA.tss@sam_data,"data.frame") %>% # make DF from metadata
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")))

#### ORDINATION
set.seed(87)
pathways.cDNA.tss.bray.ord <- metaMDS(pathways.cDNA.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

##Get ordination coordinates
pathways.cDNA.tss.bray.plot <- ordiplot(pathways.cDNA.tss.bray.ord$points)
#Extracts the scores (coordinates) of the points from the ordination plot object:
pathways.cDNA.tss.bray.scrs <- scores(pathways.cDNA.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata:
pathways.cDNA.tss.bray.scrs <- cbind(as.data.frame(pathways.cDNA.tss.bray.scrs), 
                               sample_type = pathways.cDNA.tss.df$sample_type, 
                               feedlot = factor(pathways.cDNA.tss.df$feedlot)) 
#####PERMANOVA########
#Modelling sample type and feedlot
set.seed(87)
cDNA_pathways_BC_adonis_sampletype_feedlot  <- adonis2(pathways.cDNA.tss.bray ~ sample_type+feedlot,
                                                      by = "margin",
                                                      pathways.cDNA.tss.df, permutations = 9999)
cDNA_pathways_BC_adonis_sampletype_feedlot  #Significant effect of sample type and feedlot

######SUPPLEMENTARY TABLE 5.14#######
stable5.14 <- data.frame(cDNA_pathways_BC_adonis_sampletype_feedlot, check.names = F)%>%
  mutate(Dataset = "Functional Pathways",
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
stable5.14 


#PERMDISP - Sample type
# Run the betadisper function, average distance to centroid
bray.cDNA.sampletype.disp.pathways <- betadisper(pathways.cDNA.tss.bray, pathways.cDNA.tss.df$sample_type)
bray.cDNA.sampletype.disp.pathways
##Then test by permuting
set.seed(87)
bray.cDNA.sampletype.permdisp.pathways <- permutest(bray.cDNA.sampletype.disp.pathways, permutations = 9999)
bray.cDNA.sampletype.permdisp.pathways 
##Feces vs Water p-value 1e-04 (different variances)

#PERMDISP - Feedlot
# Run the betadisper function, average distance to centroid
bray.cDNA.feedlot.disp.pathways <- betadisper(pathways.cDNA.tss.bray, pathways.cDNA.tss.df$feedlot)
bray.cDNA.feedlot.disp.pathways
##Then test by permuting
set.seed(87)
bray.cDNA.feedlot.permdisp.pathways <- permutest(bray.cDNA.feedlot.disp.pathways, permutations = 9999)
bray.cDNA.feedlot.permdisp.pathways 
##Not significant for feedlots (p = 0.9)

#####SAMPLE TYPE EFFECT#######
## BC
pathways.cDNA.tss.bray.cent.sample_type <- aggregate(cbind(MDS1,MDS2) ~ sample_type, data = pathways.cDNA.tss.bray.scrs, FUN = mean) ##Centroids according to sample type (water and feces)
pathways.cDNA.tss.bray.segs.sample_type <- merge(pathways.cDNA.tss.bray.scrs, 
                                           setNames(pathways.cDNA.tss.bray.cent.sample_type, 
                                                    c("sample_type", "cMDS1","cMDS2")), by = 'sample_type', sort = F) ##add centroids to main scrs dataframe

##Abbreviated version of sample_type for easier plotting
pathways.cDNA.tss.bray.segs.sample_type  <- pathways.cDNA.tss.bray.segs.sample_type  %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"))

# Extract R2 and p-values
R2_cDNA_adonis_sample_type_pathways <- cDNA_pathways_BC_adonis_sampletype_feedlot$R2[1] 
pvalue_cDNA_adonis_sample_type_pathways<-cDNA_pathways_BC_adonis_sampletype_feedlot$`Pr(>F)`[1]

###PLOT
cDNA_BC_pathways_beta_div_spider_sampletype <- ggplot(pathways.cDNA.tss.bray.segs.sample_type) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "Metatranscriptomic libraries (RNA (cDNA))", color = "Sample type") +
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
  guides(
    color= guide_legend(override.aes = list(size = 7)))+
  annotate("text", x = 0.15, y = 0.1, ##change coordinates as needed
           label = "Sample Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = 0.15, y = 0.1, ##change coordinates as needed
           label = paste0("R² = ", round(R2_cDNA_adonis_sample_type_pathways * 100, 1), "%",
                          "\np = ", round(pvalue_cDNA_adonis_sample_type_pathways, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values 
cDNA_BC_pathways_beta_div_spider_sampletype

#####FEEDLOT EFFECT#######
## BC
pathways.cDNA.tss.bray.scrs #Already have the ordination coordinates with metadata
pathways.cDNA.tss.bray.cent.feedlot <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = pathways.cDNA.tss.bray.scrs, 
                                           FUN = mean) ##Centroids according to feedlot
pathways.cDNA.tss.bray.segs.feedlot <- merge(pathways.cDNA.tss.bray.scrs, 
                                       setNames(pathways.cDNA.tss.bray.cent.feedlot, c("feedlot", "cMDS1","cMDS2")), 
                                       by = 'feedlot', 
                                       sort = F) ##add centroids to main scrs dataframe

# Extract R2 and p-values
R2_cDNA_adonis_feedlot_pathways <- cDNA_pathways_BC_adonis_sampletype_feedlot$R2[2] 
pvalue_cDNA_adonis_feedlot_pathways<-  cDNA_pathways_BC_adonis_sampletype_feedlot$`Pr(>F)`[2]

###PLOT
cDNA_BC_pathways_beta_div_spider_feedlot <- ggplot(pathways.cDNA.tss.bray.segs.feedlot) + theme_bw() +
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
        legend.text = element_text(colour = "black", size = 22, face = "bold"),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  annotate("text", x = -0.15, y = 0.13, ##change coordinates as needed
           label = "Feedlot",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = -0.15, y = 0.13, ##change coordinates as needed
           label = paste0("R² = ", round(R2_cDNA_adonis_feedlot_pathways * 100, 1), "%",
                          "\np = ", round(pvalue_cDNA_adonis_feedlot_pathways, 3)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") + # Annotate R² and p-values
  guides(
    color= guide_legend(override.aes = list(size = 7))
  )
cDNA_BC_pathways_beta_div_spider_feedlot

####FIGURE 4GH - FECES VS WATER EMPHASIS#####
cDNA_BC_pathways_beta_div_spider_sampletype
DNA_pathways_BC_beta_div_spider_sampletype

#cDNA and DNA
cDNAandDNA_BC_pathways_beta_div_spider_sampletype <- ggarrange(
  DNA_pathways_BC_beta_div_spider_sampletype + 
    theme(
          # plot.title = element_text(size = 35),
          plot.title = element_blank(),
          axis.title = element_text(size = 20)),
  cDNA_BC_pathways_beta_div_spider_sampletype +
    theme(      
      # plot.title = element_text(size = 35),
      plot.title = element_blank(),
      axis.title = element_text(size = 20)),
  align = "h",
  nrow = 1,
  labels = c("G", "H"),
  font.label = list(size = 30),
  common.legend = TRUE,
  legend = "none")
cDNAandDNA_BC_pathways_beta_div_spider_sampletype

##DNA and cDNA ordination plots - Figure 4GH
figure4GH <- cDNAandDNA_BC_pathways_beta_div_spider_sampletype
figure4GH
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure4GH.svg", 
       plot = figure4GH, 
       device = "svg", 
       #dpi = 600,
       width = 30, 
       height = 4.6, 
       bg = "white")


####FIGURE 7GH - EMPHASIS ON LIBRARY TYPE #######
water_genmat_pathways_BC_beta_div 
feces_genmat_pathways_BC_beta_div 

WaterandFeces_pathways_beta_div_spider_genmat <- ggarrange(
  feces_genmat_pathways_BC_beta_div + 
    theme(
      # plot.title = element_text(size = 35),
      plot.title = element_blank(),
      axis.title = element_text(size = 20)),
  water_genmat_pathways_BC_beta_div + 
    theme(
      # plot.title = element_text(size = 35),
      plot.title = element_blank(),
      axis.title = element_text(size = 20)),
  align = "h",
  nrow = 1,
  common.legend = TRUE,
  labels = c("G", "H"),
  font.label = list(size = 30),
  legend = "none")
WaterandFeces_pathways_beta_div_spider_genmat

##DNA and cDNA ordination plots - Figure 7GH
figure7GH <- WaterandFeces_pathways_beta_div_spider_genmat
figure7GH
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure7GH.svg", 
       plot = figure7GH, 
       device = "svg", 
       dpi = 600, 
       width = 30, 
       height = 4.6, 
       bg = "white")

####SUPPLEMENTARY FIGURE 6EF#######
#####Effect of feedlot On DNA and cDNA#####
sfigure6EandF <- ggarrange(DNA_pathways_BC_beta_div_spider_feedlot+
                             theme(plot.title = element_blank()), 
                           cDNA_BC_pathways_beta_div_spider_feedlot+
                             theme(plot.title = element_blank()),
                           labels = c("E", "F"), 
                           font.label = list(size = 22),
                           legend = "none",
                           common.legend = T)
  # labs(title = "FUNCTIONAL PATHWAYS")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
sfigure6EandF
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryFigure6EF.svg", 
       plot = sfigure6EandF, 
       device = "svg", 
       width = 18, 
       height = 6, 
       #dpi = 600, 
       bg = "white")

####SUPPLEMENTARY FIGURE 7EF#######
#####Effect of feedlot On Feces and CB#####
fecespathways_feedlot_BC_beta_div 
water_feedlot_pathways_BC_beta_div 

sfigure7EandF <- ggarrange(fecespathways_feedlot_BC_beta_div+
                             theme(plot.title = element_blank()),
                           water_feedlot_pathways_BC_beta_div+
                             theme(plot.title = element_blank()),
                           labels = c("E", "F"), 
                           font.label = list(size = 22),
                           legend = "bottom",
                           common.legend = T)
  # labs(title = "FUNCTIONAL PATHWAYS")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
sfigure7EandF
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryFigure7EF.svg", 
       plot = sfigure7EandF , 
       device = "svg",
       bg = "white", 
       dpi = 600,
       width = 18, 
       height = 6)


###SUPPLEMENTARY TABLE 5 - PATHWAYS SECTION#######
stable5_pathways <- bind_rows(stable5.13, 
                              stable5.14,
                              stable5.15,
                              stable5.16)%>%  
  select(Dataset, `Library Type`, `Sample Type`, `Fixed Effect`, Df, SumOfSqs, R2, `F`, `Pr(>F)`)
stable5_pathways
write_xlsx(stable5_pathways, 
           "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryTable5_3.xlsx")

#DIFFERENTIAL ABUNDANCE - Pathways####
##DNA (Feces vs Water)#######
##Getting counts in DNA samples 
pathways.DNA <- subset_samples(pathways_ps_filt, gen_material == "DNA") 
pathways.DNA <- prune_taxa(taxa_sums(pathways.DNA) > 0, pathways.DNA) 
pathways.DNA #726 taxa and 60 samples
ancom_pathways_DNA.counts <- pathways.DNA
ancom_pathways_DNA.counts@sam_data$sample_type <- factor(ancom_pathways_DNA.counts@sam_data$sample_type, levels = c("Feces", "Water"))##reorder sample_type as factor, Feces as "reference"
ancom_pathways_DNA.counts@sam_data$feedlot <- factor(ancom_pathways_DNA.counts@sam_data$feedlot, levels = c("1", "2", "3" ,"4 ","5")) ##make feedlot a factor (not num)

#Relative abundances for DNA samples 
pathways.DNA.tss <- subset_samples(pathways.tss_ps , gen_material == "DNA") 
pathways.DNA.tss <- prune_taxa(taxa_sums(pathways.DNA.tss) > 0, pathways.DNA.tss) 
pathways.DNA.tss #726 taxa and 60 samples

##Filtering out the low relative abundance (less than 0.5%) pathways (RA object sums up to 1, so 0.005 is 0.5%)
pathways.DNA.tss_filt <- filter_taxa(pathways.DNA.tss, function(x) mean(x) > 0.005, TRUE) 
pathways.DNA.tss_filt ## 66 taxa with mean RA > 0.5% across 60 samples (DNA samples) 
##Filtering those Ko groups (> 0.5% RA) on the "raw" counts phyloseq object for DNA
ancom_pathways_DNA.counts_filt <- prune_taxa(taxa_names(pathways.DNA.tss_filt), ancom_pathways_DNA.counts)
ancom_pathways_DNA.counts_filt ##66 taxa, 60 samples

###ANCOMBC####
##running ancombc on the variable of interest (sample_type)
ancombc_output_DNA.pathways <-ancombc2(data= ancom_pathways_DNA.counts_filt, 
                                 assay_name = "counts", 
                                 tax_level = NULL,
                                 fix_formula = "sample_type + feedlot",
                                 # fix_formula = "sample_type",
                                 # rand_formula =  "(1 | feedlot)",
                                 prv_cut = 0.05, 
                                 lib_cut = 0, 
                                 group= "sample_type", 
                                 struc_zero = TRUE, 
                                 neg_lb = TRUE,
                                 alpha = 0.05, #default significance
                                 n_cl = 1, verbose = TRUE)

## extract results from comparisons 
ancombc_output_DNA.pathways_res <- ancombc_output_DNA.pathways$res %>%
  select(-matches("feedlot")) #not doing feedlot comparisons

#Pivot longer the results
ancombc_output_DNA.pathways_pivot <- ancombc_output_DNA.pathways_res %>%
  mutate(across(starts_with("lfc_sample_type"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancombc_output_DNA.pathways_pivot

##Getting rid of _rounded suffix using sub command
ancombc_output_DNA.pathways_pivot$group ##want to get rid of "_rounded"
ancombc_output_DNA.pathways_pivot$group<- sub("_rounded", "", ancombc_output_DNA.pathways_pivot$group) 
ancombc_output_DNA.pathways_pivot$group #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancombc_output_DNA.pathways_pivot <- ancombc_output_DNA.pathways_pivot %>%
  mutate(group= case_when(
    group == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ group ##keeps original name for groups not specified (DNA)
  ))
ancombc_output_DNA.pathways_pivot$group ##Now the group names are shorter and more manageable
ancombc_output_DNA.pathways_pivot<- ancombc_output_DNA.pathways_pivot %>%
  rename(pathway= taxon) ##This ancombc was done at the genus level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancombc_output_DNA.pathways_pivot_2 <- ancombc_output_DNA.pathways_pivot %>%
  mutate(lower.ci = lfc_sample_typeWater - 1.96*se_sample_typeWater,
         upper.ci = lfc_sample_typeWater + 1.96*se_sample_typeWater)

##Final fix - up to make compatible with plotting
ancombc_output_DNA.pathways_pivot_3 <- ancombc_output_DNA.pathways_pivot_2 %>%
  # filter (passed_ss_sample_typeWater == 1) %>% ##Only want those that passed sensitivity testing
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
  select(pathway, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
nrow(ancombc_output_DNA.pathways_pivot_3) ##6 DA pathways between feces and water with ANCOM


###MaAsLin3#######
#Data (otu relative abundances) and metadata for MaAslin
#Relative abundances for DNA samples 
pathways.DNA.tss #726 taxa and 60 samples

##Filtering out the low relative abundance (less than 0.4%) pathways groups (RA object sums up to 1, so 0.004 is 0.4%)
pathways.DNA.tss_filt ## 66 taxa with mean RA > 0.5% across 60 samples (DNA samples) 
maaslin_pathways_DNA.tss <- pathways.DNA.tss_filt
maaslin_pathways_DNA.tss@sam_data$gen_material <- factor(maaslin_pathways_DNA.tss@sam_data$sample_type, levels = c("Feces", "Water"))##reorder sample_type as factor, Feces as "reference"
maaslin_pathways_DNA.tss@sam_data$feedlot <- factor(maaslin_pathways_DNA.tss@sam_data$feedlot, levels = c("1", "2", "3", "4", "5")) ##make feedlot a factor (not num)

##OTU RA counts for Maaslin as a data frame 
data_maaslin_DNA_pathways_tss <- data.frame(otu_table(maaslin_pathways_DNA.tss), check.names = F)

##Sample metadata
metadata_maaslin_DNA <- data.frame(sample_data(maaslin_pathways_DNA.tss), check.names = F) 

##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_DNA_pathways <-  maaslin3(
  input_data = data_maaslin_DNA_pathways_tss, 
  input_metadata = metadata_maaslin_DNA, 
  output = "MaAsLin3_DNA",  
  fixed_effect = c("sample_type", "feedlot"),
  # fixed_effect = c("sample_type"), 
  # random_effects = "feedlot",
  min_prevalence=0.05,
  median_comparison_abundance = T, # default
  median_comparison_prevalence = FALSE, #default 
  min_abundance = 0, ##Input_data has already been filtered
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_DNA_pathways$fit_data_abundance$results) #Abundance results from MaAslin - 330

#Calculate confidence intervals, add taxonomy
maaslin_DNA_pathways_2 <- maaslin_DNA_pathways$fit_data_abundance$results %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) 
#left_join(input_taxonomy_DNA, by = "feature")

##Final edits to put together for plot
maaslin_DNA_pathways_3 <- maaslin_DNA_pathways_2%>%
  #filter(qval_individual < 0.05) %>% ##Only significant taxa.qpval_individual is the corrected q-value of the individual association (only abundance. If it were qval_joint it'd be joint prevalence and abundance association)
  rename(group = name,
         pval = pval_individual,
         qval = qval_individual,
         pathway = feature)%>% ##Renaming
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
  select(pathway, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancombc_output_DNA.pathways_pivot_3)
nrow(maaslin_DNA_pathways_3) ##37 DA pathways between DNA and DNA by MaAslin

##ANCOM and MaAslin together
DA_DNA_pathways_MaAslinANCOM.data <- rbind(ancombc_output_DNA.pathways_pivot_3, maaslin_DNA_pathways_3) %>%
  filter(pathway %in% intersect(maaslin_DNA_pathways_3$pathway,
                                  ancombc_output_DNA.pathways_pivot_3$pathway)) ##Only going to plot those taxa DA by both tests

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_DNA.pathways <- maaslin_DNA_pathways$transformed_data %>% #transformed data is ra transformed otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "pathway", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  #left_join(input_taxonomy_DNA, by = "feature")%>%
  filter(pathway %in% intersect(maaslin_DNA_pathways_3$pathway,
                                  ancombc_output_DNA.pathways_pivot_3$pathway))%>%
  left_join(metadata_maaslin_DNA%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(pathway, Sample, logvalue, sample_type, plot)

#Bias-corrected abundances (ANCOM)
DNA.pathways_log_corr_abn <- ancombc_output_DNA.pathways$bias_correct_log_table %>%
  data.frame(check.names = F)%>% ##make into data frame
  rownames_to_column("pathway")%>%
  filter(pathway %in% intersect(maaslin_DNA_pathways_3$pathway,
                                  ancombc_output_DNA.pathways_pivot_3$pathway))%>% #keep only those pathways in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -pathway, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         sample_type = ifelse(grepl("W", Sample), "Water", "Feces"))

#put together objects to plot DA
DA_DNA_pathways_together <- bind_rows(DA_DNA_pathways_MaAslinANCOM.data, 
                                      DNA.pathways_log_corr_abn, 
                                      RA_MaaslinAncom_DNA.pathways) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_DNA_pathways_together$plot <- factor(DA_DNA_pathways_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_DNA_pathways_together$sample_type <- factor(DA_DNA_pathways_together$sample_type, levels = c("Feces", "Water"))

####Which are the biggest positive fold changes?######
top_positive_fold_changes_DNA <- DA_DNA_pathways_together %>%
  filter(!is.na(test)) %>%
  # get the max coef per pathway for ordering
  group_by(pathway) %>%
  mutate(max_coef = max(coef, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_coef), pathway, test) %>%
  head(n = 20) #top 3 
top_positive_fold_changes_DNA 


####Which are the biggest negative fold changes?#######
top_negative_fold_changes_DNA <- DA_DNA_pathways_together %>%
  filter(!is.na(test)) %>%
  # get the max coef per pathway for ordering
  group_by(pathway) %>%
  mutate(max_coef = max(coef, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_coef), pathway, test) %>%
  tail(n = 8) #top 4
top_negative_fold_changes_DNA 

###PLOTTING DA#
#Plotting 
DA_pathways_DNA_plot_MaAslinANCOM <-
  ggplot(data = DA_DNA_pathways_together) +
  facet_wrap(~ plot, scales='free_x',
             nrow = 1,
             strip.position = "top") +
  geom_boxplot(data=DA_DNA_pathways_together%>%filter(grepl("abundances", plot)),
               aes(x=pathway, y=logvalue, 
                   fill = sample_type, color = sample_type),
               notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(data=DA_DNA_pathways_together%>%filter(grepl("abundances", plot)),
             aes(x=pathway, y=logvalue, fill = sample_type, color = sample_type),
             size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  # facet_nested(. ~ plot, scales='free_x',
  #              space='free_y',
  #              switch='y',
  #              strip=strip_nested(text_y=list(element_text(angle=0))),
  #              labeller=labeller(group=label_wrap_gen(width=10),
  #                                sub_group=label_wrap_gen(width=10))) +
  # ggplot(data=DA_DNA_pathways_together%>%filter(grepl("abundances", plot)),
  #        aes(x=pathway, y=logvalue, fill = sample_type, color = sample_type)) +
  # geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  # geom_point(size = 1, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order = 1,title="Sample Type", title.position="top"),
         fill=guide_legend(order = 1,title="Sample Type", title.position="top"))+
  scale_fill_manual(values=sample.type.palette, 
                    labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  scale_color_manual(values = sample.type.palette, 
                     labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  new_scale_color()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_DNA_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=pathway, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_DNA_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=pathway, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_DNA_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed'
             , alpha=0.5) +
  # scale_x_discrete(position='bottom',
  #                  labels = function(x) str_wrap(x, width = 60))+
  scale_x_discrete(position='bottom', labels = function(x) str_wrap(x, width = 40))+
  scale_y_continuous(position='right') +  # Default to break_labels otherwise 
  coord_flip() +
  scale_color_manual(values=c("red", "blue")) +
  scale_linetype_manual(values=c("11", "solid")) +
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Sample Type", title.position="top"),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 1.5)),
         linetype = guide_legend(title = "Fold change source", title.position = "top",
                                 override.aes = list(linewidth = 1),
                                 theme = theme(legend.key.width = unit(1.5, "cm"))),
         pch=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  labs(title = "Metagenomic libraries (DNA)", x = "Functional Pathway")+
  theme_bw()+
  theme(legend.position="top", 
        legend.key=element_blank(),
        legend.title=element_text(size=15), legend.text=element_text(size=14),
        plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y=element_text(size= 20, angle=0, vjust= 1.045, face = "bold"), 
        axis.text.y=element_text(size=18, vjust = 0.5),
        strip.text=element_text(size=13, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_pathways_DNA_plot_MaAslinANCOM

##cDNA (Feces vs CB)#######
##Getting counts in cDNA samples 
pathways.cDNA <- subset_samples(pathways_ps_filt, gen_material == "cDNA") 
pathways.cDNA <- prune_taxa(taxa_sums(pathways.cDNA) > 0, pathways.cDNA) 
pathways.cDNA #733 taxa and 60 samples
ancom_pathways_cDNA.counts <- pathways.cDNA
ancom_pathways_cDNA.counts@sam_data$sample_type <- factor(ancom_pathways_cDNA.counts@sam_data$sample_type, levels = c("Feces", "Water"))##reorder sample_type as factor, Feces as "reference"
ancom_pathways_cDNA.counts@sam_data$feedlot <- factor(ancom_pathways_cDNA.counts@sam_data$feedlot, levels = c("1", "2", "3" ,"4 ","5")) ##make feedlot a factor (not num)

#Relative abundances for cDNA samples 
pathways.cDNA.tss <- subset_samples(pathways.tss_ps, gen_material == "cDNA") 
pathways.cDNA.tss <- prune_taxa(taxa_sums(pathways.cDNA.tss) > 0, pathways.cDNA.tss) 
pathways.cDNA.tss #733 taxa and 60 samples

##Filtering out the low relative abundance (less than 0.5 %) pathways (RA object sums up to 1, so 0.005 is 0.5%)
pathways.cDNA.tss_filt <- filter_taxa(pathways.cDNA.tss, function(x) mean(x) > 0.005, TRUE) 
pathways.cDNA.tss_filt ## 66 taxa with mean RA > 0.5% across 60 samples (cDNA samples) 
##Filtering those pathways (> 0.5% RA) on the "raw" counts phyloseq object for cDNA
ancom_pathways_cDNA.counts_filt <- prune_taxa(taxa_names(pathways.cDNA.tss_filt), ancom_pathways_cDNA.counts)
ancom_pathways_cDNA.counts_filt ##66 taxa, 60 samples

###ANCOMBC####
##running ancombc on the variable of interest (sample_type)
ancombc_output_cDNA.pathways <-ancombc2(data= ancom_pathways_cDNA.counts_filt, 
                                  assay_name = "counts", 
                                  tax_level = NULL,
                                  fix_formula = "sample_type + feedlot", 
                                  # fix_formula = "sample_type", 
                                  # rand_formula =  "(1 | feedlot)", 
                                  prv_cut = 0.05, 
                                  lib_cut = 0, 
                                  group= "sample_type", 
                                  struc_zero = TRUE, 
                                  neg_lb = TRUE,
                                  alpha = 0.05, #default significance
                                  n_cl = 1, verbose = TRUE)

## extract results from comparisons 
ancombc_output_cDNA.pathways_res <- ancombc_output_cDNA.pathways$res %>%
  select(-matches("feedlot"))

#Pivot longer the results
ancombc_output_cDNA.pathways_pivot <- ancombc_output_cDNA.pathways_res %>%
  mutate(across(starts_with("lfc_sample_type"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancombc_output_cDNA.pathways_pivot

##Getting rid of _rounded suffix using sub command
ancombc_output_cDNA.pathways_pivot$group ##want to get rid of "_rounded"
ancombc_output_cDNA.pathways_pivot$group<- sub("_rounded", "", ancombc_output_cDNA.pathways_pivot$group) 
ancombc_output_cDNA.pathways_pivot$group #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancombc_output_cDNA.pathways_pivot <- ancombc_output_cDNA.pathways_pivot %>%
  mutate(group= case_when(
    group == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ group ##keeps original name for groups not specified (DNA)
  ))
ancombc_output_cDNA.pathways_pivot$group ##Now the group names are shorter and more manageable
ancombc_output_cDNA.pathways_pivot<- ancombc_output_cDNA.pathways_pivot %>%
  rename(pathway = taxon) ##This ancombc was done at the genus level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancombc_output_cDNA.pathways_pivot_2 <- ancombc_output_cDNA.pathways_pivot %>%
  mutate(lower.ci = lfc_sample_typeWater - 1.96*se_sample_typeWater,
         upper.ci = lfc_sample_typeWater + 1.96*se_sample_typeWater)

##Final fix - up to make compatible with plotting
ancombc_output_cDNA.pathways_pivot_3 <- ancombc_output_cDNA.pathways_pivot_2 %>%
  filter (passed_ss_sample_typeWater == 1) %>%
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
  select(pathway, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
nrow(ancombc_output_cDNA.pathways_pivot_3) ##2 DA pathways between feces and water with ANCOM


###MaAsLin3#######
#Data (otu relative abundances) and metadata for MaAslin
#Relative abundances for cDNA samples 
pathways.cDNA.tss #726 taxa and 60 samples

##Filtering out the low relative abundance (less than 0.5 %) pathways (RA object sums up to 1, so 0.005 is 0.5%)
pathways.cDNA.tss_filt ## 65 taxa with mean RA > 0.5% across 60 samples (cDNA samples) 
maaslin_pathways_cDNA.tss <- pathways.cDNA.tss_filt
maaslin_pathways_cDNA.tss@sam_data$gen_material <- factor(maaslin_pathways_cDNA.tss@sam_data$sample_type, levels = c("Feces", "Water"))##reorder sample_type as factor, Feces as "reference"
maaslin_pathways_cDNA.tss@sam_data$feedlot <- factor(maaslin_pathways_cDNA.tss@sam_data$feedlot, levels = c("1", "2", "3", "4", "5")) ##make feedlot a factor (not num)

##OTU RA counts for Maaslin as a data frame 
data_maaslin_cDNA_pathways_tss <- data.frame(otu_table(maaslin_pathways_cDNA.tss), check.names = F)

##Sample metadata
metadata_maaslin_cDNA <- data.frame(sample_data(maaslin_pathways_cDNA.tss), check.names = F) 


##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_cDNA_pathways <-  maaslin3(
  input_data = data_maaslin_cDNA_pathways_tss, 
  input_metadata = metadata_maaslin_cDNA, 
  output = "MaAsLin3_cDNA",  
  fixed_effect = c("sample_type", "feedlot"),
  # fixed_effect = c("sample_type"), 
  # random_effects = "feedlot",
  min_prevalence=0.05,
  median_comparison_abundance = T, # default
  median_comparison_prevalence = FALSE, #default  
  min_abundance = 0, ##Imput_data has already been filtered
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_cDNA_pathways$fit_data_abundance$results) #Abundance results from MaAslin - 330


#Calculate confidence intervals, add taxonomy
maaslin_cDNA_pathways_2 <- maaslin_cDNA_pathways$fit_data_abundance$results %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) 
#left_join(input_taxonomy_cDNA, by = "feature")

##Final edits to put together for plot
maaslin_cDNA_pathways_3 <- maaslin_cDNA_pathways_2%>%
  #filter(qval_individual < 0.05) %>% ##Only significant taxa.qpval_individual is the corrected q-value of the individual association (only abundance. If it were qval_joint it'd be joint prevalence and abundance association)
  rename(group = name,
         pval = pval_individual,
         qval = qval_individual,
         pathway = feature)%>% ##Renaming
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
  select(pathway, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancombc_output_cDNA.pathways_pivot_3)
nrow(maaslin_cDNA_pathways_3) ##41 DA pathways between DNA and cDNA by MaAslin

##ANCOM and MaAslin together
DA_cDNA_pathways_MaAslinANCOM.data <- rbind(ancombc_output_cDNA.pathways_pivot_3, maaslin_cDNA_pathways_3) %>%
  filter(pathway %in% intersect(maaslin_cDNA_pathways_3$pathway,
                                  ancombc_output_cDNA.pathways_pivot_3$pathway)) ##Only going to plot those taxa DA by both tests

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_cDNA.pathways <- maaslin_cDNA_pathways$transformed_data %>% #transformed data is ra transformed otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "pathway", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  #left_join(input_taxonomy_cDNA, by = "feature")%>%
  filter(pathway %in% intersect(maaslin_cDNA_pathways_3$pathway,
                                  ancombc_output_cDNA.pathways_pivot_3$pathway))%>%
  left_join(metadata_maaslin_cDNA%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(pathway, Sample, logvalue, sample_type, plot)

#Bias-corrected abundances (ANCOM)
cDNA.pathways_log_corr_abn <- ancombc_output_cDNA.pathways$bias_correct_log_table %>%
  data.frame(check.names = F)%>% ##make into data frame
  rownames_to_column("pathway")%>%
  filter(pathway %in% intersect(maaslin_cDNA_pathways_3$pathway,
                                  ancombc_output_cDNA.pathways_pivot_3$pathway))%>% #keep only those pathways in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -pathway, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         sample_type = ifelse(grepl("W", Sample), "Water", "Feces"))

#put together objects to plot DA
DA_cDNA_pathways_together <- bind_rows(DA_cDNA_pathways_MaAslinANCOM.data, cDNA.pathways_log_corr_abn, RA_MaaslinAncom_cDNA.pathways) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_cDNA_pathways_together$plot <- factor(DA_cDNA_pathways_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_cDNA_pathways_together$sample_type <- factor(DA_cDNA_pathways_together$sample_type, levels = c("Feces", "Water"))

####Which are the biggest positive fold changes?#######
top_positive_fold_changes_cDNA <- DA_cDNA_pathways_together %>%
  filter(!is.na(test)) %>%
  # get the max coef per pathway for ordering
  group_by(pathway) %>%
  mutate(max_coef = max(coef, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_coef), pathway, test) %>%
  head(n = 18) #top 3 
top_positive_fold_changes_cDNA 

####Which are the biggest negative fold changes?#######
top_negative_fold_changes_cDNA <- DA_cDNA_pathways_together %>%
  filter(!is.na(test)) %>%
  # get the max coef per pathway for ordering
  group_by(pathway) %>%
  mutate(max_coef = max(coef, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_coef), pathway, test) %>%
  tail(n = 8) #top 4
top_negative_fold_changes_cDNA 


###PLOTTING DA#
#Plotting 
DA_pathways_cDNA_plot_MaAslinANCOM <-
  ggplot(data = DA_cDNA_pathways_together) +
  facet_wrap(~ plot, scales='free_x',
             nrow = 1,
             strip.position = "top") +
  geom_boxplot(data=DA_cDNA_pathways_together%>%filter(grepl("abundances", plot)),
               aes(x=pathway, y=logvalue, 
                   fill = sample_type, color = sample_type),
               notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(data=DA_cDNA_pathways_together%>%filter(grepl("abundances", plot)),
             aes(x=pathway, y=logvalue, fill = sample_type, color = sample_type),
             size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  # facet_nested(. ~ plot, scales='free_x',
  #              space='free_y',
  #              switch='y',
  #              strip=strip_nested(text_y=list(element_text(angle=0))),
  #              labeller=labeller(group=label_wrap_gen(width=10),
  #                                sub_group=label_wrap_gen(width=10))) +
  # ggplot(data=DA_cDNA_pathways_together%>%filter(grepl("abundances", plot)),
  #        aes(x=pathway, y=logvalue, fill = sample_type, color = sample_type)) +
  # geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  # geom_point(size = 1, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order = 1,title="Sample Type", title.position="top"),
         fill=guide_legend(order = 1,title="Sample Type", title.position="top"))+
  scale_fill_manual(values=sample.type.palette, 
                    labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  scale_color_manual(values = sample.type.palette, 
                     labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  new_scale_color()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_cDNA_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=pathway, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_cDNA_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=pathway, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_cDNA_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed'
             , alpha=0.5) +
  scale_x_discrete(position='bottom', labels = function(x) str_wrap(x, width = 40))+
  scale_y_continuous(position='right') +  # Default to break_labels otherwise 
  coord_flip() +
  labs(title = "Metatranscriptomic libraries (RNA (cDNA))" , x = "Functional Pathway")+
  scale_color_manual(values=c("red", "blue")) +
  scale_linetype_manual(values=c("11", "solid")) +
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Sample Type", title.position="top"),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 1.5)),
         linetype = guide_legend(title = "Fold change source", title.position = "top",
                                 override.aes = list(linewidth = 1),
                                 theme = theme(legend.key.width = unit(1.5, "cm"))),
         pch=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  theme_bw()+
  theme(legend.position="top", 
        legend.key=element_blank(),
        legend.title=element_text(size=15), legend.text=element_text(size=14),
        plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y=element_text(size= 20, angle=0, vjust= 1.045, face = "bold"), 
        axis.text.y=element_text(size=18, vjust = 0.5),
        strip.text=element_text(size=13, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_pathways_cDNA_plot_MaAslinANCOM


####FIGURE 5EandF#####
#Put together these 2 plots:
DA_pathways_cDNA_plot_MaAslinANCOM
DA_pathways_DNA_plot_MaAslinANCOM

figure5EandF <-plot_grid(DA_pathways_DNA_plot_MaAslinANCOM+
                           theme(plot.title = element_blank(),
                                 legend.position = "none",
                                 axis.title.y = element_blank(),
                                 axis.text.x=element_text(size=12)
                                 # strip.text = element_blank(),
                                 # strip.background = element_rect(fill = "white")
                                 ), 
                         DA_pathways_cDNA_plot_MaAslinANCOM+
                           theme(plot.title = element_blank(),
                                 legend.position = "none",
                                 axis.title.y = element_blank(),
                                 axis.text.x=element_text(size=12)
                                 # strip.text = element_blank(),
                                 # strip.background = element_rect(fill = "white")
                                 ), 
                         align = "v",
                         labels = c("E", "F"),
                         label_size = 32,
                         ncol = 1,
                         rel_heights = c(0.75, 0.35)
                         )
  # labs(title = "FUNCTIONAL PATHWAYS")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
figure5EandF
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure5EF.png", 
       plot = figure5EandF, 
       device = "png",
       dpi = 600,
       width = 19, 
       height = 5.5,  
       bg = "white")

##FECES (DNA vs cDNA)#####
###ANCOMBC#######
##Getting counts in fecal samples 
pathways.feces <- subset_samples(pathways_ps_filt, sample_type == "Feces") 
pathways.feces <- prune_taxa(taxa_sums(pathways.feces) > 0, pathways.feces) 
pathways.feces ## 655 taxa and 96 samples
ancom_pathways_feces.counts <- pathways.feces
ancom_pathways_feces.counts@sam_data$gen_material <- factor(ancom_pathways_feces.counts@sam_data$gen_material, levels = c("DNA", "cDNA"))##reorder gen_material as factor, DNA as "reference"
ancom_pathways_feces.counts@sam_data$feedlot <- factor(ancom_pathways_feces.counts@sam_data$feedlot, levels = c("1", "2", "3" ,"4","5")) ##make feedlot a factor (not num)

#Relative abundances for feces samples 
pathways.feces.tss <- subset_samples(pathways.tss_ps, sample_type == "Feces") 
pathways.feces.tss <- prune_taxa(taxa_sums(pathways.feces.tss) > 0, pathways.feces.tss) 
pathways.feces.tss #655 taxa and 96 samples

##Filtering out the low relative abundance (less than 0.5 %) pathways (RA object sums up to 1, so 0.005 is 0.5%)
pathways.feces.tss_filt <- filter_taxa(pathways.feces.tss, function(x) mean(x) > 0.005, TRUE) 
pathways.feces.tss_filt ## 70 taxa with mean RA > 0.5% across 96 samples (fecal samples) 
##Filtering those pathways (> 0.5% RA) on the "raw" counts phyloseq object for feces
ancom_pathways_feces.counts_filt <- prune_taxa(taxa_names(pathways.feces.tss_filt), ancom_pathways_feces.counts)
ancom_pathways_feces.counts_filt ##70 taxa, 96 samples

##running ancombc on the variable of interest (gen_material)
ancombc_output_feces.pathways <-ancombc2(data= ancom_pathways_feces.counts_filt, 
                                   assay_name = "counts", 
                                   tax_level = NULL,
                                   fix_formula = "gen_material + feedlot",
                                   rand_formula =  " (1 | original_sample)",
                                   # fix_formula = "gen_material", 
                                   # rand_formula =  "(1 | feedlot) + (1 | original_sample)", 
                                   prv_cut = 0.05, 
                                   lib_cut = 0, 
                                   group= "gen_material", 
                                   struc_zero = TRUE, 
                                   neg_lb = TRUE,
                                   alpha = 0.05, #default significance
                                   n_cl = 1, verbose = TRUE)


## extract results from comparisons 
ancombc_output_feces.pathways_res <- ancombc_output_feces.pathways$res %>%
  select(-matches("feedlot")) #not doing feedlot comparisons

#Pivot longer the results
ancombc_output_feces.pathways_pivot <- ancombc_output_feces.pathways_res %>%
  mutate(across(starts_with("lfc_gen_material"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancombc_output_feces.pathways_pivot

##Getting rid of _rounded suffix using sub command
ancombc_output_feces.pathways_pivot$group ##want to get rid of "_rounded"
ancombc_output_feces.pathways_pivot$group<- sub("_rounded", "", ancombc_output_feces.pathways_pivot$group) 
ancombc_output_feces.pathways_pivot$group #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancombc_output_feces.pathways_pivot <- ancombc_output_feces.pathways_pivot %>%
  mutate(group= case_when(
    group == "gen_materialcDNA" ~ "DNA vs RNA(cDNA)",
    TRUE ~ group ##keeps original name for groups not specified (DNA)
  ))
ancombc_output_feces.pathways_pivot$group ##Now the group names are shorter and more manageable
ancombc_output_feces.pathways_pivot<- ancombc_output_feces.pathways_pivot %>%
  rename(pathways_family = taxon) ##This ancombc was done at the genus level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancombc_output_feces.pathways_pivot_2 <- ancombc_output_feces.pathways_pivot %>%
  mutate(lower.ci = lfc_gen_materialcDNA - 1.96*se_gen_materialcDNA,
         upper.ci = lfc_gen_materialcDNA + 1.96*se_gen_materialcDNA)

##Final fix - up to make compatible with plotting
ancombc_output_feces.pathways_pivot_3 <- ancombc_output_feces.pathways_pivot_2 %>%
  filter (passed_ss_gen_materialcDNA == 1)%>%##Only want those that passed sensitivity testing
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
  select(pathways_family, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
nrow(ancombc_output_feces.pathways_pivot_3) ##30 DA pathways_families between DNA and cDNA with ANCOM


###MaAsLin3#######
#Data (otu relative abundances) and metadata for MaAslin
pathways.feces.tss #70 taxa and 96 samples
pathways.feces.tss_filt ##70 taxa with mean RA > 0.5% across 96 samples (fecal samples) 

maaslin_pathways_feces.tss <- pathways.feces.tss_filt
maaslin_pathways_feces.tss@sam_data$gen_material <- factor(maaslin_pathways_feces.tss@sam_data$gen_material, levels = c("DNA", "cDNA"))##reorder gen_material as factor, DNA as "reference"
maaslin_pathways_feces.tss@sam_data$feedlot <- factor(maaslin_pathways_feces.tss@sam_data$feedlot, levels = c("1", "2", "3", "4", "5")) ##make feedlot a factor (not num)

##OTU RA counts for Maaslin as a data frame 
data_maaslin_feces_pathways_tss <- data.frame(otu_table(maaslin_pathways_feces.tss), check.names = F)

##Sample metadata
metadata_maaslin_feces <- data.frame(sample_data(maaslin_pathways_feces.tss), check.names = F) 


##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_feces_pathwaysfamilies <-  maaslin3(
  input_data = data_maaslin_feces_pathways_tss, 
  input_metadata = metadata_maaslin_feces, 
  output = "MaAsLin3_feces",
  fixed_effect = c("gen_material", "feedlot"),
  random_effects = c("original_sample"),
  # fixed_effect = c("gen_material"),
  # random_effects = c("feedlot", "original_sample"),
  min_prevalence=0.05,  
  median_comparison_abundance = T, # default
  median_comparison_prevalence = FALSE, #default 
  min_abundance = 0, ##Not filtering by abundance
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_feces_pathwaysfamilies$fit_data_abundance$results) #Abundance results from MaAslin - 350


#Calculate confidence intervals
maaslin_feces_pathwaysfamilies_2 <- maaslin_feces_pathwaysfamilies$fit_data_abundance$results %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) 
#left_join(input_taxonomy_feces, by = "feature")

##Final edits to put together for plot
maaslin_feces_pathwaysfamilies_3 <- maaslin_feces_pathwaysfamilies_2%>%
  rename(group = name,
         pval = pval_individual,
         qval = qval_individual,
         pathways_family = feature)%>% ##Renaming
  mutate(group= case_when(
    group == "gen_materialcDNA" ~ "DNA vs RNA(cDNA)",
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
  filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  select(pathways_family, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancombc_output_feces.pathways_pivot_3)
nrow(maaslin_feces_pathwaysfamilies_3) ##45 DA pathways_families between DNA and cDNA by MaAslin

##ANCOM and MaAslin together
DA_feces_pathways_MaAslinANCOM.data <- rbind(ancombc_output_feces.pathways_pivot_3, maaslin_feces_pathwaysfamilies_3) %>%
  filter(pathways_family %in% intersect(maaslin_feces_pathwaysfamilies_3$pathways_family,
                                  ancombc_output_feces.pathways_pivot_3$pathways_family)) ##Only going to plot those taxa DA by both tests

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_feces.pathways <- maaslin_feces_pathwaysfamilies$transformed_data %>% #transformed data is ra transformed otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "pathways_family", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  #left_join(input_taxonomy_feces, by = "feature")%>%
  filter(pathways_family %in% intersect(maaslin_feces_pathwaysfamilies_3$pathways_family,
                                  ancombc_output_feces.pathways_pivot_3$pathways_family))%>%
  left_join(metadata_maaslin_feces%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(pathways_family, Sample, logvalue, gen_material, plot)

#Bias-corrected abundances (ANCOM)
feces.pathways_log_corr_abn <- ancombc_output_feces.pathways$bias_correct_log_table %>%
  data.frame(check.names = F)%>% ##make into data frame
  rownames_to_column("pathways_family")%>%
  filter(pathways_family %in% intersect(maaslin_feces_pathwaysfamilies_3$pathways_family,
                                  ancombc_output_feces.pathways_pivot_3$pathways_family))%>% #keep only those pathways in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -pathways_family, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         gen_material = ifelse(grepl("c", Sample), "cDNA", "DNA"))


#put together objects to plot DA
DA_feces_pathways_together <- bind_rows(DA_feces_pathways_MaAslinANCOM.data, feces.pathways_log_corr_abn, RA_MaaslinAncom_feces.pathways) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_feces_pathways_together$plot <- factor(DA_feces_pathways_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_feces_pathways_together$gen_material <- factor(DA_feces_pathways_together$gen_material, levels = c("DNA", "cDNA"))

#Which are the biggest fold changes?
DA_feces_pathways_together %>%
  filter(!is.na(test)) %>%
  # get the max coef per pathway for ordering
  group_by(pathways_family) %>%
  mutate(max_coef = max(coef, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_coef), pathways_family, test)


####Which are the biggest positive fold changes?#######
top_positive_fold_changes_feces <- DA_feces_pathways_together %>%
  filter(!is.na(test)) %>%
  # get the max coef per pathway for ordering
  group_by(pathways_family) %>%
  mutate(max_coef = max(coef, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_coef), pathways_family, test) %>%
  head(n = 20) #top 10 
top_positive_fold_changes_feces

####Which are the biggest negative fold changes?#######
top_negative_fold_changes_feces<- DA_feces_pathways_together %>%
  filter(!is.na(test)) %>%
  # get the max coef per pathway for ordering
  group_by(pathways_family) %>%
  mutate(max_coef = max(coef, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_coef), pathways_family, test) %>%
  tail(n = 16) #top 8
top_negative_fold_changes_feces

#Plotting 
DA_fecespathways_plot_MaAslinANCOM <-
  ggplot(data = DA_feces_pathways_together) +
  facet_wrap(~ plot, scales='free_x',
             nrow = 1,
             strip.position = "top") +
  # facet_nested(. ~ plot, scales='free_x',
  #              space='free_y',
  #              switch='y',
  #              strip=strip_nested(text_y=list(element_text(angle=0))),
  #              labeller=labeller(group=label_wrap_gen(width=10),
  #                                sub_group=label_wrap_gen(width=10))) +
  # ggplot(data=DA_feces_pathways_together%>%filter(grepl("abundances", plot)),
  #        aes(x=pathways_family, y=logvalue, fill = gen_material, color = gen_material)) +
  # geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  # geom_point(size = 1, shape = 18, position = position_dodge(width = 0.75)) +
  geom_boxplot(data=DA_feces_pathways_together%>%filter(grepl("abundances", plot)),
               aes(x=pathways_family, y=logvalue, 
                   fill = gen_material, color = gen_material),
               notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(data=DA_feces_pathways_together%>%filter(grepl("abundances", plot)),
             aes(x=pathways_family, y=logvalue, fill = gen_material, color = gen_material),
             size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order = 1,title="Library Type", title.position="top"))+
  scale_color_manual(values = gen.material.palette,
                     labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
  scale_fill_manual(values=gen.material.palette,
                    labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
  new_scale_color()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_feces_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=pathways_family, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_feces_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=pathways_family, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_feces_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed'
             , alpha=0.5) +
  # scale_x_discrete(position='bottom',
  #                  labels = function(x) str_wrap(x, width = 40))+
  scale_x_discrete(position='bottom', labels = function(x) str_wrap(x, width = 40))+
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
         shape=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  theme_bw()+
  labs(title = "FECES", x = "Functional Pathway")+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=14),
        plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size= 20, face = "bold"), 
        axis.text.y=element_text(size=16, vjust = 0.5),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text=element_text(size=16, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_fecespathways_plot_MaAslinANCOM

##CATCH BASINS - DNA vs cDNA #####
###ANCOMBC#######
##Getting counts in water samples 
pathways.water <- subset_samples(pathways_ps_filt, sample_type == "Water") 
pathways.water <- prune_taxa(taxa_sums(pathways.water) > 0, pathways.water) 
pathways.water ## 736 taxa and 24 samples
ancom_pathways_water.counts <- pathways.water
ancom_pathways_water.counts@sam_data$gen_material <- factor(ancom_pathways_water.counts@sam_data$gen_material, levels = c("DNA", "cDNA"))##reorder gen_material as factor, DNA as "reference"
ancom_pathways_water.counts@sam_data$feedlot <- factor(ancom_pathways_water.counts@sam_data$feedlot, levels = c(1, 2, 3 ,4 ,5)) ##make feedlot a factor (not num)

#Relative abundances for water samples 
pathways.water.tss <- subset_samples(pathways.tss_ps , sample_type == "Water") 
pathways.water.tss <- prune_taxa(taxa_sums(pathways.water.tss) > 0, pathways.water.tss) 
pathways.water.tss #736 taxa and 24 samples

##Filtering out the low relative abundance (less than 0.5%) pathways (RA object sums up to 1, so 0.005 is 0.5%)
pathways.water.tss_filt <- filter_taxa(pathways.water.tss, function(x) mean(x) > 0.005, TRUE) 
pathways.water.tss_filt ## 58 taxa with mean RA > 0.5% across 24 samples (water samples) 
##Filtering those pathways (> 0.5% RA) on the "raw" counts phyloseq object for water
ancom_pathways_water.counts_filt <- prune_taxa(taxa_names(pathways.water.tss_filt), ancom_pathways_water.counts)
ancom_pathways_water.counts_filt #58 taxa, 24 samples


##running ancombc on the variable of interest (gen_material)
ancombc_output_water.pathways <-ancombc2(data= ancom_pathways_water.counts_filt, 
                                   assay_name = "counts", 
                                   tax_level = NULL,
                                   fix_formula = "gen_material + feedlot",
                                   rand_formula =  " (1 | original_sample)",
                                   # fix_formula = "gen_material", 
                                   # rand_formula =  "(1 | feedlot) + (1 | original_sample)", 
                                   prv_cut = 0.05, 
                                   lib_cut = 0, 
                                   group= "gen_material", 
                                   struc_zero = TRUE, 
                                   neg_lb = TRUE,
                                   alpha = 0.05, #default significance
                                   n_cl = 1, verbose = TRUE)

## extract results from comparisons 
ancombc_output_water.pathways_res <- ancombc_output_water.pathways$res %>%
  select(-matches("feedlot")) #not doing feedlot comparisons

#Pivot longer the results
ancombc_output_water.pathways_pivot <- ancombc_output_water.pathways_res %>%
  mutate(across(starts_with("lfc_gen_material"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancombc_output_water.pathways_pivot

##Getting rid of _rounded suffix using sub command
ancombc_output_water.pathways_pivot$group ##want to get rid of "_rounded"
ancombc_output_water.pathways_pivot$group<- sub("_rounded", "", ancombc_output_water.pathways_pivot$group) 
ancombc_output_water.pathways_pivot$group #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancombc_output_water.pathways_pivot <- ancombc_output_water.pathways_pivot %>%
  mutate(group= case_when(
    group == "gen_materialcDNA" ~ "DNA vs cDNA",
    TRUE ~ group ##keeps original name for groups not specified (DNA)
  ))
ancombc_output_water.pathways_pivot$group ##Now the group names are shorter and more manageable
ancombc_output_water.pathways_pivot<- ancombc_output_water.pathways_pivot %>%
  rename(pathways_family = taxon) ##This ancombc was done at the genus level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancombc_output_water.pathways_pivot_2 <- ancombc_output_water.pathways_pivot %>%
  mutate(lower.ci = lfc_gen_materialcDNA - 1.96*se_gen_materialcDNA,
         upper.ci = lfc_gen_materialcDNA + 1.96*se_gen_materialcDNA)

##Final fix - up to make compatible with plotting
ancombc_output_water.pathways_pivot_3 <- ancombc_output_water.pathways_pivot_2 %>%
  filter (passed_ss_gen_materialcDNA == 1)%>%##Only want those that passed sensitivity testing
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
  select(pathways_family, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
nrow(ancombc_output_water.pathways_pivot_3) ##0 DA pathways_families between DNA and cDNA with ANCOM


###MaAsLin3#######
#Data (otu relative abundances) and metadata for MaAslin
pathways.water.tss ## Have this RA object 726 taxa and 24 samples
pathways.water.tss_filt ## 58 taxa with mean RA > 0.5% across 24 samples (water samples) 

maaslin_pathways_water.tss <- pathways.water.tss_filt
maaslin_pathways_water.tss@sam_data$gen_material <- factor(maaslin_pathways_water.tss@sam_data$gen_material, levels = c("DNA", "cDNA"))##reorder gen_material as factor, DNA as "reference"
maaslin_pathways_water.tss@sam_data$feedlot <- factor(maaslin_pathways_water.tss@sam_data$feedlot, levels = c("1", "2", "3", "4", "5")) ##make feedlot a factor (not num)

##OTU RA counts for Maaslin as a data frame 
data_maaslin_water_pathways_tss <- data.frame(otu_table(maaslin_pathways_water.tss), check.names = F)

##Sample metadata
metadata_maaslin_water <- data.frame(sample_data(maaslin_pathways_water.tss), check.names = F) 


##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_water_pathwaysfamilies <-  maaslin3(
  input_data = data_maaslin_water_pathways_tss, 
  input_metadata = metadata_maaslin_water, 
  output = "MaAsLin3_water",  
  fixed_effect = c("gen_material", "feedlot"),
  random_effects = c("original_sample"),
  # fixed_effect = c("gen_material"),
  # random_effects = c("feedlot", "original_sample"),
  min_prevalence=0.05,  
  median_comparison_abundance = T, # default
  median_comparison_prevalence = FALSE, #default 
  min_abundance = 0, ##Not filtering by abundance
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_water_pathwaysfamilies$fit_data_abundance$results) #Abundance results from MaAslin - 290


#Calculate confidence intervals, add taxonomy
maaslin_water_pathwaysfamilies_2 <- maaslin_water_pathwaysfamilies$fit_data_abundance$results %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) 
#left_join(input_taxonomy_water, by = "feature")

##Final edits to put together for plot
maaslin_water_pathwaysfamilies_3 <- maaslin_water_pathwaysfamilies_2%>%
  rename(group = name,
         pval = pval_individual,
         qval = qval_individual,
         pathways_family = feature)%>% ##Renaming
  mutate(group= case_when(
    group == "gen_materialcDNA" ~ "DNA vs RNA(cDNA)",
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
  select(pathways_family, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancombc_output_water.pathways_pivot_3)
nrow(maaslin_water_pathwaysfamilies_3) ##6 DA pathways between DNA and cDNA by MaAslin

##ANCOM and MaAslin together - However, since ANCOM didnt find any, wont have any
DA_water_pathways_MaAslinANCOM.data <- rbind(ancombc_output_water.pathways_pivot_3, maaslin_water_pathwaysfamilies_3) %>%
  filter(pathways_family %in% intersect(maaslin_water_pathwaysfamilies_3$pathways_family,
                                  ancombc_output_water.pathways_pivot_3$pathways_family)) ##Only going to plot those taxa DA by both tests

#STOPPING HERE SINCE THERE ARE NO DA PATHWAYS
# ##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
# ##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
# RA_MaaslinAncom_water.pathways <- maaslin_water_pathwaysfamilies$transformed_data %>% #transformed data is ra transformed otu counts from MaAslin output
#   rownames_to_column(var = "Sample")%>%
#   pivot_longer(cols = -Sample, names_to = "pathways_family", values_to = "logvalue") %>%
#   mutate(plot = 'log2(Relative abundances)') %>%
#   #left_join(input_taxonomy_water, by = "feature")%>%
#   filter(pathways_family %in% intersect(maaslin_water_pathwaysfamilies_3$pathways_family,
#                                   ancombc_output_water.pathways_pivot_3$pathways_family))%>%
#   left_join(metadata_maaslin_water%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
#   select(pathways_family, Sample, logvalue, gen_material, plot)
# 
# #Bias-corrected abundances (ANCOM)
# ##Get sampling fractions
# water.pathways_samp_frac <- ancombc_output_water.pathways$samp_frac
# 
# # Replace NA with 0 in sampling fractions
# water.pathways_samp_frac[is.na(water.pathways_samp_frac)] <-  0
# # Add pesudo-count (1) to avoid taking the log of 0,  then get log10 of feature_table
# water.pathways_log_obs_abn <-  log(ancombc_output_water.pathways$feature_table + 1)
# # Adjust the log observed abundances (for sampling fraction - bias)
# water.pathways_log_corr_abn <- t(t(water.pathways_log_obs_abn) - water.pathways_samp_frac) %>%
#   data.frame(check.names = F)%>% ##make into data frame
#   rownames_to_column("pathways_family")%>%
#   filter(pathways_family %in% intersect(maaslin_water_pathwaysfamilies_3$pathways_family,
#                                   ancombc_output_water.pathways_pivot_3$pathways_family))%>% #keep only those pathways in both the ANCOM and MaAslin outputs
#   pivot_longer(cols = -pathways_family, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
#   mutate(plot = 'log(Bias-corrected abundances)',
#          gen_material = ifelse(grepl("c", Sample), "cDNA", "DNA"))
# 
# #put together objects to plot DA
# DA_water_pathways_together <- bind_rows(DA_water_pathways_MaAslinANCOM.data, 
#                                         water.pathways_log_corr_abn, 
#                                         RA_MaaslinAncom_water.pathways) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
# DA_water_pathways_together$plot <- factor(DA_water_pathways_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
# DA_water_pathways_together$gen_material <- factor(DA_water_pathways_together$gen_material, levels = c("DNA", "cDNA"))
# 
# #Which are the biggest fold changes?
# DA_water_pathways_together %>%
#   filter(!is.na(test)) %>%
#   # get the max coef per pathway for ordering
#   group_by(pathways_family) %>%
#   mutate(max_coef = max(coef, na.rm = TRUE)) %>%
#   ungroup() %>%
#   arrange(desc(max_coef), pathways_family, test)
# 
# #Plotting 
# DA_waterpathways_plot_MaAslinANCOM <-
#   ggplot(data = DA_water_pathways_together) +
#   facet_wrap(~ plot, scales='free_x',
#              nrow = 1,
#              strip.position = "top") +
#   # facet_nested(. ~ plot, scales='free_x',
#   #              space='free_y',
#   #              switch='y',
#   #              strip=strip_nested(text_y=list(element_text(angle=0))),
#   #              labeller=labeller(group=label_wrap_gen(width=10),
#   #                                sub_group=label_wrap_gen(width=10))) +
#   # ggplot(data=DA_water_pathways_together%>%filter(grepl("abundances", plot)),
#   #        aes(x=pathways_family, y=logvalue, fill = gen_material, color = gen_material)) +
#   # geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
#   # geom_point(size = 1, shape = 18, position = position_dodge(width = 0.75)) +
#   geom_boxplot(data=DA_water_pathways_together%>%filter(grepl("abundances", plot)),
#                aes(x=pathways_family, y=logvalue, 
#                    fill = gen_material, color = gen_material),
#                notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
#   geom_point(data=DA_water_pathways_together%>%filter(grepl("abundances", plot)),
#              aes(x=pathways_family, y=logvalue, fill = gen_material, color = gen_material),
#              size = 1.5, shape = 18, position = position_dodge(width = 0.75)) +
#   guides(color=guide_legend(order = 1,title="Library Type", title.position="top"))+
#   scale_color_manual(values = gen.material.palette,
#                      labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
#   scale_fill_manual(values=gen.material.palette,
#                     labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
#   new_scale_color()+
#   geom_errorbar(inherit.aes=FALSE,
#                 data=DA_water_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
#                 aes(x=pathways_family, ymin=lower.ci, ymax=upper.ci,
#                     color=direction, linetype=test),
#                 width=0, position=position_dodge(0.75), size=0.75) +
#   geom_point(inherit.aes=FALSE,
#              data=DA_water_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
#              aes(x=pathways_family, y=coef, color=direction, pch=test),
#              position=position_dodge(0.75), size=3) +
#   geom_hline(data=DA_water_pathways_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
#              aes(yintercept=0),
#              size=0.5, linetype='dashed'
#              , alpha=0.5) +
#   # scale_x_discrete(position='bottom',
#   #                  labels = function(x) str_wrap(x, width = 40))+
#   scale_x_discrete(position='bottom')+
#   scale_y_continuous(position='right') +  # Default to break_labels otherwise 
#   coord_flip() +
#   scale_color_manual(values=c("red", "blue")) +
#   scale_linetype_manual(values=c("11", "solid")) +
#   scale_shape_manual(values=c(16, 15)) +
#   guides(fill=guide_legend(order=1, title="Library Type", title.position="top"),
#          color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 2.5)),
#          linetype = guide_legend(title = "Fold change source", title.position = "top",
#                                  override.aes = list(linewidth = 1),
#                                  theme = theme(legend.key.width = unit(1.5, "cm"))),
#          shape=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
#   theme_bw()+
#   labs(title = "CATCH BASINS", x = "Functional Pathway")+
#   theme(legend.position="top", legend.key=element_blank(),
#         legend.title=element_text(size=20), 
#         legend.text=element_text(size=14),
#         plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
#         plot.title.position = "plot",
#         axis.title.x=element_blank(), 
#         axis.text.x=element_text(size=10),
#         axis.title.y=element_text(size= 20, face = "bold"), 
#         plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
#         axis.text.y=element_text(size=16, vjust = 0.5),
#         strip.text=element_text(size=16, color = "white", face = "bold"),
#         strip.background=element_rect(fill='black'
#                                       , color='white'),
#         strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
#         panel.grid.minor = element_blank())
# DA_waterpathways_plot_MaAslinANCOM

####FIGURE 8 D - DA EMPHASIS ON DNA vs cDNA#####
#Just feces had DA pathways
DA_fecespathways_plot_MaAslinANCOM

figure8D <-plot_grid(DA_fecespathways_plot_MaAslinANCOM+
                         theme(plot.title = element_blank(),
                               legend.position = "none",
                               axis.title.y = element_blank(),
                               # strip.text = element_blank(),
                               # strip.background = element_rect(fill = "white")
                               ), 
                       align = "v", 
                       labels = c("D", " "),
                       label_size = 32,
                       ncol = 1)
# labs(title = "FUNCTIONAL PATHWAYS")+
# theme(plot.title = element_text(size = 40, face = "bold"))
figure8D
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure8D.png", 
       plot = figure8D, 
       device = "png",
       dpi = 600,
       width = 19, 
       height = 16, 
       bg = "white")

# ##UPSET PLOT - PATHWAYS #####
# #####DNA (CB vs. Feces)#####
# ##making the dataset for upset of UpSetR (binary matrix, present or absent)
# upset.pathways_ps.DNA <- MicrobiotaProcess::get_upset(pathways_ps.DNA, factorNames="sample_type") 
# upset.pathways_ps.DNA
# upset.pathways_ps.DNA <- upset.pathways_ps.DNA%>%
#   rename("CB" = "Water")
# #Plot
# upset_plot_pathways_DNA <-upset(upset.pathways_ps.DNA,
#                                 sets.bar.color = c("brown",
#                                                    "#4C72B0"),
#                                 order.by = "freq", text.scale = c(3, 2.5, 2.5, 1.5, 3, 2),
#                                 point.size = 5, line.size = 2, mainbar.y.label= "OTU count",
#                                 sets.x.label = "OTU count", 
#                                 set_size.show = F)
# upset_plot_pathways_DNA
# 
# ##Saving 
# svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure4E.svg",
#         width=7, height=6)
# upset_plot_pathways_DNA
# dev.off()
# 
# 
# #####cDNA (CB vs. Feces)#####
# ##making the dataset for upset of UpSetR (binary matrix, present or absent)
# upset.pathways_ps.cDNA <- MicrobiotaProcess::get_upset(pathways_ps.cDNA, factorNames="sample_type") 
# upset.pathways_ps.cDNA
# upset.pathways_ps.cDNA <- upset.pathways_ps.cDNA%>%
#   rename("CB" = "Water")
# 
# 
# ##Plot
# upset_plot_pathways_cDNA <-upset(upset.pathways_ps.cDNA, 
#                         sets.bar.color = c("brown",
#                                            "#4C72B0"),
#                         order.by = "freq", text.scale = c(3, 2.5, 2.5, 1.5, 3, 2),
#                         point.size = 5, line.size = 2, mainbar.y.label= "OTU count",
#                         sets.x.label = "OTU count", 
#                         set_size.show = F)
# upset_plot_pathways_cDNA
# 
# ##Saving
# svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure4F.svg",
#         width=7, height=6)
# upset_plot_pathways_cDNA
# dev.off()
# 
# #####Feces (cDNA vs. DNA)#####
# ##making the dataset for upset of UpSetR (binary matrix, present or absent)
# upset.pathways_ps.feces <- MicrobiotaProcess::get_upset(pathways_ps.feces, factorNames="gen_material") 
# upset.pathways_ps.feces
# upset.pathways_ps.feces <- upset.pathways_ps.feces%>%
#   rename("RNA (cDNA)" = "cDNA")
# 
# ##Plot
# upset_plot_pathways_feces <-UpSetR::upset(upset.pathways_ps.feces, 
#                                           sets.bar.color = c("#CC79A7", "#009E73"), 
#                                           order.by = "freq", text.scale = c(3, 2.5, 2.5, 1.5, 3, 2),
#                                           point.size = 5, line.size = 2, mainbar.y.label= "OTU count",
#                                           sets.x.label = "OTU count", 
#                                           set_size.show = F) 
# upset_plot_pathways_feces
# 
# ##Saving figure
# svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure8E.svg",
#         width=7, height=6)
# upset_plot_pathways_feces
# dev.off()
# 
# #####CATCH BASINS (cDNA vs. DNA)#####
# ##making the dataset for upset of UpSetR (binary matrix, present or absent)
# upset.pathways_ps.water <- MicrobiotaProcess::get_upset(pathways_ps.water, factorNames="gen_material") 
# upset.pathways_ps.water
# upset.pathways_ps.water <- upset.pathways_ps.water%>%
#   rename("RNA (cDNA)" = "cDNA")
# 
# ##Plot
# upset_plot_pathways_water <-grid.grabExpr(UpSetR::upset(upset.pathways_ps.water, 
#                                                sets.bar.color = c("#CC79A7", "#009E73"), 
#                                                order.by = "freq", text.scale = c(3, 2.5, 2.5, 1.5, 3, 2),
#                                                point.size = 5, line.size = 2, mainbar.y.label= "OTU count",
#                                                sets.x.label = "OTU count", 
#                                                set_size.show = F))
# upset_plot_pathways_water
# 
# ##Saving figure
# svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure8F.svg",
#         width=7, height=6)
# upset_plot_pathways_water
# dev.off()