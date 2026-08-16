#LOAD R PACKAGES ######
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


#SOURCE FUNCTIONS#####
source("Functions/MergeLowAbun_group_ARG.R")
source("Functions/top_taxa_legend_updated.R")

#IMPORT ARG COUNTS MATRIX#######
ARGcounts <- read.csv(
  'Data/SNPconfirmed_AMR_analytic_matrix.csv')

##Changing col names (sample IDs), to match them as they are in the metadata file
ARGcounts_newcolnames <- sapply(str_split(colnames(ARGcounts), "_"), `[`, 1)#Splitting col names by "_", then extracting the first part of each split column name 
colnames(ARGcounts) <- ARGcounts_newcolnames #replacing col names for new ones
colnames(ARGcounts) ##good, now sample names are OK!

#ARG ANNOTATIONS - downloaded from AMRplusplus#######
tax.table.ARG <- read.csv(
  'Data/megares_annotations_v3.00.csv',
  check.names = F)
unique(tax.table.ARG$type) 
# "Drugs"          "Multi-compound" "Biocides"       "Metals"         "Multi-compount"

#Need to fix "Multi-compount" to be Multi-compound
tax.table.ARG <- tax.table.ARG %>%
  mutate(type = str_replace_all(type, "Multi-compount", "Multi-compound"))
unique(tax.table.ARG$type) #Ok, now fixed 

tax.table.ARG <- tax.table.ARG%>%
  column_to_rownames(var = "header")%>% #make "header' rownames so it matches with the OTU table rownames
  rename_with(~ str_to_title(.))%>% #want the annotations with the first letter capitalized
  as.matrix() ##make into matrix so it is compatible with tax_table function from phyloseq

#MAKE ARG MATRIX INTO OTU TABLE########
otu_table_ARG <- ARGcounts%>%
  column_to_rownames(var = "gene")

#PHYLOSEQ########
OTU_ARG <-phyloseq::otu_table(otu_table_ARG, taxa_are_rows = TRUE)
TAX_ARG <-phyloseq::tax_table(tax.table.ARG)
phyloseq_ARG <- phyloseq(OTU_ARG, TAX_ARG) 
phyloseq_ARG ##3651 taxa and 123 samples (120 actual samples and 2 controls)

##Importing metadata##
##Metadata####
metadata <- read.csv('Data/Metadata_Feedlot_CatchBasins.csv', 
                     check.names = F,
                     row.names = "sampleID")

metadata$SampleID<- rownames(metadata) #Making a SampleID column 
metadata$original_sample <- sub("c$", "", metadata$SampleID) #Making a column for the sample that both the metagenome and metatranscriptomes come from (match 2)

##Host-free reads####
hostrem <- read.csv('Data/HostRem_Reads_Feedlot_CatchBasins.csv')

#Change zymo- to zymo. as they are in metadata
hostrem[which(hostrem$SampleID == "Zymo-1a"), "SampleID"] <- "Zymo.1a"
hostrem[which(hostrem$SampleID == "Zymo-1b"), "SampleID"] <- "Zymo.1b"

#Merge with metadata
hostrem <- hostrem %>%
  select(SampleID, Hostrem_output_total_num_seqs)%>%
  left_join(metadata, by = "SampleID")%>%
  rename(HostFree_Reads=Hostrem_output_total_num_seqs)#Hostrem_output_total_num_seqs is the total number of host free reads

#Want rownames as SampleID so I can merge into phyloseq object
rownames(hostrem) <- hostrem$SampleID

##Making into phyloseq-compatible object
sampledata_phyloseq <- sample_data(hostrem) ##use phyloseq function sample_data() to make metadata into phyloseq sample data object

# Add metadata to the phyloseq object###
ARG_data <- merge_phyloseq(phyloseq_ARG, sampledata_phyloseq)
sort(sample_names(ARG_data)) #OK
ARG_data ##3651 taxa and 123 samples (120 samples plus 2 Zymos and 1 EB. Other EBs and NTCs did not make it through)


#COLOR PALETTES#####
feedlot_palette <- c("1" = "#fcca46", 
                     "2" = "#fe7f2d", 
                     "3" = "#233d4d", 
                     "4"= "#3b9ab2", 
                     "5"= "#e1b6ff")
#Sample type
sample.type.palette <- c("Water" = "#4C72B0",
                         "Feces" = "brown") 
#Library Type
gen.material.palette <- c("cDNA" = "#009E73",  
                          "DNA"  = "#CC79A7" )  

#Sequencing batches
batch_palette <- c("no" = "#d19bac", 
                   "yes" = "#6a9c55") #sequencing batches

#Salmonella
salmonella.palette <- c("positive"= "#fc8d62", "negative" = "#8da0cb")


#PREPROCESSING#######
# some QC checks of the alligned reads per samples
min(sample_sums(ARG_data)) # 0 (EB)
max(sample_sums(ARG_data)) # 1,001,152 (Zymo.1a ) 
mean(sample_sums(ARG_data)) #  88,169.28
median(sample_sums(ARG_data)) # 79,729
sort(sample_sums(ARG_data)) #OK

### pulling out samples from ZYMOs and EBs
data2_ARG <- subset_samples(ARG_data, sample_type=="Water" | sample_type=="Feces") ##Pulling only Water and Feces samples (leaving out Zymo controls and the EB)
data2_ARG <- prune_taxa(taxa_sums(data2_ARG) > 0, data2_ARG) 
data2_ARG ###2696 taxa and 120 samples

##QC checks again
min(sample_sums(data2_ARG)) # 11,554 (F2W03c)
max(sample_sums(data2_ARG)) #162,970 (F3F07)
mean(sample_sums(data2_ARG)) #73901.33
median(sample_sums(data2_ARG)) #79357.5
sort(sample_sums(data2_ARG)) #OK

#PERCENTAGE OF NONHOST READS ALIGNED TO MEGARES GENES########
#Nonhost reads
nonhost_reads <- as.numeric(as.character(data2_ARG@sam_data$HostFree_Reads))
names(nonhost_reads) <- sample_names(data2_ARG)

#Total aligned reads
total_aligned_reads <- as.numeric(sample_sums(data2_ARG))
names(total_aligned_reads) <- sample_names(data2_ARG)

#Total aligned reads/Nonhost reads
Percentage_reads_aligned <- data.frame(
  SampleID = names(nonhost_reads),
  Percentage_reads_aligned = (total_aligned_reads / nonhost_reads) * 100)

##Add metadata 
Percentage_reads_aligned_metadata <- Percentage_reads_aligned %>%
  left_join(metadata, by = "SampleID")%>%
  mutate(gen_material= factor(gen_material, levels = c("DNA", "cDNA")))
  #filter(!Percentage_reads_aligned < 0.01)#Would you filter based on percentage reads aligned?

#Stats
summary(Percentage_reads_aligned_metadata$Percentage_reads_aligned)
sort(Percentage_reads_aligned_metadata$Percentage_reads_aligned) ##NOT DROPPING ANY 

#Descriptive stats per group
Percentage_reads_aligned_metadata %>%
  group_by(sample_type, gen_material)%>%
  summarise(mean_percentage_aligned_reads = mean(Percentage_reads_aligned),
            sd_percentage_aligned_reads = sd(Percentage_reads_aligned))
# sample_type gen_material mean_percentage_aligned_reads sd_percentage_aligned_reads
# Feces       DNA                                 0.145                       0.0234
# Feces       cDNA                                0.108                       0.0293
# Water       DNA                                 0.0450                      0.0152
# Water       cDNA                                0.0392                      0.0177

##Plots ####
###cDNA vs DNA faceted by Water and CB- SIGNIFICANT ####
percent_alignedreads_ARG.cDNAvsDNA.WaF<- ggplot(Percentage_reads_aligned_metadata%>%arrange(SampleID), 
                                            aes(x = gen_material, 
                                                y= Percentage_reads_aligned, 
                                                color = gen_material, 
                                                fill = gen_material)) +
  theme_bw() +
  labs(y= "Percentage (%) of\nAligned Reads", 
       color = "Library Type", fill = "Library Type") +
  facet_grid(~sample_type, scales = "free",  
             labeller = as_labeller(c("Feces" = "FECES", "Water" = "CATCH BASINS"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
  scale_fill_manual(values = gen.material.palette, 
                    labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  scale_color_manual(values = gen.material.palette, 
                     labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
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
percent_alignedreads_ARG.cDNAvsDNA.WaF
#Stats
##Feces
feces.Percentage_reads_aligned_metadata <- subset(Percentage_reads_aligned_metadata, 
                                                  sample_type == "Feces")
wilcox_test(feces.Percentage_reads_aligned_metadata%>%arrange(SampleID), 
            Percentage_reads_aligned ~ gen_material, 
            paired = T) #s., p = 7.11e-15
##CB
water.Percentage_reads_aligned_metadata <- subset(Percentage_reads_aligned_metadata, 
                                                  sample_type == "Water")
wilcox_test(water.Percentage_reads_aligned_metadata%>%arrange(SampleID), 
            Percentage_reads_aligned ~ gen_material, 
            paired = T
            ) #n.s., p = 0.0342

###CB vs Feces faceted by cDNA and DNA - Significant #####
percent_alignedreads_ARG_WvF.cDNAandDNA<- ggplot(Percentage_reads_aligned_metadata, 
                                                 aes(x = sample_type, 
                                                     y= Percentage_reads_aligned, 
                                                     color = sample_type, 
                                                     fill = sample_type)) +
  theme_bw() +
  labs(y= "Percentage (%) of\nAligned Reads", color = "Sample Type", fill = "Sample Type") +
  facet_grid(~gen_material,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("DNA" = "DNA", 
                                      "cDNA" = "RNA (cDNA)"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
  scale_color_manual(values = sample.type.palette, 
                     labels = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  scale_fill_manual(values = sample.type.palette, 
                    labels = c("Feces" = "Feces", "Water" = "Catch Basins")) +
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
percent_alignedreads_ARG_WvF.cDNAandDNA

#Stats
##cDNA
cDNA.Percentage_reads_aligned_metadata <- subset(Percentage_reads_aligned_metadata, 
                                                 gen_material == "cDNA")
wilcox_test(cDNA.Percentage_reads_aligned_metadata, 
            Percentage_reads_aligned ~ sample_type) #s., p = 2.22 e-09
##DNA
DNA.Percentage_reads_aligned_metadata <- subset(Percentage_reads_aligned_metadata, 
                                                gen_material == "DNA")
wilcox_test(DNA.Percentage_reads_aligned_metadata, 
            Percentage_reads_aligned~sample_type) #s., p = 1.43e-12


###Feedlot - NOT SIG ####
percent_alignedreads_ARG_feedlot <- ggplot(Percentage_reads_aligned_metadata,
                                           aes(x = feedlot, y= Percentage_reads_aligned, 
                                               color = factor(feedlot), fill = factor(feedlot))) +
  theme_bw() +
  labs(title = "NON-HOST READ ALIGNMENT TO MEGARes", 
       y= "Percentage (%) of\nAligned Reads", color = "Feedlot", fill = "Feedlot" ) +
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
percent_alignedreads_ARG_feedlot
#Stats
pairwise.wilcox.test(Percentage_reads_aligned_metadata$Percentage_reads_aligned, 
                     Percentage_reads_aligned_metadata$feedlot, p.adjust.method = "BH") #NS

###Checking for batch effects - reseq vs first run?####
percent_alignedreads_ARG_batches <- ggplot(Percentage_reads_aligned_metadata, aes(x = re_sequenced, y= Percentage_reads_aligned, color = re_sequenced, fill = re_sequenced)) +
  theme_bw() +
  labs(title = "NON-HOST READ ALIGNMENT TO MEGARes", 
       y= "Percentage (%) of\nAligned Reads", color = "Re-Sequenced", fill = "Re-Sequenced") +
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
  geom_pwc(method = "wilcox_test", p.adjust.method = "BH",
           label = "Wilcoxon, p = {p.adj}",
           hide.ns = TRUE,
           step.increase = 0.08,
           label.size = 5,
           tip.length = 0.02)
percent_alignedreads_ARG_batches

#Checking Stats
wilcox_test(Percentage_reads_aligned_metadata, Percentage_reads_aligned~re_sequenced) #S (p = 0.0336)

###SUPPLEMENTARY FIGURE?- PERCENTAGE ALIGNED READS####
percent_alignedreads_ARG_WvF.cDNAandDNA
percent_alignedreads_ARG.cDNAvsDNA.WaF

sfigureX <- plot_grid(percent_alignedreads_ARG_WvF.cDNAandDNA + 
                         theme(axis.title.y = element_text(size = 19)),
                        percent_alignedreads_ARG.cDNAvsDNA.WaF+ 
                         theme(axis.title.y = element_text(size = 19)), 
                       align = "v",
                       ncol = 1,
                       labels = "AUTO",
                       label_size = 22)
  #labs(title = "NON-HOST READ ALIGNMENT TO MEGARes")+
  #theme(plot.title = element_text(size = 30, face = "bold"))
sfigureX



#TAX GLOMMING - NOT NORMALIZED COUNTS##### 
###GROUP###
data2_ARG.group <- tax_glom(data2_ARG, taxrank = "Group", NArm = F)
data2_ARG.group #869 groups and 120 samples
###TYPE
data2_ARG.type <- tax_glom(data2_ARG.group, taxrank = "Type", NArm = F) 
data2_ARG.type #4 types (120 samples)
###CLASS
data2_ARG.class <- tax_glom(data2_ARG.group, taxrank = "Class", NArm = F) # classes
data2_ARG.class #54 classes (120 samples)
###MECHANISM
data2_ARG.mechanism <- tax_glom(data2_ARG.group, taxrank = "Mechanism", NArm = F) 
data2_ARG.mechanism #181 mechanisms (120 samples)


#TOTAL ALIGNED READS#######
sample.sums.arg <- sample_sums(data2_ARG.group) #making a sample sums object
sample.sums.df <- data.frame(SampleID = names(sample.sums.arg),
                             sample.sums.arg = sample.sums.arg) ##Turning into dataframe
data2_ARG.sampledata.df <- data.frame(phyloseq::sample_data(data2_ARG.group))

#combining sample sums with metadata
data2_ARG.df <- merge(data2_ARG.sampledata.df,
                      sample.sums.df,
                      by = "SampleID") %>%
  mutate(sample_type = factor(sample_type, levels = c("Feces", "Water")),
         gen_material = factor(gen_material, levels = c("DNA", "cDNA")))

##cDNA vs DNA faceted by CB and Feces- S. for feces####
aligned_reads_ARG.cDNAvsDNA.WaF<- ggplot(data2_ARG.df, 
                                            aes(x = gen_material, y= sample.sums.arg, color = gen_material, fill = gen_material)) +
  theme_bw() +
  labs(y= "Aligned Reads per Sample", color = "Library Type", fill = "Library Type") +
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
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))  +
  geom_pwc(method = "wilcox_test",
           label = "p = {p}",
           hide.ns = TRUE,
           method.args = list(paired = T),
           step.increase = 0.08,
           label.size = 6,
           tip.length = 0.02)
aligned_reads_ARG.cDNAvsDNA.WaF

###Stats
#CB
water.data2_ARG <- subset(data2_ARG.df, sample_type == "Water")
wilcox_test(water.data2_ARG%>%arrange(SampleID), 
            sample.sums.arg~gen_material, 
            paired = T) #s. p =0.0122
#Feces
feces.data2_ARG <- subset(data2_ARG.df, sample_type == "Feces")
wilcox_test(feces.data2_ARG%>%arrange(SampleID),
            sample.sums.arg~gen_material, 
            paired = T) #s. p = 1.55e-10 OK


##CB vs Feces faceted by cDNA and DNA - Significant for both#####
aligned_reads_ARG_WvF.cDNAandDNA<- ggplot(data2_ARG.df, aes(x = sample_type, y= sample.sums.arg, color = sample_type, fill = sample_type)) +
  theme_bw() +
  labs(y= "Aligned Reads per Sample", color = "Sample Type", fill = "Sample Type") +
  facet_grid(~gen_material,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("DNA" = "DNA", 
                                      "cDNA" = "RNA (cDNA)"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
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
aligned_reads_ARG_WvF.cDNAandDNA

#Stats
#cDNA
cDNA.data2_ARG.df <- subset(data2_ARG.df, gen_material == "cDNA")
wilcox_test(cDNA.data2_ARG.df, sample.sums.arg ~ sample_type) #s., p = 3.44 e -08
#DNA
DNA.data2_ARG.df <- subset(data2_ARG.df, gen_material == "DNA")
wilcox_test(DNA.data2_ARG.df, sample.sums.arg~sample_type) #s., p = 2.86e-12

##Feedlot - NS####
aligned_reads_ARG_feedlot <- ggplot(data2_ARG.df, aes(x = feedlot, y= sample.sums.arg, color = factor(feedlot), fill = factor(feedlot))) +
  theme_bw() +
  labs(title = "Aligned Reads per Sample", y= "Reads per Sample", color = "Feedlot", fill = "Feedlot" ) +
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
aligned_reads_ARG_feedlot
#Stats
pairwise.wilcox.test(data2_ARG.df$sample.sums.arg, data2_ARG.df$feedlot, p.adjust.method = "BH") #Ns

##Checking for batch effects - reseq vs first run ####
aligned_reads_ARG_batches <- ggplot(data2_ARG.df, aes(x = re_sequenced, y= sample.sums.arg, color = re_sequenced, fill = re_sequenced)) +
  theme_bw() +
  labs(title = "Aligned Reads per Sample", y= "Reads per Sample",color = "Re-Sequenced", fill = "Re-Sequenced") +
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
           label = "Wilcoxon, p = {p.adj}",
           hide.ns = TRUE,
           step.increase = 0.08,
           label.size = 5,
           tip.length = 0.02)
aligned_reads_ARG_batches

###Checking Stats
wilcox_test(data2_ARG.df, sample.sums.arg~re_sequenced) #S (p = 0.0168)

#TSS NORAMLIZATION (RA)####
any(sample_sums(data2_ARG.group)== 0) ## no samples with 0 OTUs

data2_ARG.tss <- transform_sample_counts(data2_ARG, function(x) x/sum(x)*100)

##TAX GLOMMING - TSS NORMALIZED COUNTS##### 
data2_ARG.tss.group <- transform_sample_counts(data2_ARG.group, function(x) x/sum(x)*100) ##Relative abundance 
data2_ARG.tss.group #869 groups

data2_ARG.tss.type <- tax_glom(data2_ARG.tss.group, taxrank = "Type", NArm = F) 
data2_ARG.tss.type #4 types (120 samples)

data2_ARG.tss.class <- tax_glom(data2_ARG.tss.group, taxrank = "Class", NArm = F) # classes
data2_ARG.tss.class #54 classes (120 samples)

data2_ARG.tss.mechanism <- tax_glom(data2_ARG.tss.group, taxrank = "Mechanism", NArm = F) 
data2_ARG.tss.mechanism #181 mechanisms (120 samples)


#ALPHA DIVERSITY - GENE GROUP LEVEL######
alpha_div1_ARG_group <- phyloseq::estimate_richness(data2_ARG.group, measures = c("Observed", 
                                                                                  "Shannon")) # richness, diversity
alpha_div2_ARG_group <- microbiome::evenness(data2_ARG.group, index = "pielou", 
                                       zeroes = FALSE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE. Keeps the focus on the taxa actually observed.
                                       detection = 0) ##evenness

# combine metrics with metadata
alpha_div_ARG_group <- cbind(alpha_div1_ARG_group, alpha_div2_ARG_group)
alpha_div_ARG_group
alpha_div_ARG_group_meta <- cbind(data2_ARG@sam_data, alpha_div_ARG_group) 
alpha_div_ARG_group_meta # metadata and div metrics
alpha_div_ARG_group_meta <- alpha_div_ARG_group_meta

#Factor variables
alpha_div_ARG_group_meta$feedlot <- factor(alpha_div_ARG_group_meta$feedlot,
                                 levels = c("1", "2", "3", "4", "5"))
alpha_div_ARG_group_meta$sample_type <- factor(alpha_div_ARG_group_meta$sample_type,
                                     levels = c("Feces", "Water"))
alpha_div_ARG_group_meta$gen_material <- factor(alpha_div_ARG_group_meta$gen_material,
                                         levels = c("DNA", "cDNA"))

#Pivot to long format 
alpha_div_ARG_group_long <- 
  alpha_div_ARG_group_meta %>%
  pivot_longer(cols = c(Observed, Shannon, pielou),  
               names_to = "alpha_div_metric", 
               values_to = "alpha_div_value") 
alpha_div_ARG_group_long


##Factoring alpha div metrics
alpha_div_ARG_group_long$alpha_div_metric<- factor(alpha_div_ARG_group_long$alpha_div_metric, levels = c("Observed","pielou", "Shannon"))

##LM MODEL ALPHA DIV INDECES#####
###Observed (Richness)#######
#How does its distribution look?
ggplot(alpha_div_ARG_group, aes (x = Observed))+
  geom_histogram()

#Model (LM)
# do the linear mixed model with random sampleID 
Observed_model_group <- lmerTest::lmer(Observed~ sample_type * gen_material * feedlot + (1 | original_sample),
                                       data = alpha_div_ARG_group_meta, 
                                       REML = TRUE)# for better estimate of the random-effects variance

#Check model
plot(Observed_model_group) 
qqnorm(residuals(Observed_model_group))
qqline(residuals(Observed_model_group))
shapiro.test(residuals(Observed_model_group)) #Not really normally distributed 

#Summary
summary(Observed_model_group)
#Confidence Intervals
confint(Observed_model_group)

#Anova type3 - instead of testing each coefficient individually, Type III ANOVA tests the factor as a whole.
Anova(Observed_model_group, type = "III")
Anova(Observed_model_group, type = "III", test.statistic = "F") 
#Interaction between sample_type and gen_material
#Interaction between sample_type and feedlot
#Interaction between sample_type, gen_material, and feedlot

#####Emmeans
##Interaction between sample_type and gen_material
emmeans(Observed_model_group, pairwise~sample_type|gen_material) ###No significant effect of sample type fr either DNA or RNA
emmeans(Observed_model_group, pairwise~gen_material|sample_type) #Feces DNA higher richness than feces cDNA 

##Interaction between sample_type and feedlot
emmeans(Observed_model_group, pairwise~feedlot|sample_type) 
#Differences between feedlots: Feedlot 1 vs feedlot 2 in fecal samples. For CB, feedlot 1 vs feedlot 4 and feedlot 2 vs feedlot 4.
emmeans(Observed_model_group, pairwise~sample_type|feedlot)
#Only feedlot = 4 had fecal richness higher than holding CB. Feedlot = 1 had lower richness in feces than CB

##Sample_type:gen_material:feedlot interaction.
emmeans(Observed_model_group, pairwise~feedlot|sample_type + gen_material) 
#For feces DNA and cDNA, feedlot 1 vs feedlot 2 (contrast estimate -49.33, -42.67)
#For CB DNA and cDNA, feedlot 2 vs feedlot 4 (contrast estimate 120.67, 115.33)
#For CB cDNA, feedlot 1 vs feedlot 4 (contrast estimate 141.00)
#For CB DNA, feedlot 2 vs feedlot 3 (contrast estimate 97.00)

emmeans(Observed_model_group, pairwise~sample_type|feedlot + gen_material)
#Within DNA and cDNA samples, in feedlot 4 (contrast estimate 70.3, 91.2) feces had higher richness
#In cDNA samples, in feedlot 1 feces had lower richness (contrast estimate -80.4)

emmeans(Observed_model_group, pairwise~gen_material|feedlot + sample_type)
#Within Fecal samples, in feedlot 2 (contrast estimate 18.2) and feedlot 4 (20.1), DNA samples had higher richness
#Within water samples, in feedlot 1 (contrast estimate -38.3) DNA samples had lower richness
#Within water samples, in feedlots 2 (contrast estimate 46.3) and 4 (41.0) DNA samples had higher richness

###Evenness#######
#How does its distribution look?
ggplot(alpha_div_ARG_group, aes (x = pielou))+
  geom_histogram()

#Model (LM)
Evenness_model_group <- lmerTest::lmer(pielou~ sample_type * gen_material * feedlot + (1 | original_sample),
                                       data = alpha_div_ARG_group_meta, 
                                       REML = TRUE)# for better estimate of the random-effects variance
#Check model
#plot(Evenness_model_group) 
qqnorm(residuals(Evenness_model_group))
qqline(residuals(Evenness_model_group))
shapiro.test(residuals(Evenness_model_group)) #Not really normally distributed 

#Summary
summary(Evenness_model_group)
#Confidence Intervals
confint(Evenness_model_group)
#####Anova type3
Anova(Evenness_model_group, type = "III")#Interaction between sample_type and gen_material, 
#Interaction between gen_material and feedlot
#Three way interaction as well
#Anova(Evenness_model_group, type = "III", test.statistic = "F") 

#####Emmeans
emmeans(Evenness_model_group, pairwise~gen_material|sample_type+feedlot) ##Sample type : gen_material: feedlot interaction 
emmeans(Evenness_model_group, pairwise~sample_type|gen_material+feedlot) ##Sample type : gen_material : feedlot interaction 
emmeans(Evenness_model_group, pairwise~feedlot|gen_material+sample_type) ##Sample type : gen_material : feedlot interaction 

emmeans(Evenness_model_group, pairwise~sample_type|gen_material) ##Sample type: gen_material interaction 
emmeans(Evenness_model_group, pairwise~gen_material|sample_type) ##Sample type: gen_material interaction 

emmeans(Evenness_model_group, pairwise~gen_material|feedlot) ##gen_material:feedlot interaction 
emmeans(Evenness_model_group, pairwise~feedlot|gen_material) ##gen_material:feedlot interaction 


###Shannon#######
#How does its distribution look?
ggplot(alpha_div_ARG_group, aes (x = Shannon))+
  geom_histogram()
alpha_div_ARG_group
#Model
Shannons_model_group <- lmerTest::lmer(Shannon ~ sample_type * gen_material * feedlot + (1 | original_sample),
                                       data = alpha_div_ARG_group_meta,
                                       REML = TRUE) # for better estimate of the random-effects variances

#Check model
plot(Shannons_model_group) 
qqnorm(residuals(Shannons_model_group))
qqline(residuals(Shannons_model_group))
shapiro.test(residuals(Shannons_model_group)) #Not really normally distributed 

#Summary
summary(Shannons_model_group) #Interaciton between sample_type and gen_material
#Confidence Intervals
confint(Shannons_model_group)

#Anova type3 - instead of testing each coefficient individually, Type III ANOVA tests the factor as a whole#
Anova(Shannons_model_group, type = "III")
##Interaction between sample_type and gen_material. 
##Interaction bwtween gen_material:feedlot
##Three way interaction sample_type:gen_material:feedlot

#Emmeans
emmeans(Shannons_model_group, pairwise~gen_material|sample_type+feedlot) ##Sample type : gen_material: feedlot interaction 
emmeans(Shannons_model_group, pairwise~sample_type|gen_material+feedlot) ##Sample type : gen_material : feedlot interaction 
emmeans(Shannons_model_group, pairwise~feedlot|gen_material+sample_type) ##Sample type : gen_material : feedlot interaction 

emmeans(Shannons_model_group, pairwise~sample_type|gen_material) ##Sample type: gen_material interaction 
emmeans(Shannons_model_group, pairwise~gen_material|sample_type) ##Sample type: gen_material interaction 

emmeans(Shannons_model_group, pairwise~gen_material|feedlot) ##gen_material:feedlot interaction 
emmeans(Shannons_model_group, pairwise~feedlot|gen_material) ##gen_material:feedlot interaction 

##ANOVA TABLES TOGETHER######
###RICHNESS#####
richness_anovaIII <- data.frame(
  Anova(Observed_model_group, type = "III", test.statistic = "F"),
  check.names = F)%>%
  mutate(Metric = "Richness (Observed)",
         Dataset = "Resistome",
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
  
  

###EVENNESS#####
evenness_anovaIII <- data.frame(
  Anova(Evenness_model_group, type = "III", test.statistic = "F"),
  check.names = F)%>%
  mutate(Metric = "Evenness (pielou)",
         Dataset = "Resistome",
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

###SHANNON#####
shannons_anovaIII <- data.frame(
  Anova(Shannons_model_group, type = "III", test.statistic = "F"),
  check.names = F)%>%
  mutate(Metric = "Shannon",
         Dataset = "Resistome",
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

#Bind dataframes
alpha_div_anovaIII_resistome <- bind_rows(richness_anovaIII,
          evenness_anovaIII,
          shannons_anovaIII)
alpha_div_anovaIII_resistome


####SUPPLEMENTARY TABLE 6_2######
stable6.2 <- alpha_div_anovaIII_resistome%>%
  select(Dataset, Metric, `Fixed Effects`, `F`, Df, `Pr(>F)`)

##EMMEANS TABLES#######
###RICHNESS####
####Pairwise~gen_material|sample_type+feedlot
richness_emmeans_gm_st_feedlot <- emmeans(Observed_model_group, 
                                          pairwise~gen_material|sample_type+feedlot) ##Sample type : gen_material: feedlot interaction 
richness_emmeans_gm_st_feedlot_df <- data.frame(richness_emmeans_gm_st_feedlot$contrasts, 
                                                check.names = F)%>%
  mutate(Metric = "Richness (Observed)",
         Dataset = "Resistome")

####Pairwise~sample_type|gen_material+feedlot
richness_emmeans_st_gm_feedlot <- emmeans(Observed_model_group, 
                                          pairwise~sample_type|gen_material+feedlot) ##Sample type : gen_material : feedlot interaction 
richness_emmeans_st_gm_feedlot_df <- data.frame(richness_emmeans_st_gm_feedlot$contrasts,
                                                check.names = F)%>%
  mutate(Metric = "Richness (Observed)",
         Dataset = "Resistome")

####Pairwise~feedlot|gen_material+sample_type
richness_emmeans_feedlot_gm_st <- emmeans(Observed_model_group, 
                                          pairwise~feedlot|gen_material+sample_type)
richness_emmeans_feedlot_gm_st_df <- data.frame(richness_emmeans_feedlot_gm_st$contrasts,
                                                check.names = F)%>%
  mutate(Metric = "Richness (Observed)",
         Dataset = "Resistome")


###EVENNESS####
####Pairwise~gen_material|sample_type+feedlot
evenness_emmeans_gm_st_feedlot <- emmeans(Evenness_model_group, 
                                       pairwise~gen_material|sample_type+feedlot) ##Sample type : gen_material: feedlot interaction 
evenness_emmeans_gm_st_feedlot_df <- data.frame(evenness_emmeans_gm_st_feedlot$contrasts, 
                                                check.names = F)%>%
  mutate(Metric = "Evenness (pielou)",
         Dataset = "Resistome")

####Pairwise~sample_type|gen_material+feedlot
evenness_emmeans_st_gm_feedlot <- emmeans(Evenness_model_group, 
                                       pairwise~sample_type|gen_material+feedlot) ##Sample type : gen_material : feedlot interaction 
evenness_emmeans_st_gm_feedlot_df <- data.frame(evenness_emmeans_st_gm_feedlot$contrasts,
                                             check.names = F)%>%
  mutate(Metric = "Evenness (pielou)",
         Dataset = "Resistome")
  
####Pairwise~feedlot|gen_material+sample_type
evenness_emmeans_feedlot_gm_st <- emmeans(Evenness_model_group, 
                                       pairwise~feedlot|gen_material+sample_type)
evenness_emmeans_feedlot_gm_st_df <- data.frame(evenness_emmeans_feedlot_gm_st$contrasts,
                                                check.names = F)%>%
  mutate(Metric = "Evenness (pielou)",
         Dataset = "Resistome")


###SHANNONS#####
####Pairwise~gen_material|sample_type+feedlot
shannon_emmeans_gm_st_feedlot <- emmeans(Shannons_model_group, 
                                          pairwise~gen_material|sample_type+feedlot) ##Sample type : gen_material: feedlot interaction 
shannon_emmeans_gm_st_feedlot_df <- data.frame(shannon_emmeans_gm_st_feedlot$contrasts, 
                                                check.names = F)%>%
  mutate(Metric = "Shannon",
         Dataset = "Resistome")

####Pairwise~sample_type|gen_material+feedlot
shannon_emmeans_st_gm_feedlot <- emmeans(Shannons_model_group, 
                                          pairwise~sample_type|gen_material+feedlot) ##Sample type : gen_material : feedlot interaction 
shannon_emmeans_st_gm_feedlot_df <- data.frame(shannon_emmeans_st_gm_feedlot$contrasts,
                                                check.names = F)%>%
  mutate(Metric = "Shannon",
         Dataset = "Resistome")

####Pairwise~feedlot|gen_material+sample_type
shannon_emmeans_feedlot_gm_st <- emmeans(Shannons_model_group, 
                                          pairwise~feedlot|gen_material+sample_type)
shannon_emmeans_feedlot_gm_st_df <- data.frame(shannon_emmeans_feedlot_gm_st$contrasts,
                                                check.names = F)%>%
  mutate(Metric = "Shannon",
         Dataset = "Resistome")

#Bind data frames for emmeans
alpha_div_emmeans_resistome <- bind_rows(
  richness_emmeans_st_gm_feedlot_df,
  richness_emmeans_gm_st_feedlot_df,
  richness_emmeans_feedlot_gm_st_df,
  evenness_emmeans_st_gm_feedlot_df,
  evenness_emmeans_gm_st_feedlot_df,
  evenness_emmeans_feedlot_gm_st_df,
  shannon_emmeans_st_gm_feedlot_df,
  shannon_emmeans_gm_st_feedlot_df,
  shannon_emmeans_feedlot_gm_st_df
  )
alpha_div_emmeans_resistome

#Editing
alpha_div_emmeans_resistome <- alpha_div_emmeans_resistome%>%
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
  
####SUPPLEMENTARY TABLE 7_2 ###### 
stable7.2 <- alpha_div_emmeans_resistome


##EMMEANS PLOTS######
###Plotting data - facet just by gen material ######
Shannon_group_WvF.cDNAandDNA_emmeans <- emmip(Shannons_model_group,~sample_type|gen_material,
      CIs = T, 
      type = "response",
      nesting.order = F, 
      plotit = F)%>%
  mutate(alpha_div_metric = "Shannon")

Evenness_group_WvF.cDNAandDNA_emmean <- emmip(Evenness_model_group,~sample_type|gen_material,
                                                    CIs = T, 
                                                    type = "response",
                                                    nesting.order = F, 
                                                    plotit = F)%>%
  mutate(alpha_div_metric = "pielou")
Observed_group_WvF.cDNAandDNA_emmean <- emmip(Observed_model_group,~sample_type|gen_material,
                                              CIs = T, 
                                              type = "response",
                                              nesting.order = F, 
                                              plotit = F)%>%
  mutate(alpha_div_metric = "Observed")

##Put them together 
alpha_div_emmeans_data <- bind_rows(Shannon_group_WvF.cDNAandDNA_emmeans, 
                                    Evenness_group_WvF.cDNAandDNA_emmean, 
                                    Observed_group_WvF.cDNAandDNA_emmean)

###Plotting data - facet by gen material and feedlot####
Shannon_group_WvF.cDNAandDNA_feedlot_emmeans <- emmip(Shannons_model_group,~sample_type|gen_material + feedlot,
                                                      CIs = T, 
                                                      type = "response",
                                                      nesting.order = F, 
                                                      plotit = F)%>%
  mutate(alpha_div_metric = "Shannon")

Evenness_group_WvF.cDNAandDNA_feedlot_emmean <- emmip(Evenness_model_group,~sample_type|gen_material + feedlot,
                                                      CIs = T, 
                                                      type = "response",
                                                      nesting.order = F, 
                                                      plotit = F)%>%
  mutate(alpha_div_metric = "pielou")
Observed_group_WvF.cDNAandDNA_feedlot_emmean <- emmip(Observed_model_group,~sample_type|gen_material + feedlot,
                                                      CIs = T, 
                                                      type = "response",
                                                      nesting.order = F, 
                                                      plotit = F)%>%
  mutate(alpha_div_metric = "Observed")

##Put them together 
alpha_div_emmeans_data_ARGgroup_WvsF_cDNAandDNA_feedlot <- bind_rows(Shannon_group_WvF.cDNAandDNA_feedlot_emmeans, 
                                                                     Evenness_group_WvF.cDNAandDNA_feedlot_emmean, 
                                                                     Observed_group_WvF.cDNAandDNA_feedlot_emmean)

####PLOTS#######
#####Feces vs CB facet by cDNA and DNA#####
alpha_div_emmeans_WvsF_DNAandcDNA_resistome <- alpha_div_emmeans_data %>%
  ggplot(aes(x = sample_type, y = yvar, color = sample_type)) +
  geom_jitter(aes(x = sample_type,
                  y = alpha_div_value,
                  shape = feedlot),
              alpha = 0.35,
              width = 0.2,
              size = 3,
              data = alpha_div_ARG_group_long) +#raw data
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
  labs(title= "RESISTOME") +
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
alpha_div_emmeans_WvsF_DNAandcDNA_resistome

#####Feces vs CB faceted by cDNA and DNA as well as feedlot#####
alpha_div_emmeans_WvsF_DNAandcDNA_feedlot_resistome <- alpha_div_emmeans_data_ARGgroup_WvsF_cDNAandDNA_feedlot %>%
  ggplot(aes(x = sample_type, y = yvar, color = sample_type)) +
  geom_jitter(aes(x = sample_type,
                  y = alpha_div_value),
              alpha = 0.35,
              width = 0.2,
              size = 3,
              data = alpha_div_ARG_group_long) +#raw data
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
  labs(title= "RESISTOME") +
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
alpha_div_emmeans_WvsF_DNAandcDNA_feedlot_resistome


######FIGURE 3B#####
figure3B <- plot_grid(alpha_div_emmeans_WvsF_DNAandcDNA_feedlot_resistome+
                        theme(plot.title = element_blank()), 
                      labels = c("B"), 
                      label_size = 22)
figure3B



#####DNA vs cDNA facetted by Feces and CB#####
alpha_div_emmeans_cDNAvsDNA_WandF_resistome <- alpha_div_emmeans_data %>%
  ggplot(aes(x = gen_material, y = yvar, color = gen_material)) +
  geom_jitter(aes(x = gen_material,
                  y = alpha_div_value,
                  shape = feedlot),
              alpha = 0.5,
              width = 0.2,
              size = 3,
              data = alpha_div_ARG_group_long) +#raw data
  geom_point(size = 4, shape = 20) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), 
                position = position_dodge(width = 0.5), 
                width = 0.2,
                linewidth = 1.5) +
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  scale_shape_manual(name = "Feedlot",
                     values = c(15,17,18,19, 12))+
  theme_bw() +
  labs(title= "RESISTOME", color = "Library Type") +
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
           label = "p.signif",
           step.increase = 0.1,
           size = 0.5,
           label.size = 5,
           tip.length = 0.02,
           hide.ns = T) 
alpha_div_emmeans_cDNAvsDNA_WandF_resistome

#####DNA vs cDNA facet by Feces and CB and feedlot#####
alpha_div_emmeans_cDNAvsDNA_WandF_feedlot_resistome <- 
  alpha_div_emmeans_data_ARGgroup_WvsF_cDNAandDNA_feedlot %>%
  ggplot(aes(x = gen_material, y = yvar, color = gen_material)) +
  geom_jitter(aes(x = gen_material,
                  y = alpha_div_value),
              alpha = 0.5,
              width = 0.2,
              size = 3,
              data = alpha_div_ARG_group_long) +#raw data
  geom_point(size = 4, shape = 20) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), 
                position = position_dodge(width = 0.5), 
                width = 0.2,
                linewidth = 1.5) +
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+
  theme_bw() +
  labs(title= "RESISTOME", color = "Library Type") +
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
alpha_div_emmeans_cDNAvsDNA_WandF_feedlot_resistome


######FIGURE 6B - FOCUS ON LIBRARY TYPE#####
figure6B <- plot_grid(alpha_div_emmeans_cDNAvsDNA_WandF_feedlot_resistome+
                        theme(plot.title = element_blank()), 
                      labels = c("B"), 
                      label_size = 22)



##BOX PLOTS ALPHA DIV ######
###Feces vs CB faceted by cDNA and DNA#####
alpha_div_ARG_group_WvF.cDNAandDNA <- ggplot(alpha_div_ARG_group_long, aes(x = sample_type,
                                                                    y= alpha_div_value, 
                                                                    fill= sample_type, colour = sample_type)) +
  theme_bw() +
  labs(
    title= "RESISTOME\nALPHA DIVERSITY", 
    color = "Sample Type", fill = "Sample Type") +
  facet_grid(alpha_div_metric~gen_material,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "DNA" = "DNA", 
                                      "cDNA" = "RNA (cDNA)",
                                      "1" = "1",
                                    "2" = "2", 
                                    "3" = "3",
                                    "4" = "4", 
                                    "5" = "5"))) +
  geom_jitter(alpha = 0.5,
              width = 0.2,
              shape = 18,
              size = 4) +#raw data
  geom_boxplot(alpha = 0.1) +
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
  # #For easy plotting of p values (I know it's not the same as ANOVA III, but keeps the significant results)
  stat_compare_means(method = "anova",
                     label.y.npc = "top",
                     hide.ns = TRUE,
                     show.legend = F,
                     label = "p.signif",
                     size = 8)
alpha_div_ARG_group_WvF.cDNAandDNA


###cDNA vs DNA faceted by CB and Feces#####
alpha_div_ARG_group_cDNAvsDNA.WandF  <- ggplot(alpha_div_ARG_group_long, aes(x = gen_material, y= alpha_div_value, fill= gen_material, colour = gen_material)) +
  theme_bw() +
  labs(title= "RESISTOME\nALPHA DIVERSITY", color = "Library Type", fill = "Library Type") +
  facet_grid(alpha_div_metric~sample_type,
             scales = "free",
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "Water" = "CATCH BASINS", 
                                      "Feces" = "FECES"))) +
  geom_jitter(alpha = 0.5,
              width = 0.2,
              shape = 18,
              size = 4) +#raw data
  geom_boxplot(alpha = 0.1) +
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
  stat_compare_means(method = "anova",
                     label.y.npc = "top",
                     hide.ns = TRUE,
                     show.legend = F,
                     label = "p.signif",
                     size = 8)
alpha_div_ARG_group_cDNAvsDNA.WandF 

###Reseq vs not reseq (batch effect check)#####
alpha_div_ARG_group_reseq <- ggplot(alpha_div_ARG_group_long, aes(x = re_sequenced, 
                                                                  y= alpha_div_value, 
                                                                  fill= re_sequenced, colour = re_sequenced)) +
  theme_bw() +
  labs(title= "ALPHA DIVERSITY") +
  facet_grid(alpha_div_metric ~ gen_material,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS", 
                                      "cDNA" = "RNA (cDNA)", 
                                      "DNA" = "DNA"))) +
  geom_jitter(alpha = 0.5,
              width = 0.2,
              shape = 18,
              size = 4) +#raw data
  geom_boxplot(alpha = 0.1) +
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
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
  # geom_pwc (method = "wilcox_test",
  #           label = "p = {p}",
  #           step.increase = 0.1,
  #           size = 0.5,
  #           label.size = 5,
  #           tip.length = 0.02,
  #           hide.ns = T)
alpha_div_ARG_group_reseq


# BETA DIVERSITY #####
##BRAY CURTIS####
###ONLY cDNA SAMPLES (METATRANSCRIPTOMIC)#############
##Subsetting only cDNA samples
data2_ARG.cDNA.tss <- subset_samples(data2_ARG.tss, gen_material=="cDNA")
data2_ARG.cDNA.tss <- prune_taxa(taxa_sums(data2_ARG.cDNA.tss) > 0, data2_ARG.cDNA.tss) 

##Distance matrix
data2_ARG.cDNA.tss.bray <- vegdist(t(data2_ARG.cDNA.tss@otu_table), method = "bray") 
data2_ARG.cDNA.tss.df <- as(data2_ARG.cDNA.tss@sam_data,"data.frame") %>% # make DF from metadata
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")))

#### ORDINATION
set.seed(87)
data2_ARG.cDNA.tss.bray.ord <- metaMDS(data2_ARG.cDNA.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

##Get ordination coordinates
data2_ARG.cDNA.tss.bray.plot <- ordiplot(data2_ARG.cDNA.tss.bray.ord$points)
#Extracts the scores (coordinates) of the points from the ordination plot object:
data2_ARG.cDNA.tss.bray.scrs <- scores(data2_ARG.cDNA.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata:
data2_ARG.cDNA.tss.bray.scrs <- cbind(as.data.frame(data2_ARG.cDNA.tss.bray.scrs), 
                                  sample_type = data2_ARG.cDNA.tss.df$sample_type, 
                                  feedlot = factor(data2_ARG.cDNA.tss.df$feedlot),
                                  SampleID = data2_ARG.cDNA.tss.df$SampleID) 
####PERMANOVA########
#Modelling sample type and feedlot 
set.seed(87)
cDNA_BC_adonis_sampletype_feedlot  <- adonis2(data2_ARG.cDNA.tss.bray ~ sample_type+feedlot,
                                       by = "margin",
                                       data2_ARG.cDNA.tss.df, permutations = 9999)
cDNA_BC_adonis_sampletype_feedlot #Significance of both sample type and feedlot

######SUPPLEMENTARY TABLE 5.10 #######
stable5.10 <- data.frame(cDNA_BC_adonis_sampletype_feedlot, check.names = F)%>%
  mutate(Dataset = "Resistome",
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
stable5.10

#PERMDISP - Sample type
# Run the betadisper function, average distance to centroid
bray.cDNA.sampletype.disp <- betadisper(data2_ARG.cDNA.tss.bray, 
                                        data2_ARG.cDNA.tss.df$sample_type)
bray.cDNA.sampletype.disp
##Then test by permuting
set.seed(87)
bray.cDNA.sampletype.permdisp <- permutest(bray.cDNA.sampletype.disp, permutations = 9999)
bray.cDNA.sampletype.permdisp ##Significant Feces vs Water p-value 1e-04

#PERMDISP - Feedlot
# Run the betadisper function, average distance to centroid
bray.cDNA.feedlot.disp <- betadisper(data2_ARG.cDNA.tss.bray, data2_ARG.cDNA.tss.df$feedlot)
bray.cDNA.feedlot.disp
##Then test by permuting
set.seed(87)
bray.cDNA.feedlot.permdisp <- permutest(bray.cDNA.feedlot.disp, permutations = 9999)
bray.cDNA.feedlot.permdisp ##Not significant for feedlots (p = 0.9)

####SAMPLE TYPE EFFECT#######
## BC
data2_ARG.cDNA.tss.bray.cent.sample_type <- aggregate(cbind(MDS1,MDS2) ~ sample_type, 
                                                      data = data2_ARG.cDNA.tss.bray.scrs,
                                                      FUN = mean) ##Centroids according to sample type (water and feces)
data2_ARG.cDNA.tss.bray.segs.sample_type <- merge(data2_ARG.cDNA.tss.bray.scrs, 
                                                  setNames(data2_ARG.cDNA.tss.bray.cent.sample_type, c("sample_type", "cMDS1","cMDS2")), by = 'sample_type', sort = F) ##add centroids to main scrs dataframe

##Abbreviated version of sample_type for easier plotting
data2_ARG.cDNA.tss.bray.segs.sample_type  <- data2_ARG.cDNA.tss.bray.segs.sample_type  %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"))

# Extract R2 and p-values
R2_cDNA_adonis_sample_type <- cDNA_BC_adonis_sampletype_feedlot$R2[1] 
pvalue_cDNA_adonis_sample_type<-  cDNA_BC_adonis_sampletype_feedlot$`Pr(>F)`[1]

#####PLOT #####
cDNA_BC_beta_div_spider_sampletype <- ggplot(data2_ARG.cDNA.tss.bray.segs.sample_type)+
                                               #filter(!SampleID %in% c("F2F02c", "F4W01c", "F4W02c"))) + 
  theme_bw() +
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
  annotate("text", x = -2.2, y = 0.55, ##change coordinates as needed
           label = "Sample type",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = -2.2, y = 0.55, ##change coordinates as needed
           label = paste0("R² = ", round(R2_cDNA_adonis_sample_type * 100, 1), "%",
                          "\np = ", round(pvalue_cDNA_adonis_sample_type, 4)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
cDNA_BC_beta_div_spider_sampletype

####FEEDLOT EFFECT#######
## BC
data2_ARG.cDNA.tss.bray.scrs #Already have the ordination coordinates with metadata
data2_ARG.cDNA.tss.bray.cent.feedlot <- aggregate(cbind(MDS1,MDS2) ~ feedlot, 
                                                  data = data2_ARG.cDNA.tss.bray.scrs,
                                                    #filter(!SampleID %in% c("F2F02","F2F02c", "F4W01c", "F4W02c")), 
                                              FUN = mean) ##Centroids according to feedlot
data2_ARG.cDNA.tss.bray.segs.feedlot <- merge(data2_ARG.cDNA.tss.bray.scrs, 
                                          setNames(data2_ARG.cDNA.tss.bray.cent.feedlot, c("feedlot", "cMDS1","cMDS2")), 
                                          by = 'feedlot', 
                                          sort = F) ##add centroids to main scrs dataframe

# Extract R2 and p-values
R2_cDNA_adonis_feedlot <- cDNA_BC_adonis_sampletype_feedlot$R2[2] 
pvalue_cDNA_adonis_feedlot<-  cDNA_BC_adonis_sampletype_feedlot$`Pr(>F)`[2]

#####PLOT#####
cDNA_BC_beta_div_spider_feedlot <- ggplot(data2_ARG.cDNA.tss.bray.segs.feedlot)+
                                            #filter(!SampleID %in% c("F2F02c", "F4W01c", "F4W02c"))) + 
  theme_bw() +
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
        legend.title = element_text(colour = "black", size = 22,face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
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
  annotate("text", x = -2.2, y = 0.6, ##change coordinates as needed
           label = "Feedlot",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x = -2.2, y = 0.6, ##change coordinates as needed
           label = paste0("R² = ", round(R2_cDNA_adonis_feedlot * 100, 1), "%",
                          "\np = ", round(pvalue_cDNA_adonis_feedlot, 3)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
cDNA_BC_beta_div_spider_feedlot

###ONLY DNA SAMPLES (METAGENOMICS)#############
##Subsetting only DNA samples
data2_ARG.DNA.tss <- subset_samples(data2_ARG.tss, gen_material=="DNA")
data2_ARG.DNA.tss <- prune_taxa(taxa_sums(data2_ARG.DNA.tss) > 0, data2_ARG.DNA.tss) 

##Distance matrix
data2_ARG.DNA.tss.bray <- vegdist(t(data2_ARG.DNA.tss@otu_table), method = "bray") 
data2_ARG.DNA.tss.df <- as(data2_ARG.DNA.tss@sam_data,"data.frame") %>% # make DF from metadata
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")))

#### ORDINATION
set.seed(87)
data2_ARG.DNA.tss.bray.ord <- metaMDS(data2_ARG.DNA.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

##Get ordination coordinates
data2_ARG.DNA.tss.bray.plot <- ordiplot(data2_ARG.DNA.tss.bray.ord$points)
#Extracts the scores (coordinates) of the points from the ordination plot object:
data2_ARG.DNA.tss.bray.scrs <- scores(data2_ARG.DNA.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata:
data2_ARG.DNA.tss.bray.scrs <- cbind(as.data.frame(data2_ARG.DNA.tss.bray.scrs), 
                                      sample_type = data2_ARG.DNA.tss.df$sample_type, 
                                      feedlot = factor(data2_ARG.DNA.tss.df$feedlot),
                                      SampleID = data2_ARG.DNA.tss.df$SampleID) 
####PERMANOVA########
#Modelling sample type and feedlot 
set.seed(87)
DNA_BC_adonis_sampletype_feedlot  <- adonis2(data2_ARG.DNA.tss.bray ~ sample_type+feedlot,
                                             by = "margin",
                                              data2_ARG.DNA.tss.df, permutations = 9999)
DNA_BC_adonis_sampletype_feedlot #Significant effects of sample type and feedlot


######SUPPLEMENTARY TABLE 5.9 #######
stable5.9 <- data.frame(DNA_BC_adonis_sampletype_feedlot, check.names = F)%>%
  mutate(Dataset = "Resistome",
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
stable5.9

#PERMDISP - Sample type
# Run the betadisper function, average distance to centroid
bray.DNA.sampletype.disp <- betadisper(data2_ARG.DNA.tss.bray, 
                                       data2_ARG.DNA.tss.df$sample_type)
bray.DNA.sampletype.disp
##Then test by permuting
set.seed(87)
bray.DNA.sampletype.permdisp <- permutest(bray.DNA.sampletype.disp, permutations = 9999)
bray.DNA.sampletype.permdisp 
##Feces vs Water p-value 1e-04

#PERMDISP - Feedlot
# Run the betadisper function, average distance to centroid
bray.DNA.feedlot.disp <- betadisper(data2_ARG.DNA.tss.bray, data2_ARG.DNA.tss.df$feedlot)
bray.DNA.feedlot.disp
##Then test by permuting
set.seed(87)
bray.DNA.feedlot.permdisp <- permutest(bray.DNA.feedlot.disp, permutations = 9999)
bray.DNA.feedlot.permdisp 
##Not significant for feedlots (p = 0.9)

####SAMPLE TYPE EFFECT#######
## BC
data2_ARG.DNA.tss.bray.cent.sample_type <- aggregate(cbind(MDS1,MDS2) ~ sample_type, 
                                                      data = data2_ARG.DNA.tss.bray.scrs,
                                                       #filter(!SampleID %in% c("F2F02", "F4W01", "F4W02")),
                                                      FUN = mean) ##Centroids according to sample type (water and feces)
data2_ARG.DNA.tss.bray.segs.sample_type <- merge(data2_ARG.DNA.tss.bray.scrs, 
                                                  setNames(data2_ARG.DNA.tss.bray.cent.sample_type, c("sample_type", "cMDS1","cMDS2")), by = 'sample_type', sort = F) ##add centroids to main scrs dataframe

##Abbreviated version of sample_type for easier plotting
data2_ARG.DNA.tss.bray.segs.sample_type  <- data2_ARG.DNA.tss.bray.segs.sample_type  %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F"))

# Extract R2 and p-values
R2_DNA_adonis_sample_type <- DNA_BC_adonis_sampletype_feedlot$R2[1] 
pvalue_DNA_adonis_sample_type<-  DNA_BC_adonis_sampletype_feedlot$`Pr(>F)`[1]

#####PLOT#####
DNA_BC_beta_div_spider_sampletype <- ggplot(data2_ARG.DNA.tss.bray.segs.sample_type)+
                                              #filter(!SampleID %in% c("F2F02c", "F4W01c", "F4W02c"))) + 
  theme_bw() +
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
  annotate("text", x= 0.1, y = 0.2, ##change coordinates as needed
           label = "Sample type",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text", x= 0.1, y = 0.2, ##change coordinates as needed
           label = paste0("R² = ", round(R2_DNA_adonis_sample_type * 100, 1), "%",
                          "\np = ", round(pvalue_DNA_adonis_sample_type, 4)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
DNA_BC_beta_div_spider_sampletype

####FEEDLOT EFFECT#######
## BC
data2_ARG.DNA.tss.bray.scrs #Already have the ordination coordinates with metadata
data2_ARG.DNA.tss.bray.cent.feedlot <- aggregate(cbind(MDS1,MDS2) ~ feedlot, 
                                                  data = data2_ARG.DNA.tss.bray.scrs,
                                                    #filter(!SampleID %in% c("F2F02","F2F02c", "F4W01c", "F4W02c")), 
                                                  FUN = mean) ##Centroids according to feedlot
data2_ARG.DNA.tss.bray.segs.feedlot <- merge(data2_ARG.DNA.tss.bray.scrs, 
                                              setNames(data2_ARG.DNA.tss.bray.cent.feedlot, c("feedlot", "cMDS1","cMDS2")), 
                                              by = 'feedlot', 
                                              sort = F) ##add centroids to main scrs dataframe

# Extract R2 and p-values
R2_DNA_adonis_feedlot <- DNA_BC_adonis_sampletype_feedlot$R2[2] 
pvalue_DNA_adonis_feedlot<-  DNA_BC_adonis_sampletype_feedlot$`Pr(>F)`[2]

#####PLOT - JUST feedlot #####
DNA_BC_beta_div_spider_feedlot <- ggplot(data2_ARG.DNA.tss.bray.segs.feedlot)+
                                            #filter(!SampleID %in% c("F2F02","F2F02c", "F4W01c", "F4W02c"))) +
  theme_bw() +
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
        legend.title = element_text(colour = "black", size = 22,face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
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
  annotate("text",  x = 1.7, y = -0.1, ##change coordinates as needed
           label = "Feedlot",
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable
  annotate("text",  x = 1.7, y = -0.1, ##change coordinates as needed
           label = paste0("R² = ", round(R2_DNA_adonis_feedlot * 100, 1), "%",
                          "\np = ", round(pvalue_DNA_adonis_feedlot, 3)),
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
DNA_BC_beta_div_spider_feedlot

###ONLY CATCH BASIN SAMPLES#######
##Subsetting only Water samples
data2_ARG.water.tss <- subset_samples(data2_ARG.tss, sample_type=="Water")
data2_ARG.water.tss <- prune_taxa(taxa_sums(data2_ARG.water.tss) > 0, data2_ARG.water.tss) 
data2_ARG.water.tss #2129 taxa and 24 samples

##Distance matrix
data2_ARG.water.tss.bray <- vegdist(t(data2_ARG.water.tss@otu_table), method = "bray") 
data2_ARG.water.tss.df <- as(data2_ARG.water.tss@sam_data,"data.frame")

##Adding a column concatenating "gen_material" and "feedlot", then making feedlot a factor
data2_ARG.water.tss.df<- data2_ARG.water.tss.df %>%
  mutate (gen_material_spec_2 = paste(gen_material, feedlot, sep = '_'),
          feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")))

####PERMANOVA#####
#Is there an interaction? 
set.seed(87)
water_BC_adonis_interaction <- adonis2(data2_ARG.water.tss.bray ~ gen_material * feedlot,
                                       by = "margin",
                                       data2_ARG.water.tss.df, 
                                       p.adjust.methods = "BH", permutations = 9999)
water_BC_adonis_interaction #No

##Model nucleic acid and feedlot, stratified by original sample from which DNA and RNA were extracted 
set.seed(87)
water_BC_adonis <- adonis2(data2_ARG.water.tss.bray ~ gen_material + feedlot, 
                           strata = data2_ARG.water.tss.df$original_sample,
                           by = "margin",
                           data2_ARG.water.tss.df, 
                           p.adjust.methods = "BH", 
                           permutations = 9999)
water_BC_adonis #11.6% of variation is due to gen_material (cDNA vs DNA) p= 0.0004883
#58.7% of variation is due to feedlot p= 0.0004883

######SUPPLEMENTARY TABLE 5.12 #######
stable5.12 <- data.frame(water_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Resistome",
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
stable5.12

#PERMDISP- Library type
# Run the betadisper function, average distance to centroid
bray.water.genmat.disp <- betadisper(data2_ARG.water.tss.bray, 
                                     data2_ARG.water.tss.df$gen_material)
bray.water.genmat.disp
##Then test by permuting
set.seed(87)
bray.water.genmat.permdisp <- permutest(bray.water.genmat.disp, permutations = 9999)
bray.water.genmat.permdisp #No difference in variance between DNA and cDNA (p = 0.57)

#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.water.feedlot.disp <- betadisper(data2_ARG.water.tss.bray, 
                                      data2_ARG.water.tss.df$feedlot)
bray.water.feedlot.disp
##Then test by permuting
set.seed(87)
bray.water.feedlot.permdisp <- permutest(bray.water.feedlot.disp, permutations = 9999, pairwise = 1)
bray.water.feedlot.permdisp #No difference in variance between feedlots (p = 0.2)

#### ORDINATION
set.seed(87)
data2_ARG.water.tss.bray.ord <- metaMDS(data2_ARG.water.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

#Extract ordination coordinates, add metadata
data2_ARG.water.tss.bray.plot <- ordiplot(data2_ARG.water.tss.bray.ord$points)
#Extract the scores (coordinates) of the points from the ordination plot object 
data2_ARG.water.tss.bray.scrs <- scores(data2_ARG.water.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
data2_ARG.water.tss.bray.scrs <- cbind(as.data.frame(data2_ARG.water.tss.bray.scrs), 
                                   gen_material = data2_ARG.water.tss.df$gen_material, 
                                   feedlot = factor(data2_ARG.water.tss.df$feedlot), 
                                   gen_material_spec_2 = data2_ARG.water.tss.df$gen_material_spec_2,
                                   sampleID = rownames(data2_ARG.water.tss.df))

####FEEDLOT EFFECT#####
## BC
data2_ARG.water.tss.bray.scrs #have coodrinates and metadata
##Centroids according to feedlot
data2_ARG.water.feedlot.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = data2_ARG.water.tss.bray.scrs, FUN = mean) 
data2_ARG.water.feedlot.tss.bray.segs <- merge(data2_ARG.water.tss.bray.scrs, 
                                           setNames(data2_ARG.water.feedlot.tss.bray.cent, c("feedlot","cMDS1","cMDS2")),
                                           by = 'feedlot', sort = F)

# Extract R2 and p-values
##Model gen_material and feedlot
water_BC_adonis #have the model
# Extract R2 and p-values
R2_water_BC_adonis_feedlot <- water_BC_adonis$R2[2] 
pvalue_water_BC_adonis_feedlot<-  water_BC_adonis$`Pr(>F)`[2]

#### PLOT
water_feedlot_BC_beta_div <- ggplot(data2_ARG.water.feedlot.tss.bray.segs) + theme_bw() +
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
  annotate("text", x = -0.9, y = 0.7, ##change coordinates as needed
           label = "Feedlot", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (feedlot)
  annotate("text", x = -0.9, y = 0.7, ##change coordinates as needed
           label = paste("R² = ", round(R2_water_BC_adonis_feedlot * 100, 1), "%",
                         "\np = ", round(pvalue_water_BC_adonis_feedlot, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black")+ # Annotate R² and p-values
  guides(color = guide_legend(override.aes = list(size = 7)))
water_feedlot_BC_beta_div 

####LIBRARY TYPE EFFECT#####
## BC
data2_ARG.water.tss.bray.scrs #have coordinates and metadata
##Centroids according to gen_material (cDNA/DNA)
data2_ARG.water.genmat.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ gen_material, 
                                              data = data2_ARG.water.tss.bray.scrs, FUN = mean) 
data2_ARG.water.genmat.tss.bray.segs <- merge(data2_ARG.water.tss.bray.scrs, 
                                          setNames(data2_ARG.water.genmat.tss.bray.cent, c("gen_material","cMDS1","cMDS2")), 
                                          by = 'gen_material', sort = F)%>%
  mutate(gen_material_plot = dplyr::recode(gen_material, 
                                           "DNA" = "DNA", 
                                           "cDNA" = "RNA"))

# Extract R2 and p-values for genmat
water_BC_adonis #have the model
# Extract R2 and p-values
R2_water_BC_adonis_gen_material <- water_BC_adonis$R2[1] 
pvalue_water_BC_adonis_gen_material<-  water_BC_adonis$`Pr(>F)`[1] 

#### PLOT
water_genmat_BC_beta_div <- ggplot(data2_ARG.water.genmat.tss.bray.segs) + theme_bw() +
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
  annotate("text", x = 0.5, y = 0.55, ##change coordinates as needed
           label = "Library Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Library Type)
  annotate("text", x = 0.5, y = 0.55, ##change coordinates as needed
           label = paste("R² = ", round(R2_water_BC_adonis_gen_material * 100, 1), "%",
                         "\np = ", round(pvalue_water_BC_adonis_gen_material, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
water_genmat_BC_beta_div 

###ONLY FECES#######
##Subsetting only feces samples
data2_ARG.feces.tss <- subset_samples(data2_ARG.tss, sample_type=="Feces")
data2_ARG.feces.tss <- prune_taxa(taxa_sums(data2_ARG.feces.tss) > 0, data2_ARG.feces.tss) 

##Distance matrix
data2_ARG.feces.tss.bray <- vegdist(t(data2_ARG.feces.tss@otu_table), method = "bray") 
data2_ARG.feces.tss.df <- as(data2_ARG.feces.tss@sam_data,"data.frame") # make DF from metadata

##Adding a column concatenating "gen_material" and "feedlot", then making feedlot a factor
data2_ARG.feces.tss.df<- data2_ARG.feces.tss.df %>%
  mutate (gen_material_spec_2 = paste(gen_material, feedlot, sep = '_'),
          feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")))

####PERMANOVA#####
#Is there an interaction? 
set.seed(87)
feces_BC_adonis_interaction <- adonis2(data2_ARG.feces.tss.bray ~ gen_material * feedlot, 
                                       by = "margin",data2_ARG.feces.tss.df, 
                                       p.adjust.methods = "BH", permutations = 9999)
feces_BC_adonis_interaction #No

##Model library type and feedlot, stratify by original sample from which DNA and RNA were extracted from 
set.seed(87)
feces_BC_adonis <- adonis2(data2_ARG.feces.tss.bray ~ gen_material + feedlot,
                           strata = data2_ARG.feces.tss.df$original_sample,
                           by = "margin",
                           data2_ARG.feces.tss.df, 
                           p.adjust.methods = "BH", permutations = 9999)
feces_BC_adonis 
#25.6% of variation is due to gen_material (cDNA vs DNA) p= 0.0001
#13.2% of variation is due to feedlot p= 0.0001

######SUPPLEMENTARY TABLE 5.11 #######
stable5.11 <- data.frame(feces_BC_adonis, check.names = F)%>%
  mutate(Dataset = "Resistome",
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
stable5.11

#PERMDISP- Feedlot
# Run the betadisper function, average distance to centroid
bray.feces.feedlot.disp <- betadisper(data2_ARG.feces.tss.bray, data2_ARG.feces.tss.df$feedlot)
bray.feces.feedlot.disp
##Then test by permuting
set.seed(87)
bray.feces.feedlot.permdisp <- permutest(bray.feces.feedlot.disp, permutations = 9999, pairwise = 1)
bray.feces.feedlot.permdisp 
##No Differentf variance (p=0.27)

#PERMDISP- Library table
# Run the betadisper function, average distance to centroid
bray.feces.genmat.disp <- betadisper(data2_ARG.feces.tss.bray, data2_ARG.feces.tss.df$gen_material)
bray.feces.genmat.disp
##Then test by permuting
set.seed(87)
bray.feces.genmat.permdisp <- permutest(bray.feces.genmat.disp, permutations = 9999)
bray.feces.genmat.permdisp 
#Difference in dispersions - variance between DNA and cDNA (p = 0.0284)

#### ORDINATION
set.seed(87)
data2_ARG.feces.tss.bray.ord <- metaMDS(data2_ARG.feces.tss.bray, k= 2, try= 20, trymax= 1000, autotransform = F)

#Extract ordination coordinates, add metadata
data2_ARG.feces.tss.bray.plot <- ordiplot(data2_ARG.feces.tss.bray.ord$points)
#Extract the scores (coordinates) of the points from the ordination plot object 
data2_ARG.feces.tss.bray.scrs <- scores(data2_ARG.feces.tss.bray.plot, display = "sites") 
#This step creates a new data frame that includes both the ordination scores (MDS) and additional metadata
data2_ARG.feces.tss.bray.scrs <- cbind(as.data.frame(data2_ARG.feces.tss.bray.scrs), 
                                   gen_material = data2_ARG.feces.tss.df$gen_material, 
                                   feedlot = factor(data2_ARG.feces.tss.df$feedlot), 
                                   gen_material_spec_2 = data2_ARG.feces.tss.df$gen_material_spec_2,
                                   sampleID = rownames(data2_ARG.feces.tss.df))

####FEEDLOT EFFECT#####
## BC
data2_ARG.feces.tss.bray.scrs #have coodrinates and metadata
##Centroids according to feedlot
data2_ARG.feces.feedlot.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ feedlot, data = data2_ARG.feces.tss.bray.scrs, FUN = mean) 
data2_ARG.feces.feedlot.tss.bray.segs <- merge(data2_ARG.feces.tss.bray.scrs, 
                                           setNames(data2_ARG.feces.feedlot.tss.bray.cent, c("feedlot","cMDS1","cMDS2")),
                                           by = 'feedlot', sort = F)

# Extract R2 and p-values
##Model gen_material and feedlot
feces_BC_adonis #have the model
R2_feces_BC_adonis_feedlot <- feces_BC_adonis$R2[2]
pvalue_feces_BC_adonis_feedlot<- feces_BC_adonis$`Pr(>F)`[2] #p value

#### PLOT
feces_feedlot_BC_beta_div <- ggplot(data2_ARG.feces.feedlot.tss.bray.segs) + theme_bw() +
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
  annotate("text", x = 0.8, y = -0.5, ##change coordinates as needed
           label = "Feedlot", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (feedlot)
  annotate("text", x = 0.8, y = -0.5, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces_BC_adonis_feedlot * 100, 1), "%",
                         "\np = ", round(pvalue_feces_BC_adonis_feedlot, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black")+ # Annotate R² and p-values
  guides(color = guide_legend(override.aes = list(size = 7)))
feces_feedlot_BC_beta_div 

####LIBRARY TYPE EFFECT#####
## BC
data2_ARG.feces.tss.bray.scrs #have coordinates and metadata
##Centroids according to gen_material (cDNA/DNA)
data2_ARG.feces.genmat.tss.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ gen_material, 
                                              data = data2_ARG.feces.tss.bray.scrs, FUN = mean) 
data2_ARG.feces.genmat.tss.bray.segs <- merge(data2_ARG.feces.tss.bray.scrs, 
                                          setNames(data2_ARG.feces.genmat.tss.bray.cent, c("gen_material","cMDS1","cMDS2")), 
                                          by = 'gen_material', sort = F)%>%
  mutate(gen_material_plot = dplyr::recode(gen_material, 
                                    "DNA" = "DNA", 
                                    "cDNA" = "RNA"))

# Extract R2 and p-values for genmat
feces_BC_adonis #have the model
R2_feces_BC_adonis_gen_material <- feces_BC_adonis$R2[1] #R2
pvalue_feces_BC_adonis_gen_material <- feces_BC_adonis$`Pr(>F)`[1]#p value

#### PLOT
feces_genmat_BC_beta_div <- ggplot(data2_ARG.feces.genmat.tss.bray.segs) + theme_bw() +
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
  geom_text(aes (x= cMDS1, y = cMDS2, label= gen_material_plot), colour= "white", size = 3, fontface = "bold") +
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
  annotate("text", x = 0.8, y = -0.35, ##change coordinates as needed
           label = "Library Type", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Library Type)
  annotate("text", x = 0.8, y = -0.35, ##change coordinates as needed
           label = paste("R² = ", round(R2_feces_BC_adonis_gen_material * 100, 1), "%",
                         "\np = ", round(pvalue_feces_BC_adonis_gen_material, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values
feces_genmat_BC_beta_div 

####SUPPLEMENTARY FIGURE 6CD#######
#####Effect of feedlot On DNA and cDNA#####
sfigure6CandD <- ggarrange(DNA_BC_beta_div_spider_feedlot+
                             theme(plot.title = element_blank()), 
                           cDNA_BC_beta_div_spider_feedlot+
                             theme(plot.title = element_blank()),
                           labels = c("C", "D"), 
                           font.label = list(size = 22),
                           legend = "none",
                           common.legend = T)
  # labs(title = "RESISTOME")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
sfigure6CandD


####SUPPLEMENTARY FIGURE 7CD#######
#####Effect of feedlot On Feces and CB#####
sfigure7CandD <- ggarrange(feces_feedlot_BC_beta_div+
                             theme(plot.title = element_blank()), 
                           water_feedlot_BC_beta_div+
                             theme(plot.title = element_blank()),
                           labels = c("C", "D"), 
                           font.label = list(size = 22),
                           legend = "none",
                           common.legend = T)
  # labs(title = "RESISTOME")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
sfigure7CandD


###SUPPLEMENTARY TABLE 5- RESISTOME SECTION#######
stable5_resistome <- bind_rows(stable5.9,
                               stable5.10,
                               stable5.11,
                               stable5.12)%>%  
  select(Dataset, `Library Type`, `Sample Type`, `Fixed Effect`, Df, SumOfSqs, R2, `F`, `Pr(>F)`)
stable5_resistome

#DENDROGRAMS #####
##All samples####
data2_ARG.tss.bray <- vegdist(t(data2_ARG.tss.group@otu_table), method = "bray") 
data2_ARG.bray.hclust <- hclust(data2_ARG.tss.bray, method = "ward.D2")
plot(data2_ARG.bray.hclust, hang = -1)
data2_ARG.bray.dendro <- as.dendrogram(data2_ARG.bray.hclust) # Build dendrogram object from hclust results
data2_ARG.bray.dendro.data <- dendro_data(data2_ARG.bray.dendro, type = "rectangle") # Extract the dendrogram plot data. Type wil draw rectangular lines
data2_ARG.tss.group@sam_data$sampleID <- rownames(data2_ARG.tss.group@sam_data) ##need a column with sample_ID
data2_ARG.bray.dendro.metadata <- 
  data.frame(data2_ARG.tss.group@sam_data) %>%
  mutate(sample_type.abbrv = dplyr::recode(sample_type, "Water"= "CB", "Feces"= "F")) %>%
  mutate(gen_material.abbrv = dplyr::recode(gen_material, "cDNA"= "T", "DNA"= "G"))
data2_ARG.bray.dendro.data$labels <- data2_ARG.bray.dendro.data$labels %>%
  left_join(data2_ARG.bray.dendro.metadata, by = c("label" = "sampleID")) 


##Plot
dendro.bray.ARG.plot <- ggplot(data2_ARG.bray.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  #Sample type
  geom_point(data = data2_ARG.bray.dendro.data$labels, 
             aes(x = x, y = y, colour = sample_type),size = 5, shape = 15, 
             position = position_nudge(y = -0.2))+
  scale_color_manual(name = "Sample Type", values = sample.type.palette,
                     label = c("Feces" = "Feces", "Water" = "Catch Basins")) +
  guides(color=guide_legend(title.position="top"))+
  new_scale_color()+ 
  #Feedlots
  geom_point(
    data = data2_ARG.bray.dendro.data$labels, 
    aes(x = x, y = y, colour = factor(feedlot)),size = 5, shape = 15, 
    position = position_nudge(y = -0.5))+
  scale_color_manual(name = "Feedlot", values = feedlot_palette) +
  guides(color=guide_legend(title.position="top"))+
  new_scale_color()+
  #Nucleic acid
  geom_point(data = data2_ARG.bray.dendro.data$labels,  
             aes(x = x, y = y, color = gen_material),
             size = 5, shape = 15, 
             position = position_nudge(y = -0.8)) +
  guides(color=guide_legend(title.position="top"))+
  scale_color_manual(name = "Library Type",values = gen.material.palette,  
                     label = c("DNA" = "DNA", "cDNA" = "RNA(cDNA)")) +
  #Sample type
  geom_text(data = data2_ARG.bray.dendro.data$labels, 
            aes(x=x, y=y, label = sample_type.abbrv), 
            colour = "white", size =3, 
            position = position_nudge(y=-0.2), 
            fontface = "bold")+
  #Feedlot
  geom_text(data = data2_ARG.bray.dendro.data$labels, 
            aes(x=x, y=y, label = factor(feedlot)), 
            colour = "white", size =4, 
            position = position_nudge(y=-0.5), fontface = "bold") +
  #Library Type
  geom_text(data = data2_ARG.bray.dendro.data$labels, 
            aes(x=x, y=y, label = gen_material.abbrv), 
            colour = "white", size =4, 
            position = position_nudge(y=-0.8), fontface = "bold") +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
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
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
dendro.bray.ARG.plot

##RA PLOT - TYPE####
dendro_bray_data2_ARG_order <- data2_ARG.bray.dendro.data$labels$label
data2_ARG.tss.type.melt<- psmelt(data2_ARG.tss.type)
type.filt.palette <- distinctColorPalette(length(unique(data2_ARG.tss.type.melt$Type)))

dendroRA.ARG.type.plot <- ggplot(data2_ARG.tss.type.melt, aes(x=Sample, y= Abundance, fill = Type)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_bray_data2_ARG_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =type.filt.palette) +
  guides(fill=guide_legend(title.position="top", nrow = 1, keywidth = 1, keyheight = 1))+
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
        axis.text.x = element_blank()) 
dendroRA.ARG.type.plot
#RA plot together with dendrogram
dendroRA.ARG.type.plot.2 <- plot_grid(dendro.bray.ARG.plot, dendroRA.ARG.type.plot, 
                                      align = "v", 
                                      ncol = 1,
                                      rel_heights = c(0.3, 0.7))
dendroRA.ARG.type.plot.2

##RA PLOT - CLASS####
data2_ARG.class.filt <- merge_low_abundance_ARG_ra(data2_ARG.tss.class, 
                                                   variable = "sample_type", 
                                                   level = "Class",
                                                   threshold = 1.5)
data2_ARG.class.filt #13 classes with RA over 1.5% (120 samples)
data2_ARG_class.filt.melt <- psmelt(data2_ARG.class.filt)%>%
  mutate(Class = factor(Class, 
                         levels = c(setdiff(Class, 
                                            unique(grep("Others", Class, value = TRUE))), 
                                    unique(grep("Others", Class, value = TRUE)))))##Factoring the order column so that "Others.." is the last category
#Color palette
class.filt.palette <- c(
  "Tetracyclines" = "#AA4499",
  "Drug and biocide resistance" = "#E69F00", 
  "MLS" = "#009E73",
  "Iron resistance" = "#0072B2",
  "betalactams" = "#D55E00",
  "Copper resistance" = "darkred",
  "Aminoglycosides" = "#332288",
  "Multi-metal resistance" = "#88CCEE",
  "Sulfonamides" = "maroon",
  "Rifampin" = "#117733",
  "Multi-biocide resistance" = "#DDCC77",
  "Multi-drug resistance" = "dodgerblue")
class.filt.palette$'Others <1.5% RA' <- "grey95"

##Which are the top most abundant classes by group?
data2_ARG_class.filt.melt %>%
  group_by(sample_type, Class) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(sample_type,  desc(mean_abun))%>%
  print(n=40)

#Per gen_material:
most_abun_classes_ARG <- data2_ARG_class.filt.melt %>%
  group_by(gen_material,sample_type, Class) %>%
  summarise(
  mean_relative_abundance = mean(Abundance, na.rm = TRUE),
  sd = sd(Abundance, na.rm = TRUE),
  .groups = "drop") %>%
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
most_abun_classes_ARG

#####SUPPLEMENTARY TABLE 8.2####
stable8.2 <- most_abun_classes_ARG%>%
  mutate(Dataset = "Resistome")%>%
  select(Dataset, `Library Type`, `Sample Type`, `Class`, `Mean Relative Abundance (%) ± SD`)

#Plot
dendroRA.ARG.class.plot <- ggplot(data2_ARG_class.filt.melt, aes(x=Sample, y= Abundance, fill = Class)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_bray_data2_ARG_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =class.filt.palette) +
  guides(fill = guide_legend(title.position="top", nrow = 2, 
                             keywidth = 1, keyheight = 1))+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 22),
        legend.text = element_text(size = 18),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
dendroRA.ARG.class.plot

#RA and Dendrogram together
dendroRA.ARG.class.plot.2 <- plot_grid(dendro.bray.ARG.plot,
                                       dendroRA.ARG.class.plot, 
                                       align = "v", 
                                       ncol = 1,
                                       rel_heights = c(0.35, 0.65))
dendroRA.ARG.class.plot.2


###FIGURE 4D-E-F - FECES VS CB EMPHASIS#####
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
  labels = c("D", "E"),
  font.label = list(size = 30),
  legend = "none"
)
cDNAandDNA_BC_beta_div_spider_sampletype

##DNA and cDNA ordination plots, on top of dendrogram at the phylum level - Figure 5B&C
figure4DEF <- plot_grid(cDNAandDNA_BC_beta_div_spider_sampletype, 
                          dendroRA.ARG.class.plot.2, 
                      align = "v", 
                      ncol = 1,
                      labels = c(" ", "F"),
                      label_size = 30,
                      rel_heights = c(0.25, 0.75))
  # labs(title = "RESISTOME")+
  # theme(plot.title = element_text(size = 30, face = "bold"))
figure4DEF



##RA PLOT - MECHANISM####
data2_ARG.mechanism.filt <- merge_low_abundance_ARG_ra(data2_ARG.tss.mechanism, 
                                                       variable = "sample_type", 
                                                       level = "Class", 
                                                       threshold = 0.5)
data2_ARG.mechanism.filt #55 mechanisms over 0.5% RA
data2_ARG_mechanism.filt.melt <- psmelt(data2_ARG.mechanism.filt)%>%
  mutate(Mechanism = factor(Mechanism, 
                        levels = c(setdiff(Mechanism, 
                                           unique(grep("Others", Mechanism, value = TRUE))), 
                                   unique(grep("Others", Mechanism, value = TRUE)))))##Factoring the order column so that "Others.." is the last category

#Create base colors based on class palette
base_color_class <- c(class.filt.palette, 
                      "Mupirocin" = "#E51932", #added missing classes 
                      "Tellurium resistance" = "#CCBFFF", #added missing classes 
                      "Others <0.5% RA" = "grey95")

#Make hues based on families within each ammonia-nitrite oxidizing group
palette_ARG_mechanism_df <- data2_ARG_mechanism.filt.melt %>% 
  distinct(Mechanism, Class) %>%
  group_by(Class) %>%
  arrange(Class) %>%   
  mutate(
    base_color = base_color_class[Class],
    shade = seq(0, 0.9, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()

#Set up final palette
mechanism.filt.palette <- setNames(
  palette_ARG_mechanism_df$color,
  palette_ARG_mechanism_df$Mechanism)
mechanism.filt.palette


##Which are the top most abundant taxa by group? 
data2_ARG_mechanism.filt.melt %>%
  group_by(sample_type, Mechanism) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(sample_type,  desc(mean_abun))%>%
  print(n=40)

#Per library type:
most_abun_mechanism <- data2_ARG_mechanism.filt.melt %>%
  group_by(sample_type, gen_material, Mechanism) %>%
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
most_abun_mechanism

# .groups = "drop_last"
# Drops the last level of grouping (default).
# If you grouped by sample_type and Phylum,
# the result will stay grouped by sample_type.

##Apply the function to obtain top mechanisms (n=15)
top_mechanisms <- top_taxa_legend(data2_ARG_mechanism.filt.melt, 
                              taxlevel = "Mechanism", n = 15)
top_mechanisms

#Plot RA Mechanisms
dendroRA.ARG.mechanism.plot <- ggplot(data2_ARG_mechanism.filt.melt, aes(x=Sample, y= Abundance, fill = Mechanism)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_bray_data2_ARG_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =mechanism.filt.palette, 
                    breaks = top_mechanisms) +
  guides(fill = guide_legend(title.position="top", nrow = 6, keywidth = 1, keyheight = 1))+
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
dendroRA.ARG.mechanism.plot

#RA with dendrogram
dendroRA.ARG.mechanism.plot.2 <- plot_grid(dendro.bray.ARG.plot, 
                                           dendroRA.ARG.mechanism.plot, 
                                           align = "v", 
                                           ncol = 1,
                                           rel_heights = c(0.3, 0.7))
dendroRA.ARG.mechanism.plot.2


##RA PLOT - GROUP#### 
data2_ARG.group.filt <- merge_low_abundance_ARG_ra(data2_ARG.tss.group, 
                                                   variable = "sample_type", 
                                                   level = "Group", 
                                                   threshold = 0.4)
data2_ARG.group.filt #59 groups over 0.4% RA in both sample types

#Factoring Group by class, and leaving "Others" last)
data2_ARG_group.filt.melt <- psmelt(data2_ARG.group.filt) %>%
  mutate(Group = factor(Group, 
                        levels = c(
                          unlist(
                            lapply(unique(Class), function(c) {
                              setdiff(unique(Group[Class == c]), unique(grep("Others", Group, value = TRUE)))
                            })
                          ), 
                          unique(grep("Others", Group, value = TRUE))
                        )))

levels(data2_ARG_group.filt.melt$Group) #OK

##Which are the top most abundant taxa by group? 
data2_ARG_group.filt.melt %>%
  group_by(sample_type, Group) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(sample_type,  desc(mean_abun))%>%
  print(n=40)

#Per library typ (40 most abundant ARG groups)
most_abun_group <-data2_ARG_group.filt.melt %>%
  group_by(sample_type, gen_material, Group) %>%
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
most_abun_group
  
#####SUPPLEMENTARY TABLE 9_2####
stable9.2 <- most_abun_group%>%
  mutate(Dataset = "Resistome")%>%
  select(Dataset, `Library Type`, `Sample Type`, `Mean Relative Abundance (%) ± SD`, everything()) 

#Create base colors based on class palette
base_color_class_2 <- c(class.filt.palette, 
                      "Mupirocin" = "#E51932", #added missing classes 
                      "Tellurium resistance" = "purple", #added missing classes 
                      # "Tellurium resistance" = "#CCBFFF", #added missing classes 
                      "Aminocoumarins" = "orange",#added missing classes 
                      "Biocide and metal resistance" = "#FFFF32",#added missing classes 
                      "Elfamycins" = "green",#added missing classes 
                      "Mercury resistance" = "cyan",#added missing classes 
                      "Peroxide resistance" = "darkblue", #added missing classes 
                      "Others <0.4% RA" = "grey95")

#Make hues based on families within each ammonia-nitrite oxidizing group
palette_ARG_group_df <- data2_ARG_group.filt.melt %>% 
  distinct(Group, Class) %>%
  group_by(Class) %>%
  arrange(Class) %>%   
  mutate(
    base_color = base_color_class_2[Class],
    shade = seq(-0.7, 0.5, length.out = n()),
    # shade = (seq(0, 1, length.out = n()) ^ 0.6) * 0.9,
    color = darken(base_color, amount = shade))%>%
  ungroup()

#Set up final palette
group.filt.palette <- setNames(
  palette_ARG_group_df$color,
  palette_ARG_group_df$Group)
group.filt.palette

#PLOT
dendroRA.ARG.group.plot <- ggplot(data2_ARG_group.filt.melt, aes(x=Sample, y= Abundance, fill = Group)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_bray_data2_ARG_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = group.filt.palette)+
  guides(fill = guide_legend(title.position="top", nrow = 3, keywidth = 1, keyheight = 1))+
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
        axis.text.x = element_blank()) 
dendroRA.ARG.group.plot

#RA with dendrogram
dendroRA.ARG.group.plot.2 <- plot_grid(dendro.bray.ARG.plot, 
                                       dendroRA.ARG.group.plot, 
                                       align = "v", 
                                       ncol = 1,
                                       rel_heights = c(0.35, 0.65))
dendroRA.ARG.group.plot.2


###FIGURE 7DEF - EMPHASIS ON NUCLEIC ACID #######
water_genmat_BC_beta_div 
feces_genmat_BC_beta_div 

WaterandFeces_ARG_beta_div_spider_genmat <- ggarrange(
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
  labels = c("D", "E"),
  font.label = list(size = 30),
  legend = "none"
)
WaterandFeces_ARG_beta_div_spider_genmat

##Add to dendrogram
figure7DEF <- plot_grid(WaterandFeces_ARG_beta_div_spider_genmat, 
                           dendroRA.ARG.group.plot.2,
                        align = "v", 
                        ncol = 1,
                        labels = c(" ", "F"),
                        label_size = 30,
                        rel_heights = c(0.25, 0.75))
  # labs(title = "RESISTOME")+
  # theme(plot.title = element_text(size = 30, face = "bold"))
figure7DEF



#DIFFERNTIAL ABUNDANCE ####
##FECES (DNA vs cDNA)#######
###ANCOMBC#######
##Getting untransformed (raw) counts in fecal samples 
data2_ARG.feces <- subset_samples(data2_ARG, sample_type == "Feces") 
data2_ARG.feces <- prune_taxa(taxa_sums(data2_ARG.feces) > 0, data2_ARG.feces) 
data2_ARG.feces ## 2040 taxa and 96 samples 
ancombc_feces.counts <-data2_ARG.feces 
sample_data(ancombc_feces.counts)$feedlot <- factor(sample_data(ancombc_feces.counts)$feedlot, levels = c("1", "2", "3", "4", "5"))
sample_data(ancombc_feces.counts)$gen_material <- factor(sample_data(ancombc_feces.counts)$gen_material, levels = c("DNA", "cDNA"))##reorder gen_material as factor, DNA as "reference"


##Group
# #Preprocessing
data2_ARG.feces.tss <- subset_samples(data2_ARG.tss, sample_type == "Feces") ##Only fecal samples (RA)
data2_ARG.feces.tss <- prune_taxa(taxa_sums(data2_ARG.feces.tss) > 0, data2_ARG.feces.tss)
data2_ARG.feces.tss ##2040 taxa and 96 samples
data2_ARG.tss.group.feces <- tax_glom(data2_ARG.feces.tss, taxrank = "Group", NArm = F) ##Glom to the group level
data2_ARG.tss.group.feces ##672 groups (RA) and 96 samples


##Filtering out the low relative abundance (less than 0.01%) groups
data2_ARG.feces_group.ra.filt <- filter_taxa(data2_ARG.tss.group.feces, 
                                             function(x) mean(x) > 0.01, TRUE)
data2_ARG.feces_group.ra.filt ## 174 genera with mean RA > 0.4% across 96 samples (fecal samples)
##Filtering those genera (> 1% RA) on the raw counts phyloseq object for feces
feces_group.counts_filtered <- subset_taxa(ancombc_feces.counts, 
                                           Group %in% tax_table(data2_ARG.feces_group.ra.filt)[,"Group"])
feces_group.counts_filtered #879 taxa for those 174 genera (96 samples)
##To include feedlot as a random effect, making sure it is a factor (not continuous variable)
sample_data(feces_group.counts_filtered)$feedlot <- factor(sample_data(feces_group.counts_filtered)$feedlot, 
                                                           levels = c("1", "2", "3", "4", "5"))
sample_data(feces_group.counts_filtered)$gen_material <- factor(sample_data(feces_group.counts_filtered)$gen_material,
                                                                levels = c("DNA", "cDNA"))##reorder gen_material as factor, DNA as "reference"


##running ancombc on the variable of interest (gen_material)
ancombc_output_feces.group <-ancombc2(data= feces_group.counts_filtered, 
                                      assay_name = "counts", 
                                      tax_level = "Group",
                                      fix_formula = "gen_material+feedlot",
                                      rand_formula =  "(1 | original_sample)",
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
res.feces.group <- ancombc_output_feces.group$res %>%
  select(-matches("feedlot")) #not doing feedlot comparisons

#Pivot longer the results
ancom_gen_material_feces.group <- res.feces.group %>%
  mutate(across(starts_with("lfc_gen_material"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "comparison", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancom_gen_material_feces.group

##Getting rid of _rounded suffix using sub command
ancom_gen_material_feces.group$comparison ##want to get rid of "_rounded"
ancom_gen_material_feces.group$comparison<- sub("_rounded", "", ancom_gen_material_feces.group$comparison) 
ancom_gen_material_feces.group$comparison #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_gen_material_feces.group <- ancom_gen_material_feces.group %>%
  mutate(comparison= case_when(
    comparison == "gen_materialcDNA" ~ "DNA vs cDNA",
    TRUE ~ comparison ##keeps original name for comparisons not specified (DNA)
  ))
ancom_gen_material_feces.group$comparison ##Now the comparison names are shorter and more manageable
ancom_gen_material_feces.group<- ancom_gen_material_feces.group %>%
  rename(Group = taxon) ##This ancombc was done at the Group level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancom_gen_material_feces.group_2 <- ancom_gen_material_feces.group %>%
  mutate(lower.ci = lfc_gen_materialcDNA - 1.96*se_gen_materialcDNA,
         upper.ci = lfc_gen_materialcDNA + 1.96*se_gen_materialcDNA)

##Final fix - up to make compatible with plotting
ancom_gen_material_feces.group_3 <- ancom_gen_material_feces.group_2 %>%
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
  filter(DA == "q ≤ 0.05")%>%
  select(Group, comparison, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns
nrow(ancom_gen_material_feces.group_3) ##52 DA genera between DNA and cDNA with ANCOM

###MaAsLin3#######
data2_ARG.tss.group.feces ## 672 groups in 96 fecal samples (Relative abundances, non filtered)
data2_ARG.feces_group.ra.filt ## 174 groups in 96 fecal samples (Relative abundances, filtered)
data_maaslin_feces_ra  <- data.frame(t(data2_ARG.feces_group.ra.filt@otu_table),
                                     check.rows = F,
                                     check.names = F)

##Sample metadata
metadata_maaslin_feces <- data.frame(data2_ARG.feces_group.ra.filt@sam_data) %>%
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")),#making feedlot a factor since I'll be adding it as a random effect on the MaAslin model
         gen_material = factor(gen_material, levels = c("DNA", "cDNA"))) #gen_material factor for DNA to be "reference"



##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_feces_group <-  maaslin3(
  input_data = data_maaslin_feces_ra, 
  input_metadata = metadata_maaslin_feces, 
  output = "MaAsLin3_feces",  
  fixed_effect = c("gen_material", 'feedlot'), 
  random_effects = c("original_sample"),
  # random_effects = c("feedlot", "original_sample"),
  # fixed_effect = c("gen_material"), 
  min_prevalence=0.05,
  median_comparison_abundance = T, #Default
  median_comparison_prevalence = FALSE, #default  
  min_abundance = 0, #not filtering by abundance
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_feces_group$fit_data_abundance$results) #Abundance results from MaAslin - 870


##Taxonomy of feces AMR groups
input_taxonomy_feces <- data.frame(data2_ARG.tss.group.feces@tax_table) %>%
  rownames_to_column(var = "feature")

#Calculate confidence intervals, add taxonomy
maaslin_feces_group_2 <- maaslin_feces_group$fit_data_abundance$results %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) %>%
  left_join(input_taxonomy_feces, by = "feature")

##Final edits to put together for plot
maaslin_feces_group_3 <- maaslin_feces_group_2%>%
  rename(comparison = name,
         pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(comparison= case_when(
    comparison == "gen_materialcDNA" ~ "DNA vs cDNA",
    TRUE ~ comparison ##keeps original name for groups not specified 
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
  #filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Group, comparison, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancom_sample_type_DNA.group_3)
nrow(maaslin_feces_group_3) ##68 DA genera between Feces and Water by MaAslin


##ANCOM and MaAslin together
DA_feces_plot_MaAslinANCOM.data <- rbind(ancom_gen_material_feces.group_3, maaslin_feces_group_3) %>%
  filter(Group %in% intersect(maaslin_feces_group_3$Group,
                              ancom_gen_material_feces.group_3$Group)) ##Only going to plot those taxa DA by both tests

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_feces <- maaslin_feces_group$transformed_data %>% #transformed data is log2 transformed TSS normalized otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "feature", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy_feces, by = "feature")%>%
  filter(Group %in% intersect(maaslin_feces_group_3$Group,
                              ancom_gen_material_feces.group_3$Group))%>%
  left_join(metadata_maaslin_feces%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(Group, Sample, logvalue, gen_material, plot)
RA_MaaslinAncom_feces 

#Bias-corrected abundances (ANCOM)
feces.group_log_corr_abn <- ancombc_output_feces.group$bias_correct_log_table %>%
  data.frame()%>% ##make into data frame
  rownames_to_column("Group")%>%
  filter(Group %in% intersect(maaslin_feces_group_3$Group,
                              ancom_gen_material_feces.group_3$Group))%>% #keep only those genera in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -Group, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         gen_material = ifelse(grepl("c", Sample), "cDNA", "DNA"))

#put together objects to plot DA
DA_feces_plot_together <- bind_rows(DA_feces_plot_MaAslinANCOM.data, feces.group_log_corr_abn, RA_MaaslinAncom_feces) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_feces_plot_together$plot <- factor(DA_feces_plot_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_feces_plot_together$gen_material <- factor(DA_feces_plot_together$gen_material, levels = c("DNA", "cDNA"))

###PLOTTING DA#
##Ordering how I want the "group" taxlevel to show up on the plot 
input_taxonomy_feces ##dataframe object for Taxonomy of feces AMR gene groups 

# Create the taxonomy plot data and modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_fecal <- DA_feces_plot_together %>%
  left_join(input_taxonomy_feces, by= "Group")%>%
  distinct(Type, Class, Group) %>%
  arrange(Type, Class) %>%
  mutate(Group = factor(Group, levels = rev(Group)))%>% ##Since I arranged by Class, this is the order I want the genera to show up
  dplyr::group_by(Class) %>%
  dplyr::mutate(label_Class = ifelse(row_number() == 1, Class, "")) %>%  # Create 'label_Class' with only the first occurrence of each class
  ungroup() 

##Factor "Group" level by the order I want (taxonomy_plot_data_fecal$group)
DA_feces_plot_together$Group <- factor(DA_feces_plot_together$Group, levels = rev(taxonomy_plot_data_fecal$Group))


# Create the updated taxonomy plot
taxonomy_plot_fecal <- ggplot(taxonomy_plot_data_fecal) +
  geom_text(aes(x =0, y = Group, 
                label = label_Class), 
            hjust = 0, vjust = 0.5, size = 6.43, family = "sans") +  # Move text left by adjusting x
  labs(title = "Class\n \n ") +
  theme_void() +
  #scale_y_discrete(limits =  rev(taxonomy_plot_data_fecal$group)) +  # Ensure y-axis matches the taxon order
  theme(plot.title = element_text(hjust = 0, size = 20, vjust = -0.8, face = "bold", family = "sans"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )+
  coord_cartesian(xlim = c(0.05, 1)) #had to add this to move the geom_text more to the left
taxonomy_plot_fecal


#Plotting 
DA_feces_plot_MaAslinANCOM <-
  ggplot(data = DA_feces_plot_together) +
  facet_wrap(~ plot, scales='free_x',
             nrow = 1,
             strip.position = "top") +
  geom_boxplot(data=DA_feces_plot_together%>%filter(grepl("abundances", plot)),
               aes(x=Group, y=logvalue, 
                   fill = gen_material, color = gen_material),
               notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(data=DA_feces_plot_together%>%filter(grepl("abundances", plot)),
             aes(x=Group, y=logvalue, fill = gen_material, color = gen_material),
             size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = gen.material.palette,
                     labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
  scale_fill_manual(values=gen.material.palette,
                    labels = c("DNA" = "DNA", "cDNA" = "RNA (cDNA)"))+ 
  guides(color=guide_legend(order = 1,title="Library Type", title.position="top"), 
         fill=guide_legend(order = 1,title="Library Type", title.position="top"))+
  new_scale_color()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_feces_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=Group, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_feces_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=Group, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_feces_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed', 
             alpha=0.5) +
  scale_x_discrete(position='bottom') +
  scale_y_continuous(position='right') +  # Default to break_labels otherwise 
  coord_flip() +
  scale_color_manual(values = c("depleted" = "red", "elevated" = "blue")) +
  scale_linetype_manual(values=c("11", "solid")) +
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Library Type", title.position="top"),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 2.5)),
         linetype = guide_legend(title = "Fold change source", title.position = "top",
                                 override.aes = list(linewidth = 1),
                                 theme = theme(legend.key.width = unit(1.5, "cm"))),
         pch=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  labs(title = "FECES")+
  theme_bw()+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=15), legend.text=element_text(size=14),
        plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y=element_text(size= 20, angle=0, vjust= 1.045, face = "bold"), 
        axis.text.y=element_text(size=18, vjust = 0.5),
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
            aes(x = Group, y = coef, label = DA),
            position = position_dodge2(width = 0.75),
            vjust = -0.6, size = 3)

##Putting together DA plot with taxonomy (family) plot
combined_plot_feces <- plot_grid(
  taxonomy_plot_fecal+ 
    theme(plot.margin = unit(c(0.28, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_feces_plot_MaAslinANCOM + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                     plot.title = element_blank(),
                                     legend.position = "none",
                                     # strip.text = element_blank(),
                                     # strip.background = element_rect(fill = "white")
                                     ),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb")
  # labs(title = "FECES")+
  # theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_feces

##Adding q values 
combined_plot_feces_q <- plot_grid(
  taxonomy_plot_fecal+ theme(plot.margin = unit(c(0.25, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_feces_plot_MaAslinANCOM_q +  theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"), 
                                         plot.title = element_blank(),
                                         legend.position = "none",
                                         strip.text = element_blank(),
                                         strip.background = element_rect(fill = "white")),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") 
  # labs(title = "FECES")+
  # theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_feces_q

##CATCH BASINS (DNA vs cDNA)#####
###ANCOMBC#######
##Getting untransformed (raw) counts in water samples 
data2_ARG.water <- subset_samples(data2_ARG, sample_type == "Water") 
data2_ARG.water <- prune_taxa(taxa_sums(data2_ARG.water) > 0, data2_ARG.water) 
data2_ARG.water ## 2129 taxa and 24 samples 
ancombc_water.counts <-data2_ARG.water 
sample_data(ancombc_water.counts)$feedlot <- factor(sample_data(ancombc_water.counts)$feedlot, levels = c("1", "2", "3", "4", "5"))
sample_data(ancombc_water.counts)$gen_material <- factor(sample_data(ancombc_water.counts)$gen_material, levels = c("DNA", "cDNA"))##reorder gen_material as factor, DNA as "reference"


##Group
# #Preprocessing
data2_ARG.water.tss <- subset_samples(data2_ARG.tss, sample_type == "Water") ##Only water samples (RA)
data2_ARG.water.tss <- prune_taxa(taxa_sums(data2_ARG.water.tss) > 0, data2_ARG.water.tss)
data2_ARG.water.tss ## 2129 taxa and 24 samples  (RA)
data2_ARG.water.tss.group <- tax_glom(data2_ARG.water.tss, taxrank = "Group", NArm = F) ##Glom to the group level
data2_ARG.water.tss.group ##747 groups and 24 samples

##running ancombc on the variable of interest (gen_material)
ancombc_output_water.group <-ancombc2(data= ancombc_water.counts, 
                                      assay_name = "counts", 
                                      tax_level = "Group",
                                      fix_formula = "gen_material+feedlot",
                                      # rand_formula =  "(1 | original_sample)",
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
res.water.group <- ancombc_output_water.group$res %>%
  select(-matches("feedlot")) #not doing feedlot comparisons

#Pivot longer the results
ancom_gen_material_water.group <- res.water.group %>%
  mutate(across(starts_with("lfc_gen_material"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "comparison", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancom_gen_material_water.group

##Getting rid of _rounded suffix using sub command
ancom_gen_material_water.group$comparison ##want to get rid of "_rounded"
ancom_gen_material_water.group$comparison<- sub("_rounded", "", ancom_gen_material_water.group$comparison) 
ancom_gen_material_water.group$comparison #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_gen_material_water.group <- ancom_gen_material_water.group %>%
  mutate(comparison= case_when(
    comparison == "gen_materialcDNA" ~ "DNA vs cDNA",
    TRUE ~ comparison ##keeps original name for comparisons not specified (DNA)
  ))
ancom_gen_material_water.group$comparison ##Now the comparison names are shorter and more manageable
ancom_gen_material_water.group<- ancom_gen_material_water.group %>%
  rename(Group = taxon) ##This ancombc was done at the Group level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancom_gen_material_water.group_2 <- ancom_gen_material_water.group %>%
  mutate(lower.ci = lfc_gen_materialcDNA - 1.96*se_gen_materialcDNA,
         upper.ci = lfc_gen_materialcDNA + 1.96*se_gen_materialcDNA)

##Final fix - up to make compatible with plotting
ancom_gen_material_water.group_3 <- ancom_gen_material_water.group_2 %>%
  filter(passed_ss_gen_materialcDNA==1)%>% #Only want those that passed sensitivity testing
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
      TRUE ~ "Not significant")) %>%
  #filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Group, comparison, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns
nrow(ancom_gen_material_water.group_3) ##4 DA genera between DNA and Water cDNA ANCOM

###MaAsLin3#######
#data2_ARG.water_group.ra.filt #Will feed it this filtered ps object (filtered for those genera with mean RA > 0.3% across 96 samples (water samples))
data2_ARG.water.tss.group  ##747 groups and 24 samples (Relative abundances, non filtered)
data_maaslin_water_ra  <- data.frame(t(data2_ARG.water.tss.group@otu_table),
                                     check.rows = F,
                                     check.names = F)

##Sample metadata
metadata_maaslin_water <- data.frame(data2_ARG.water.tss@sam_data) %>%
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")),#making feedlot a factor since I'll be adding it as a random effect on the MaAslin model
         gen_material = factor(gen_material, levels = c("DNA", "cDNA"))) #gen_material factor for DNA to be "reference"


##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_water_group <-  maaslin3(
  input_data = data_maaslin_water_ra, 
  input_metadata = metadata_maaslin_water, 
  output = "MaAsLin3_water",  
  fixed_effect = c("gen_material", 'feedlot'), 
  random_effects = c("original_sample"),
  # random_effects = c("feedlot", "original_sample"),
  # fixed_effect = c("gen_material"), 
  min_prevalence=0.05,
  median_comparison_abundance = T, #Default
  median_comparison_prevalence = FALSE, #default  
  min_abundance = 0, #not filtering by abundance 
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_water_group$fit_data_abundance$results) #Abundance results from MaAslin - 2955

##Taxonomy of water OTUs 
input_taxonomy_water <- data.frame(data2_ARG.water.tss.group@tax_table) %>%
  rownames_to_column(var = "feature")

#Calculate confidence intervals, add taxonomy
maaslin_water_group_2 <- maaslin_water_group$fit_data_abundance$results %>%
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  left_join(input_taxonomy_water, by = "feature")

##Final edits to put together for plot
maaslin_water_group_3 <- maaslin_water_group_2%>%
  rename(comparison = name,
         pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(comparison= case_when(
    comparison == "gen_materialcDNA" ~ "DNA vs cDNA",
    TRUE ~ comparison ##keeps original name for groups not specified 
  ))%>%
  mutate(direction = ifelse(coef > 0, "elevated", "depleted"))%>%
  mutate(
    plot = "Log2 Fold change with 95%CI", 
    test = "MaAsLin3",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant"))%>%
  # filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Group, comparison, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancom_sample_type_DNA.group_3)
nrow(maaslin_water_group_3) ##12 DA genera between DNA vs cDNA by MaAslin


##ANCOM and MaAslin together
DA_water_plot_MaAslinANCOM.data <- rbind(ancom_gen_material_water.group_3, 
                                         maaslin_water_group_3) %>%
  filter(Group %in% intersect(maaslin_water_group_3$Group,
                              ancom_gen_material_water.group_3$Group)) ##Only going to plot those taxa DA by both tests
DA_water_plot_MaAslinANCOM.data #None

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_water <- maaslin_water_group$transformed_data %>% #transformed data is log2 transformed TSS normalized otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "feature", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy_water, by = "feature")%>%
  filter(Group %in% intersect(maaslin_water_group_3$Group,
                              ancom_gen_material_water.group_3$Group))%>%
  left_join(metadata_maaslin_water%>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(Group, Sample, logvalue, gen_material, plot)

#Bias-corrected abundances (ANCOM)
water.group_log_corr_abn <- ancombc_output_water.group$bias_correct_log_table %>%
  data.frame()%>% ##make into data frame
  rownames_to_column("Group")%>%
  filter(Group %in% intersect(maaslin_water_group_3$Group,
                              ancom_gen_material_water.group_3$Group))%>% #keep only those genera in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -Group, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         gen_material = ifelse(grepl("c", Sample), "cDNA", "DNA"))

####FIGURE 8C - DA EMPHASIS ON DNA vs cDNA#####
#Only feces had Differentially abundant ARG groups
combined_plot_feces

figure8C <-plot_grid(combined_plot_feces, 
                            labels = c("C", " "),
                            label_size = 32,
                            ncol = 1)
  #labs(title = "RESISTOME")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
figure8C


##DNA (CB vs Feces)#####
###ANCOMBC#######
##Getting untransformed (raw) counts in DNA samples 
data2_ARG.DNA <- subset_samples(data2_ARG, gen_material == "DNA") 
data2_ARG.DNA <- prune_taxa(taxa_sums(data2_ARG.DNA) > 0, data2_ARG.DNA) 
data2_ARG.DNA ## 2445 taxa and 60 samples
ancombc_DNA.counts <-data2_ARG.DNA 
##Factoring varible of interest (sample_type), so Feces is the reference
sample_data(ancombc_DNA.counts)$sample_type <- factor(sample_data(ancombc_DNA.counts)$sample_type, levels = c("Feces", "Water"))
##To include feedlot as a random effect, making sure it is a factor (not continuous variable)
sample_data(ancombc_DNA.counts)$feedlot <- factor(sample_data(ancombc_DNA.counts)$feedlot, levels = c("1", "2", "3", "4", "5"))

##GROUP
# #Preprocessing
data2_ARG.DNA.tss <- subset_samples(data2_ARG.tss, gen_material == "DNA") ##Only fecal samples (RA)
data2_ARG.DNA.tss <- prune_taxa(taxa_sums(data2_ARG.DNA.tss) > 0, data2_ARG.DNA.tss)
data2_ARG.DNA.tss ##2445 taxa in 60 DNA samples (RA)
data2_ARG.tss.group.DNA <- tax_glom(data2_ARG.DNA.tss, taxrank = "Group", NArm = F) ##Glom to the group level
data2_ARG.tss.group.DNA ##808 groups (RA) in 60 DNA samples

##Filtering out the low relative abundance (less than 0.01 %) groups
data2_ARG_DNA.group.tss.filt <- filter_taxa(data2_ARG.tss.group.DNA, function(x) mean(x) > 0.01, TRUE)
data2_ARG_DNA.group.tss.filt  ## 220 groups with mean tss > 0.01 % across 60 samples (DNA samples)
##Filtering those genes (> 1% tss) on the raw counts phyloseq object for DNA
DNA_group.counts_filtered <- subset_taxa(ancombc_DNA.counts, 
                                         Group %in% phyloseq::tax_table(data2_ARG_DNA.group.tss.filt)[,"Group"])
DNA_group.counts_filtered ##1333 taxa for those 220 groups (60 samples)
##Factoring varible of interest (sample_type), so Feces is the reference
sample_data(DNA_group.counts_filtered)$sample_type <- factor(sample_data(DNA_group.counts_filtered)$sample_type, levels = c("Feces", "Water"))
##To include feedlot as a tssndom effect, making sure it is a factor (not continuous variable)
sample_data(DNA_group.counts_filtered)$feedlot <- factor(sample_data(DNA_group.counts_filtered)$feedlot, levels = c("1", "2", "3", "4", "5"))

##running ancombc on the variable of interest (sample_type)
ancombc_output_DNA.group <-ancombc2(data= DNA_group.counts_filtered, 
                                    assay_name = "counts", 
                                    tax_level = "Group",
                                    fix_formula = "sample_type + feedlot", 
                                    # rand_formula =  "(1 | feedlot)", 
                                    prv_cut = 0.06, 
                                    lib_cut = 0, 
                                    group= "sample_type", 
                                    struc_zero = T, 
                                    neg_lb = T,
                                    alpha = 0.05, #default significance
                                    n_cl = 1, verbose = TRUE)

## extract results from comparisons 
res.DNA_AMR.group <- ancombc_output_DNA.group$res %>%
  select(-matches("feedlot")) #not doing feedlot comparisons

#Pivot longer the results
ancom_sample_type_DNA.group<- res.DNA_AMR.group %>%
  mutate(across(starts_with("lfc_sample_type"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "comparison", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancom_sample_type_DNA.group

##Getting rid of _rounded suffix using sub command
ancom_sample_type_DNA.group$comparison ##want to get rid of "_rounded"
ancom_sample_type_DNA.group$comparison<- sub("_rounded", "", ancom_sample_type_DNA.group$comparison) 
ancom_sample_type_DNA.group$comparison #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_sample_type_DNA.group <- ancom_sample_type_DNA.group %>%
  mutate(comparison= case_when(
    comparison == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ comparison ##keeps original name for groups not specified (DNA)
  ))
ancom_sample_type_DNA.group$comparison ##Now the group names are shorter and more manageable
ancom_sample_type_DNA.group <- ancom_sample_type_DNA.group %>%
  rename(Group = taxon) ##This ancombc was done at the group level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancom_sample_type_DNA.group_2 <- ancom_sample_type_DNA.group %>%
  mutate(lower.ci = lfc_sample_typeWater - 1.96*se_sample_typeWater,
         upper.ci = lfc_sample_typeWater + 1.96*se_sample_typeWater)


##Final fix - up to make compatible with plotting
ancom_sample_type_DNA.group_3 <- ancom_sample_type_DNA.group_2 %>%
  filter (passed_ss_sample_typeWater == 1)%>%##Only want those that passed sensitivity testing
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
  #filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Group, comparison, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
nrow(ancom_sample_type_DNA.group_3) ##92 DA genera between Feces and Water with ANCOM

###MaAsLin3#######
#Data (otu counts) and metadata for MaAslin
#Counts for data_maaslin_DNA_ra are tss transformed
#data2_ARG.tss.group.DNA.filt  ## 49 groups with mean RA > 0.05% across 60 samples (DNA samples) 
data2_ARG.tss.group.DNA ##808 groups and 60 DNA samples (RA)
data2_ARG_DNA.group.tss.filt##220 groups and 60 samples (RA filtered)

data_maaslin_DNA_ra  <- data.frame(t(data2_ARG_DNA.group.tss.filt@otu_table),
                                   check.rows = F,
                                   check.names = F) 

##Sample metadata
metadata_maaslin_DNA <- data.frame(data2_ARG_DNA.group.tss.filt@sam_data) %>%
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")),
         sample_type = factor(sample_type, levels = c("Feces", "Water"))) #making feedlot a factor since I'll be adding it as a random effect on the MaAslin model


##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_DNA_group <-  maaslin3(
  input_data = data_maaslin_DNA_ra, 
  input_metadata = metadata_maaslin_DNA, 
  output = "MaAsLin3_DNA",  
  fixed_effect = c("sample_type", "feedlot"),
  # fixed_effect = "sample_type",
  # random_effects = "feedlot",
  min_prevalence=0.05,
  median_comparison_abundance = F, #not default
  median_comparison_prevalence = F, #default 
  min_abundance = 0, ##input_data not going to filter
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_DNA_group$fit_data_abundance$results) #Abundance results from MaAslin - 220


##Taxonomy of DNA OTUs 
input_taxonomy_DNA <- data.frame(data2_ARG.tss.group.DNA@tax_table) %>%
  rownames_to_column(var = "feature")

#Calculate confidence intervals, add taxonomy
maaslin_DNA_group_2 <- maaslin_DNA_group$fit_data_abundance$results %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) %>%
  left_join(input_taxonomy_DNA, by = "feature")

##Final edits to put together for plot
maaslin_DNA_group_3 <- maaslin_DNA_group_2%>%
  rename(comparison = name,
         pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(comparison= case_when(
    comparison == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ comparison ##keeps original name for groups not specified 
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
  #filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Group, comparison, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancom_sample_type_DNA.group_3)
nrow(maaslin_DNA_group_3) ##143 DA genera between Feces and Water by MaAslin

##ANCOM and MaAslin together
DA_DNA_plot_MaAslinANCOM.data <- rbind(ancom_sample_type_DNA.group_3, maaslin_DNA_group_3) %>%
  filter(Group %in% intersect(maaslin_DNA_group_3$Group,
                              ancom_sample_type_DNA.group_3$Group))

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_DNA <- maaslin_DNA_group$transformed_data %>% #transformed data is TSS (ra) log2 transformed otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "feature", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy_DNA, by = "feature")%>%
  filter(Group %in% intersect(maaslin_DNA_group_3$Group,
                              ancom_sample_type_DNA.group_3$Group))%>%
  left_join(metadata_maaslin_DNA %>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(Group, Sample, logvalue, sample_type, plot)
RA_MaaslinAncom_DNA
#Bias-corrected abundances (ANCOM)
DNA.group_log_corr_abn <- ancombc_output_DNA.group$bias_correct_log_table %>%
  data.frame()%>% ##make into data frame
  rownames_to_column("Group")%>%
  filter(Group %in% intersect(maaslin_DNA_group_3$Group,
                              ancom_sample_type_DNA.group_3$Group))%>% #keep only those genera in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -Group, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         sample_type = ifelse(grepl("W", Sample), "Water", "Feces"))

#put together objects to plot DA
DA_DNA_plot_together <- bind_rows(DA_DNA_plot_MaAslinANCOM.data, DNA.group_log_corr_abn, RA_MaaslinAncom_DNA) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_DNA_plot_together$plot <- factor(DA_DNA_plot_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_DNA_plot_together$sample_type <- factor(DA_DNA_plot_together$sample_type, 
                                           levels = c("Feces", "Water"))
#Pick the top 60 coef 
top_DA_DNA <- DA_DNA_plot_together %>%
  filter(test %in% c("ANCOM-BC", "MaAsLin3")) %>%
  group_by(test) %>%
  arrange(desc(coef)) %>%
  slice_head(n = 20) %>%        # top positive
  bind_rows(DA_DNA_plot_together %>%
      filter(test %in% c("ANCOM-BC", "MaAsLin3")) %>%
      group_by(test) %>%
      arrange(coef) %>%
      slice_head(n = 20)         # top negative
  ) %>%
  distinct() %>%                # safety check
  ungroup()
top_DA_DNA

#Keep only those
DA_DNA_plot_together <- DA_DNA_plot_together%>%
  filter(Group %in% unique(top_DA_DNA$Group))

###PLOTTING DA#
##Ordering how I want the "Group" taxlevel to show up on the plot 
input_taxonomy_DNA ##dataframe object for Taxonomy of DNA OTUs 

# Create the taxonomy plot data and modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_DNA <- DA_DNA_plot_together %>%
  left_join(input_taxonomy_DNA, by= "Group")%>%
  distinct(Type, Class, Group) %>%
  arrange(Type, Class) %>%
  mutate(Group = factor(Group, levels = rev(Group)))%>% ##Since I arranged by Class, this is the order I want the groups to show up
  dplyr::group_by(Class) %>%
  dplyr::mutate(label_Class = ifelse(row_number() == 1, Class, "")) %>%  # Create 'label_Class' with only the first occurrence of each class
  ungroup() 

levels(taxonomy_plot_data_DNA$Group)

##Factor "Group" level by the order I want (taxonomy_plot_data_DNA$Group)
DA_DNA_plot_together$Group <- factor(DA_DNA_plot_together$Group, levels = rev(taxonomy_plot_data_DNA$Group))

# Create the updated taxonomy plot
taxonomy_plot_DNA <- ggplot(taxonomy_plot_data_DNA) +
  geom_text(aes(x =0, y = Group, label = label_Class), hjust = 0, vjust = 0.5, size = 7.14, family = "sans") +  # Move text left by adjusting x
  labs(title = "Class") +
  theme_void() +
  #scale_y_discrete(limits =  rev(taxonomy_plot_data_DNA$Group)) +  # Ensure y-axis matches the taxon order
  theme(plot.title = element_text(hjust = 0, size = 22, vjust = -0.3, face = "bold", family = "sans"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )+
  coord_cartesian(xlim = c(0.05, 1)) #had to add this to move the geom_text more to the left
taxonomy_plot_DNA

#Plotting 
DA_DNA_plot_MaAslinANCOM <-
  ggplot(data = DA_DNA_plot_together) +
  facet_wrap(~ plot, scales='free_x',
             nrow = 1,
             strip.position = "top") +
  geom_boxplot(data=DA_DNA_plot_together%>%filter(grepl("abundances", plot)),
               aes(x=Group, y=logvalue, 
                   fill = sample_type, color = sample_type),
               notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(data=DA_DNA_plot_together%>%filter(grepl("abundances", plot)),
             aes(x=Group, y=logvalue, fill = sample_type, color = sample_type),
             size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  # facet_nested(. ~ plot, scales='free_x',
  #              space='free_y',
  #              switch='y',
  #              strip=strip_nested(text_y=list(element_text(angle=0))),
  #              labeller=labeller(group=label_wrap_gen(width=10),
  #                                sub_group=label_wrap_gen(width=10))) +
  # ggplot(data=DA_DNA_plot_together%>%filter(grepl("abundances", plot)),
  #        aes(x=Group, y=logvalue, fill = sample_type, color = sample_type)) +
  # geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  # geom_point(size = 1.5, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order = 1,title="Sample type", title.position="top"),
         fill=guide_legend(order = 1,title="Sample type", title.position="top"))+
  scale_fill_manual(values=sample.type.palette, 
                    labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  scale_color_manual(values = sample.type.palette, 
                     labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  new_scale_color()+
  new_scale_fill()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_DNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=Group, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_DNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=Group, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_DNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed'
             , alpha=0.5) +
  scale_x_discrete(position='bottom') +
  scale_y_continuous(position='right') +  # Default to break_labels otherwise 
  coord_flip() +
  scale_fill_manual(values=sample.type.palette) +
  scale_color_manual(values=c("red", "blue")) +
  scale_linetype_manual(values=c("11", "solid")) +
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Sample Type", title.position="top"),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 2.5)),
         linetype = guide_legend(title = "Fold change source", title.position = "top",
                                 override.aes = list(linewidth = 1),
                                 theme = theme(legend.key.width = unit(1.5, "cm"))),
         pch=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  labs(title = "Metagenomic libraries (DNA)")+
  theme_bw()+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=15), legend.text=element_text(size=14),
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

DA_DNA_plot_MaAslinANCOM
#There's a warning here. It is because of those samples where the group is not present at all. 
#For example: data2_ARG.group@otu_table["MEG_7112", "F2F02"]. MEG_71112 (representative of the TETM group) is not present at all in F2F02. So won't get an abundance value out of it

##Adding the q values
DA_DNA_plot_MaAslinANCOM_q <- DA_DNA_plot_MaAslinANCOM +
  geom_text(inherit.aes=FALSE,
            data = DA_DNA_plot_together %>% filter(plot == 'Log2 Fold change with 95%CI'),
            aes(x = Group, y = coef, label = DA),
            position = position_dodge2(width = 0.75),
            vjust = -0.6, size = 3)

##Putting together DA plot with taxonomy (family) plot
combined_plot_DNA <- plot_grid(
  taxonomy_plot_DNA+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_DNA_plot_MaAslinANCOM + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                   plot.title = element_blank(),
                                   legend.position = "none",
                                   # strip.text = element_blank(),
                                   # strip.background = element_rect(fill = "white")
                                   ),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "Metagenomic libraries (DNA)")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_DNA

##Adding q values 
combined_plot_DNA_q <- plot_grid(
  taxonomy_plot_DNA+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_DNA_plot_MaAslinANCOM_q + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                     plot.title = element_blank(),
                                     strip.text=element_text(size=16, color = "white", face = "bold")),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "Metagenomic libraries (DNA)")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_DNA_q


##cDNA (CB vs Feces)#####
###ANCOMBC#######
##Getting untransformed (raw) counts in cDNA samples 
data2_ARG.cDNA <- subset_samples(data2_ARG, gen_material == "cDNA") 
data2_ARG.cDNA <- prune_taxa(taxa_sums(data2_ARG.cDNA) > 0, data2_ARG.cDNA) 
data2_ARG.cDNA ## 2287 taxa and 60 samples (cDNA)
ancombc_cDNA.counts <-data2_ARG.cDNA 
##Factoring varible of interest (sample_type), so Feces is the reference
sample_data(ancombc_cDNA.counts)$sample_type <- factor(sample_data(ancombc_cDNA.counts)$sample_type, levels = c("Feces", "Water"))
##To include feedlot as a random effect, making sure it is a factor (not continuous variable)
sample_data(ancombc_cDNA.counts)$feedlot <- factor(sample_data(ancombc_cDNA.counts)$feedlot, levels = c("1", "2", "3", "4", "5"))

# #Preprocessing
data2_ARG.cDNA.tss <- subset_samples(data2_ARG.tss, gen_material == "cDNA") ##Only fecal samples (RA)
data2_ARG.cDNA.tss <- prune_taxa(taxa_sums(data2_ARG.cDNA.tss) > 0, data2_ARG.cDNA.tss)
data2_ARG.cDNA.tss ##2287 taxa and 60 cDNA samples  (RA)
data2_ARG.tss.group.cDNA <- tax_glom(data2_ARG.cDNA.tss, taxrank = "Group", NArm = F) ##Glom to the group level
data2_ARG.tss.group.cDNA ##761 groups (RA) in 60 cDNA samples (RA)

##Filtering out the low relative abundance (less than 0.05 %) groups
data2_ARG_cDNA.group.tss.filt <- filter_taxa(data2_ARG.tss.group.cDNA, function(x) mean(x) > 0.05, TRUE)
data2_ARG_cDNA.group.tss.filt  ## 104 groups with mean tss > 0.05% across 60 samples (cDNA samples)
##Filtering those genes (> 0.05% tss) on the raw counts phyloseq object for cDNA
cDNA_group.counts_filtered <- subset_taxa(ancombc_cDNA.counts, Group %in% phyloseq::tax_table(data2_ARG_cDNA.group.tss.filt)[,"Group"])
cDNA_group.counts_filtered ##848 taxa for those 104 groups (60 samples)
##Factoring varible of interest (sample_type), so Feces is the reference
sample_data(cDNA_group.counts_filtered)$sample_type <- factor(sample_data(cDNA_group.counts_filtered)$sample_type, levels = c("Feces", "Water"))
##To include feedlot as a tssndom effect, making sure it is a factor (not continuous variable)
sample_data(cDNA_group.counts_filtered)$feedlot <- factor(sample_data(cDNA_group.counts_filtered)$feedlot, levels = c("1", "2", "3", "4", "5"))

##running ancombc on the variable of interest (sample_type)
ancombc_output_cDNA.group <-ancombc2(data= cDNA_group.counts_filtered, 
                                    assay_name = "counts", 
                                    tax_level = "Group",
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
res.cDNA_AMR.group <- ancombc_output_cDNA.group$res%>%
  select(-matches("feedlot")) #not doing feedlot comparisons 

#Pivot longer the results
ancom_sample_type_cDNA.group<- res.cDNA_AMR.group %>%
  mutate(across(starts_with("lfc_sample_type"),
                ~round(.x, 2), .names = "{.col}_rounded"))%>% ##Round the log-fold change of abundance to 2 decimal places
  pivot_longer(cols=contains("_rounded"), names_to = "comparison", values_to = "value", names_prefix= "lfc_")%>% ##Pivot to long format
  arrange(taxon)
ancom_sample_type_cDNA.group

##Getting rid of _rounded suffix using sub command
ancom_sample_type_cDNA.group$comparison ##want to get rid of "_rounded"
ancom_sample_type_cDNA.group$comparison<- sub("_rounded", "", ancom_sample_type_cDNA.group$comparison) 
ancom_sample_type_cDNA.group$comparison #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_sample_type_cDNA.group <- ancom_sample_type_cDNA.group %>%
  mutate(comparison= case_when(
    comparison == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ comparison ##keeps original name for groups not specified (cDNA)
  ))
ancom_sample_type_cDNA.group$comparison ##Now the group names are shorter and more manageable
ancom_sample_type_cDNA.group <- ancom_sample_type_cDNA.group %>%
  rename(Group = taxon) ##This ancombc was done at the group level, changing the column name to that

##Obtaining the confidence intervals for the log fold change
ancom_sample_type_cDNA.group_2 <- ancom_sample_type_cDNA.group %>%
  mutate(lower.ci = lfc_sample_typeWater - 1.96*se_sample_typeWater,
         upper.ci = lfc_sample_typeWater + 1.96*se_sample_typeWater)

##Final fix - up to make compatible with plotting
ancom_sample_type_cDNA.group_3 <- ancom_sample_type_cDNA.group_2 %>%
  filter(passed_ss_sample_typeWater==1)%>% #Only want those that passed sensitivity testing
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
  #filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Group, comparison, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
nrow(ancom_sample_type_cDNA.group_3) ##52 DA genera between Feces and Water with ANCOM

###MaAsLin3#######
#Data (otu counts) and metadata for MaAslin
#Counts for data_maaslin_cDNA_ra are tss transformed
# data2_ARG_cDNA.group.ra.filt  ## 72 groups with mean RA > 0.05% across 60 samples (cDNA samples) 
data2_ARG.tss.group.cDNA #761 groups and 60 samples  cDNA samples (RA)
data2_ARG_cDNA.group.tss.filt #104 groups and 60 samples  cDNA samples (RA)
data_maaslin_cDNA_ra  <- data.frame(t(data2_ARG_cDNA.group.tss.filt@otu_table), 
                                    check.rows = F,
                                    check.names = F)  


##Sample metadata
metadata_maaslin_cDNA <- data.frame(data2_ARG_cDNA.group.tss.filt@sam_data) %>%
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")),
         sample_type = factor(sample_type, levels = c("Feces", "Water"))) #making feedlot a factor since I'll be adding it as a random effect on the MaAslin model

##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_cDNA_group <-  maaslin3(
  input_data = data_maaslin_cDNA_ra, 
  input_metadata = metadata_maaslin_cDNA, 
  output = "MaAsLin3_cDNA",  
  fixed_effect = c("sample_type", "feedlot"),
  # fixed_effect = "sample_type",
  # random_effects = "feedlot",
  min_prevalence=0.05,
  median_comparison_abundance = T, #Default
  median_comparison_prevalence = FALSE, #default 
  min_abundance = 0, ##not filtering input data
  normalization="NONE", #input_data is already normalized
  transform='LOG', #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
nrow(maaslin_cDNA_group$fit_data_abundance$results) #Abundance results from MaAslin - 520

##Taxonomy of cDNA OTUs 
input_taxonomy_cDNA <- data.frame(data2_ARG_cDNA.group.tss.filt@tax_table) %>%
  rownames_to_column(var = "feature")

#Calculate confidence intervals, add taxonomy
maaslin_cDNA_group_2 <- maaslin_cDNA_group$fit_data_abundance$results %>%
  filter(metadata != "feedlot")%>% #not doing feedlot comparisons
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr) %>%
  left_join(input_taxonomy_cDNA, by = "feature")

##Final edits to put together for plot
maaslin_cDNA_group_3 <- maaslin_cDNA_group_2%>%
  rename(comparison = name,
         pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(comparison= case_when(
    comparison == "sample_typeWater" ~ "Feces vs Water",
    TRUE ~ comparison ##keeps original name for groups not specified 
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
  #filter(DA == "q ≤ 0.05" | DA == "q ≤ 0.1")%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Group, comparison, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) ##Only need these columns (same as ancom_sample_type_cDNA.group_3)
nrow(maaslin_cDNA_group_3) ##61 DA genera between Feces and Water by MaAslin

##ANCOM and MaAslin together
DA_cDNA_plot_MaAslinANCOM.data <- rbind(ancom_sample_type_cDNA.group_3, maaslin_cDNA_group_3) %>%
  filter(Group %in% intersect(maaslin_cDNA_group_3$Group,
                              ancom_sample_type_cDNA.group_3$Group)) ##Only going to plot those taxa DA by both tests

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_cDNA <- maaslin_cDNA_group$transformed_data %>% #transformed data is TSS (ra) log2 transformed otu counts from MaAslin output
  rownames_to_column(var = "Sample")%>%
  pivot_longer(cols = -Sample, names_to = "feature", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy_cDNA, by = "feature")%>%
  filter(Group %in% intersect(maaslin_cDNA_group_3$Group,
                              ancom_sample_type_cDNA.group_3$Group))%>%
  left_join(metadata_maaslin_cDNA %>% rownames_to_column(var = "Sample"), by = "Sample")%>%
  select(Group, Sample, logvalue, sample_type, plot)
RA_MaaslinAncom_cDNA
#Bias-corrected abundances (ANCOM)
cDNA.group_log_corr_abn <- ancombc_output_cDNA.group$bias_correct_log_table %>%
  data.frame()%>% ##make into data frame
  rownames_to_column("Group")%>%
  filter(Group %in% intersect(maaslin_cDNA_group_3$Group,
                              ancom_sample_type_cDNA.group_3$Group))%>% #keep only those genera in both the ANCOM and MaAslin outputs
  pivot_longer(cols = -Group, names_to = "Sample", values_to = "logvalue")%>% ##Pivot to longer data format
  mutate(plot = 'log(Bias-corrected abundances)',
         sample_type = ifelse(grepl("W", Sample), "Water", "Feces"))

#put together objects to plot DA
DA_cDNA_plot_together <- bind_rows(DA_cDNA_plot_MaAslinANCOM.data, cDNA.group_log_corr_abn, RA_MaaslinAncom_cDNA) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_cDNA_plot_together$plot <- factor(DA_cDNA_plot_together$plot, levels = c("Log2 Fold change with 95%CI", "log(Bias-corrected abundances)", 'log2(Relative abundances)'))
DA_cDNA_plot_together$sample_type <- factor(DA_cDNA_plot_together$sample_type, levels = c("Feces", "Water"))

#Pick the top 60 coef 
top_DA_cDNA <- DA_cDNA_plot_together %>%
  filter(test %in% c("ANCOM-BC", "MaAsLin3")) %>%
  group_by(test) %>%
  arrange(desc(coef)) %>%
  slice_head(n = 20) %>%        # top positive
  bind_rows(DA_cDNA_plot_together %>%
              filter(test %in% c("ANCOM-BC", "MaAsLin3")) %>%
              group_by(test) %>%
              arrange(coef) %>%
              slice_head(n = 20)         # top negative
  ) %>%
  distinct() %>%                # safety check
  ungroup()
top_DA_cDNA
#Keep only those
DA_cDNA_plot_together <- DA_cDNA_plot_together%>%
  filter(Group %in% unique(top_DA_cDNA$Group))

###PLOTTING DA#
##Ordering how I want the "Group" taxlevel to show up on the plot 
input_taxonomy_cDNA ##dataframe object for Taxonomy of cDNA OTUs 

# Create the taxonomy plot data and modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_cDNA <- DA_cDNA_plot_together %>%
  left_join(input_taxonomy_cDNA, by= "Group")%>%
  distinct(Type, Class, Group) %>%
  arrange(Type, Class) %>%
  mutate(Group = factor(Group, levels = rev(Group)))%>% ##Since I arranged by Class, this is the order I want the groups to show up
  dplyr::group_by(Class) %>%
  dplyr::mutate(label_Class = ifelse(row_number() == 1, Class, "")) %>%  # Create 'label_Class' with only the first occurrence of each class
  ungroup() 

##Factor "Group" level by the order I want (taxonomy_plot_data_cDNA$Group)
DA_cDNA_plot_together$Group <- factor(DA_cDNA_plot_together$Group, levels = rev(taxonomy_plot_data_cDNA$Group))

# Create the updated taxonomy plot
taxonomy_plot_cDNA <- ggplot(taxonomy_plot_data_cDNA) +
  geom_text(aes(x =0, y = Group, label = label_Class), hjust = 0, vjust = 0.5,size = 7.14, family = "sans") +  # Move text left by adjusting x
  labs(title = "Class") +
  theme_void() +
  #scale_y_discrete(limits =  rev(taxonomy_plot_data_cDNA$Group)) +  # Ensure y-axis matches the taxon order
  theme(plot.title = element_text(hjust = 0, size = 22, vjust = -0.3, face = "bold", family = "sans"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )+
  coord_cartesian(xlim = c(0.05, 1)) #had to add this to move the geom_text more to the left
taxonomy_plot_cDNA

#Plotting 
DA_cDNA_plot_MaAslinANCOM <-
  ggplot(data = DA_cDNA_plot_together) +
  facet_wrap(~ plot, scales='free_x',
             nrow = 1,
             strip.position = "top") +
  geom_boxplot(data=DA_cDNA_plot_together%>%filter(grepl("abundances", plot)),
               aes(x=Group, y=logvalue, 
                   fill = sample_type, color = sample_type),
               notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  geom_point(data=DA_cDNA_plot_together%>%filter(grepl("abundances", plot)),
             aes(x=Group, y=logvalue, fill = sample_type, color = sample_type),
             size =2, shape = 18, position = position_dodge(width = 0.75)) +
  # facet_nested(. ~ plot, scales='free_x',
  #              space='free_y',
  #              switch='y',
  #              strip=strip_nested(text_y=list(element_text(angle=0))),
  #              labeller=labeller(group=label_wrap_gen(width=10),
  #                                sub_group=label_wrap_gen(width=10))) +
  # ggplot(data=DA_cDNA_plot_together%>%filter(grepl("abundances", plot)),
  #        aes(x=Group, y=logvalue, fill = sample_type, color = sample_type)) +
  # geom_boxplot(notch=FALSE, outlier.size=0.5, size = 0.5, alpha = 0.3) +
  # geom_point(size = 1.5, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order = 1,title="Sample type", title.position="top"),
         fill=guide_legend(order = 1,title="Sample type", title.position="top"))+
  scale_fill_manual(values=sample.type.palette, 
                    labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  scale_color_manual(values = sample.type.palette, 
                     labels = c("Water" = "Catch Basins", "Feces" = "Feces")) +
  new_scale_color()+
  new_scale_fill()+
  geom_errorbar(inherit.aes=FALSE,
                data=DA_cDNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
                aes(x=Group, ymin=lower.ci, ymax=upper.ci,
                    color=direction, linetype=test),
                width=0, position=position_dodge(0.75), size=0.75) +
  geom_point(inherit.aes=FALSE,
             data=DA_cDNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(x=Group, y=coef, color=direction, pch=test),
             position=position_dodge(0.75), size=3) +
  geom_hline(data=DA_cDNA_plot_together%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(yintercept=0),
             size=0.5, linetype='dashed'
             , alpha=0.5) +
  scale_x_discrete(position='bottom') +
  scale_y_continuous(position='right') +  # Default to break_labels otherwise 
  coord_flip() +
  scale_fill_manual(values=sample.type.palette) +
  scale_color_manual(values=c("red", "blue")) +
  scale_linetype_manual(values=c("11", "solid")) +
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Sample Type", title.position="top"),
         color=guide_legend(order=2, title="Fold change direction", 
                            title.position="top", override.aes = list(size = 2.5)),
         linetype = guide_legend(title = "Fold change source", title.position = "top",
                                 override.aes = list(linewidth = 1),
                                 theme = theme(legend.key.width = unit(1.5, "cm"))),
         pch=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 3))) +
  labs(title = "Metatranscriptomic libraries (RNA (cDNA))")+
  theme_bw()+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=15), legend.text=element_text(size=14),
        plot.title = element_text(size = 30, colour = "black", face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        axis.title.x=element_blank(), 
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size= 22, angle=0, vjust= 1.045, face = "bold"), 
        axis.text.y=element_text(size=20, vjust = 0.5),
        strip.text=element_text(size=16, color = "white", face = "bold"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_cDNA_plot_MaAslinANCOM

##Adding the q values
DA_cDNA_plot_MaAslinANCOM_q <- DA_cDNA_plot_MaAslinANCOM +
  geom_text(inherit.aes=FALSE,
            data = DA_cDNA_plot_together %>% filter(plot == 'Log2 Fold change with 95%CI'),
            aes(x = Group, y = coef, label = DA),
            position = position_dodge2(width = 0.75),
            vjust = -0.6, size = 3)


##Putting together DA plot with taxonomy plot
combined_plot_cDNA <- plot_grid(
  taxonomy_plot_cDNA+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_cDNA_plot_MaAslinANCOM  + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                     plot.title = element_blank(),
                                     legend.position = "none",
                                     # strip.text = element_blank(),
                                     # strip.background = element_rect(fill = "white")
                                     ),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "Metatranscriptomic libraries (RNA (cDNA))")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_cDNA


##Adding q values 
combined_plot_cDNA_q <- plot_grid(
  taxonomy_plot_cDNA+ theme(plot.margin = unit(c(0.2, 0, 0.2, 0.5), "cm")),  # Change margins, need less right margin
  DA_cDNA_plot_MaAslinANCOM_q + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm"),
                                      plot.title = element_blank(),
                                      strip.text=element_text(size=16, color = "white", face = "bold")),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb") +
  labs(title = "Metatranscriptomic libraries (RNA (cDNA))")+
  theme(plot.title = element_text(size = 30, colour = "black", face = "bold"))
combined_plot_cDNA_q



####FIGURE 5CandD#####
#Put together these 2 plots:
combined_plot_DNA
combined_plot_cDNA

figure5CandD <-plot_grid(combined_plot_DNA+
                           theme(plot.title = element_blank()), 
                         combined_plot_cDNA+
                           theme(plot.title = element_blank()), 
                         align = "v",
                         labels = c("C", "D"),
                         # rel_heights = c(1.05,0.95),
                         label_size = 32,
                         ncol = 1)
  # labs(title = "RESISTOME")+
  # theme(plot.title = element_text(size = 40, face = "bold"))
figure5CandD


#UPSET PLOT #####
##DNA (CB vs. Feces)-Figure 4C#####
##making the dataset for upset of UpSetR (binary matrix, present or absent)
data2_ARG.DNA.group <- tax_glom (data2_ARG.DNA, taxrank = "Group", NArm = F)
data2_ARG.DNA.group ##808 groups  are present in 60 DNA samples
upset.data2_ARG.DNA.group <- MicrobiotaProcess::get_upset(data2_ARG.DNA.group, factorNames="sample_type") 
upset.data2_ARG.DNA.group

upset.data2_ARG.DNA.group <- upset.data2_ARG.DNA.group%>%
  rename("CB" = "Water")

#Plot
upset_plot_ARG_DNA <- UpSetR::upset(upset.data2_ARG.DNA.group,
              sets.bar.color = c("#4C72B0", 
                                 "brown"),
              order.by = "freq", text.scale = c(3, 2.5, 2.5, 1.5, 3, 2), 
              point.size = 5, line.size = 2, mainbar.y.label= "Group count",
              sets.x.label = "Group count", 
              set_size.show = F) 
upset_plot_ARG_DNA

##saving the upset plot
svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure4C.svg",
        width=7, height=6)
upset_plot_ARG_DNA
dev.off()

##cDNA (CB vs. Feces)-Figure 4D#####
##making the dataset for upset of UpSetR (binary matrix, present or absent)
data2_ARG.cDNA.group <- tax_glom (data2_ARG.cDNA, taxrank = "Group", NArm = F)
data2_ARG.cDNA.group ##300 groups are present in cDNA samples
upset.data2_ARG.cDNA.group <- MicrobiotaProcess::get_upset(data2_ARG.cDNA.group, factorNames="sample_type") 
upset.data2_ARG.cDNA.group
upset.data2_ARG.cDNA.group <- upset.data2_ARG.cDNA.group%>%
  rename("CB" = "Water")

#Plot
upset_plot_ARG_cDNA <- UpSetR::upset(
  upset.data2_ARG.cDNA.group,
  sets.bar.color = c("#4C72B0", 
                     "brown"),
  order.by = "freq",
  text.scale = c(3, 2.5, 2.5, 1.5, 3, 2),
  point.size = 5,
  line.size = 2,
  mainbar.y.label = "Group count",
  sets.x.label = "Group count",
  set_size.show = FALSE)
upset_plot_ARG_cDNA 

##saving the upset plot
svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure4D.svg",
        width=7, height=6)
upset_plot_ARG_cDNA 
dev.off()


##Feces (cDNA vs. DNA) - Figure 8C#####
##making the dataset for upset of UpSetR (binary matrix, present or absent)
data2_ARG.feces.group <- tax_glom (data2_ARG.feces, taxrank = "Group", NArm = F)
data2_ARG.feces.group ##279 groups are present in feces samples (96 samples)
upset.data2_ARG.feces.group <- MicrobiotaProcess::get_upset(data2_ARG.feces.group, factorNames="gen_material") 
upset.data2_ARG.feces.group

upset.data2_ARG.feces.group <- upset.data2_ARG.feces.group %>%
  rename("RNA(cDNA)" = cDNA)

#Plot
upset_plot_ARG_feces<- UpSetR::upset(upset.data2_ARG.feces.group,
              sets.bar.color = c("#CC79A7", "#009E73"), 
              order.by = "freq", text.scale = c(3, 2.5, 2.5, 2, 3, 2), 
              point.size = 5, line.size = 2, mainbar.y.label= "Group count",
              sets.x.label = "Group count", 
              set_size.show = F) 
upset_plot_ARG_feces

##saving the upset plot - Figure 8C
svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure8C.svg",
        width=7, height=6)
upset_plot_ARG_feces
dev.off()


##CATCH BASINS (cDNA vs. DNA) - Figure 8D#####
##making the dataset for upset of UpSetR (binary matrix, present or absent)
data2_ARG.water.group <- tax_glom (data2_ARG.water, taxrank = "Group", NArm = F)
data2_ARG.water.group ##groups are present in water samples (12 samples)
upset.data2_ARG.water.group <- MicrobiotaProcess::get_upset(data2_ARG.water.group, factorNames="gen_material") 
upset.data2_ARG.water.group

upset.data2_ARG.water.group <- upset.data2_ARG.water.group %>%
  rename("RNA(cDNA)" = cDNA)

#Plot
upset_plot_ARG_CB <- UpSetR::upset(upset.data2_ARG.water.group,
              sets.bar.color = c("#CC79A7", "#009E73"), 
              order.by = "freq", text.scale = c(3, 2.5, 2.5, 2, 3, 2), 
              point.size = 5, line.size = 2, mainbar.y.label= "Group count",
              sets.x.label = "Group count", 
              set_size.show = F) 
upset_plot_ARG_CB

##saving the upset plot - Figure 8D
svglite("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/Figure8D.svg",
        width=7, height=6)
upset_plot_ARG_CB
dev.off()
