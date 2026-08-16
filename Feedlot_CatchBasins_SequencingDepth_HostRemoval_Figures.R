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
  "tidyverse", "dplyr", "stringr", "paletteer", "vegan", "cowplot",
  "ggplot2",
  "ggdendro", "randomcoloR", "ggpubr", "ggsignif", 
  "ggtext", "ggnewscale", "rstatix", "ggrepel", "ggh4x",
  "writexl", 
  "Polychrome", "colorspace"
)

load_packages(cran_pkgs)


#IMPORT METADATA####
metadata <- read.csv('Data/Metadata_Feedlot_CatchBasins.csv', 
                     check.names = F)%>%
  mutate(SampleID=sampleID)%>%
  column_to_rownames(var="sampleID") #making sampleID rownames

##Host-free reads####
hostrem <- read.csv('Data/HostRem_Reads_Feedlot_CatchBasins.csv')

#Change zymo- to zymo. as they are in metadata
hostrem[128, "SampleID"] <- "Zymo.1a"
hostrem[129, "SampleID"] <- "Zymo.1b"

#Merge with metadata
hostrem <- hostrem %>%
  filter(!(grepl("F5W02", SampleID)))%>% #Not using F5W02 or F5W02c
  select(SampleID, Hostrem_output_total_num_seqs, Hostrem_percentage_total_seqs_removed)%>%
  left_join(metadata, by = "SampleID")%>%
  rename(HostFree_Reads=Hostrem_output_total_num_seqs)#Hostrem_output_total_num_seqs is the total number of host free reads

#Color Palettes
feedlot_palette <- c("1" = "#fcca46", 
                     "2" = "#fe7f2d", 
                     "3" = "#233d4d", 
                     "4"= "#3b9ab2", 
                     "5"= "#e1b6ff")
###Sample type
sample.type.palette <- c("Water" = "#4C72B0",
                         "Feces" = "brown") 
#Library type
gen.material.palette <- c("cDNA" = "#009E73",  
                          "DNA"  = "#CC79A7" )  

#sequencing batches
batch_palette <- c("no" = "#d19bac", 
                   "yes" = "#6a9c55") #sequencing batches

##Salmonella color palette
salmonella.palette <- c("positive"= "#fc8d62", "negative" = "#8da0cb")



#SEQUENCING DEPTH #####
##Trimmomatic stats###
trimmomatic_reads <- read.csv('Data/Trimmomatic_Reads_Feedlot_CatchBasins.csv')


#Merge with metadata
seqdepth <- trimmomatic_reads %>%
  select(SampleID, Trimm_NumberOfInputReads_Paired, Resequenced)%>%
  left_join(metadata, by = "SampleID")%>%
  mutate(SequencingDepth_Reads=Trimm_NumberOfInputReads_Paired*2)#Trimmomatic input is the total amount of reads sequenced

##Check out sequencing depths, take out those with low reads
sort(seqdepth$SequencingDepth_Reads) #Low ones belong to blanks, NTCs and F5W02 and F5W02c

#Only samples and zymos 
seqdepth.samples_zymos <- seqdepth%>%
  filter(!grepl("NTC|EB|F5W02", SampleID))
nrow(seqdepth.samples_zymos) #Left with 122 samples

#Only samples
seqdepth.samples <- seqdepth%>%
  filter(!grepl("NTC|EB|Zymo|F5W02", SampleID))
nrow(seqdepth.samples) #Left with 120 samples

##Factor order 
seqdepth.samples <- seqdepth.samples%>%
  mutate(sample_type = factor (sample_type, levels = c("Feces", "Water")))
seqdepth.samples <- seqdepth.samples%>%
  mutate(gen_material = factor (gen_material, levels = c("DNA", "cDNA")))

######cDNA vs DNA faceted by CB and Feces- N.S.####
sequencing_depth_cDNAvDNA.WandF<- ggplot(seqdepth.samples, aes(x = gen_material, y= SequencingDepth_Reads, color = gen_material, fill = gen_material)) +
  theme_bw() +
  labs(y= "Reads", color = "Library Type", fill = "Library Type") +
  facet_grid(~sample_type, scales = "free",  labeller = as_labeller(c("Feces" = "FECES", "Water" = "CATCH BASIN"))) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  scale_fill_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA(cDNA)"))+
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA(cDNA)"))+
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
sequencing_depth_cDNAvDNA.WandF

#Checking Stats
water.seqdepth.samples <- subset(seqdepth.samples, sample_type == "Water")
wilcox_test(water.seqdepth.samples%>%arrange(SampleID), #Have to make sure they are arranged by sampleID so 'paired' works correclty
            SequencingDepth_Reads~gen_material,
            paired = TRUE) #n.s. p = 0.519
#cDNA and DNA n = 12
feces.seqdepth.samples <- subset(seqdepth.samples, sample_type == "Feces")
wilcox_test(feces.seqdepth.samples%>%arrange(SampleID), 
            SequencingDepth_Reads~gen_material,
            paired = TRUE) #n.s. p = 0.137
#cDNA and DNA n = 48

######CB vs Feces faceted by cDNA and DNA - Not significant#####
sequencing_depth_WvF.cDNAandDNA<- ggplot(seqdepth.samples, aes(x = sample_type, y= SequencingDepth_Reads, color = sample_type, fill = sample_type)) +
  theme_bw() +
  labs(y= "Reads", color = "Sample Type", fill = "Sample Type") +
  facet_wrap(~gen_material, labeller = as_labeller(c("DNA" = "DNA", "cDNA" = "RNA(cDNA)"))) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +  scale_fill_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch basin"))+
  scale_color_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch basin"))+
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
# geom_pwc(method = "wilcox_test",
#          label = "p = {p.adj}",
#          hide.ns = TRUE,
#          step.increase = 0.08,
#          label.size = 6,
#          tip.length = 0.02)
sequencing_depth_WvF.cDNAandDNA


#Checking Stats
cDNA.seqdepth.samples <- subset(seqdepth.samples, gen_material == "cDNA")
wilcox_test(cDNA.seqdepth.samples, SequencingDepth_Reads~ sample_type) #ns, p = 0.0875
#Feces n=48,Water  n=12 
DNA.seqdepth.samples <- subset(seqdepth.samples, gen_material == "DNA")
wilcox_test(DNA.seqdepth.samples, SequencingDepth_Reads~ sample_type) #ns p = 0.82

#######Feedlot - Not significant ####
sequencing_depth_feedlot <- ggplot(seqdepth.samples, aes(x = feedlot, y= SequencingDepth_Reads, color = factor(feedlot), fill = factor(feedlot))) +
  theme_bw() +
  labs(title = "SEQUENCING DEPTH", y= "Reads", color = "Feedlot", fill = "Feedlot" ) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +  scale_fill_manual(values = feedlot_palette) +
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
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) 
# geom_pwc(method = "wilcox_test", p.adjust.method = "BH",label = "p = {p.adj}",
#          hide.ns = TRUE,
#          step.increase = 0.08,
#          label.size = 5,
#          tip.length = 0.02)
sequencing_depth_feedlot

#Checking Stats
wilcox_test(seqdepth.samples, 
            SequencingDepth_Reads~feedlot, 
            p.adjust.method = "BH") #Not significant

#######Checking for batch effects - reseq vs first run - Not significant ####
sequencing_depth_batches <- ggplot(seqdepth.samples, aes(x = re_sequenced, y= SequencingDepth_Reads, color = re_sequenced, fill = re_sequenced)) +
  theme_bw() +
  labs(title = "SEQUENCING DEPTH", y= "Reads",color = "Re-Sequenced", fill = "Re-Sequenced") +
  geom_boxplot(alpha = 0.1) +
  #geom_boxplot(aes(group = factor(feedlot)), alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
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
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))
sequencing_depth_batches

#Checking Stats
wilcox_test(seqdepth.samples%>%arrange(SampleID), 
            SequencingDepth_Reads~re_sequenced) #NS (p = 0.814)
#Not resequenced n=114, resequenced n=6  

######SUPPLEMENTARY FIGURE 1#####
sfigure1 <- cowplot::plot_grid(
  sequencing_depth_cDNAvDNA.WandF + theme(axis.title.y = element_text(size = 20)), 
  sequencing_depth_WvF.cDNAandDNA + theme(axis.title.y = element_text(size = 20)), 
  sequencing_depth_batches + theme(plot.title = element_blank(),
                                   axis.title.y = element_text(size = 20)),
  labels = "AUTO", 
  label_size = 22,
  ncol = 1)+
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
sfigure1
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryFigure1.png", 
       plot = sfigure1, 
       device = "png", width = 10, height = 10, dpi = 600)

#HOST REMOVAL STATS#####
hostrem

#Only samples and zymos (also taking out F5W02 and F5W02c cuz of low sequencing depth)
hostrem.samples_zymos <- hostrem%>%
  filter(!grepl("NTC|EB|F5W02", SampleID))
nrow(hostrem.samples_zymos) #Left with 122 samples

#Only analyzing samples (also taking out F5W02 and F5W02c cuz of low sequencing depth)
hostrem.samples <- hostrem%>%
  filter(!grepl("NTC|EB|Zymo|F5W02", SampleID))

#Factor order
hostrem.samples <- hostrem.samples%>%
  mutate(gen_material = factor (gen_material, levels = c("DNA", "cDNA")))

###cDNA vs DNA faceted by CB and Feces####
####Reads - S. for feces #####
host_free_reads_cDNAvDNA.WandF<- ggplot(hostrem.samples%>%arrange(SampleID), 
                                        aes(x = gen_material, y= HostFree_Reads, color = gen_material, fill = gen_material)) +
  theme_bw() +
  labs(y= "Reads", color = "Library Type", fill = "Library Type") +
  facet_grid(~sample_type, scales = "free",  labeller = as_labeller(c("Feces" = "FECES", "Water" = "CATCH BASINS"))) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  scale_fill_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA(cDNA)"))+
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA(cDNA)"))+
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
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc(method = "wilcox_test",
           label = "Wilcoxon, p = {p}",
           method.args = list(paired = TRUE),
           hide.ns = TRUE,
           step.increase = 0.08,
           label.size = 6,
           tip.length = 0.02)
host_free_reads_cDNAvDNA.WandF

#Checking Stats
water.hostrem.samples <- subset(hostrem.samples, sample_type == "Water")
wilcox_test(water.hostrem.samples%>%arrange(SampleID), 
            HostFree_Reads~gen_material, 
            paired = T) #n.s. p = 0.85
#cDNA and DNA n = 12

feces.hostrem.samples <- subset(hostrem.samples, sample_type == "Feces")
wilcox_test(feces.hostrem.samples%>%arrange(SampleID),
            HostFree_Reads~gen_material,
            paired = T) #S. p = 0.0221
#cDNA and DNA n = 48

####Percentage - SIGNIFICANT####
host_free_percentage_cDNAvDNA.WandF<- ggplot(hostrem.samples%>%arrange(SampleID), 
                                             aes(x = gen_material, y= Hostrem_percentage_total_seqs_removed, 
                                                                  color = gen_material, fill = gen_material)) +
  theme_bw() +
  labs(y= "Percentage (%)\n of reads", color = "Library Type", fill = "Library Type") +
  facet_grid(~sample_type, scales = "free",  labeller = as_labeller(c("Feces" = "FECES", "Water" = "CATCH BASINS"))) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  scale_fill_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA(cDNA)"))+
  scale_color_manual(values = gen.material.palette, labels = c("DNA" = "DNA", "cDNA" = "RNA(cDNA)"))+
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
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc(method = "wilcox_test",
           label = "Wilcoxon, p = {p}",
           hide.ns = TRUE,
           method.args = list(paired = TRUE),
           step.increase = 0.08,
           label.size = 6,
           tip.length = 0.02)
host_free_percentage_cDNAvDNA.WandF

#Checking Stats
#CB
wilcox_test(water.hostrem.samples%>%arrange(SampleID), 
            Hostrem_percentage_total_seqs_removed~gen_material, 
            paired = T) #s. p = 0.000977
#cDNA and DNA n = 12
#Feces
wilcox_test(feces.hostrem.samples%>%arrange(SampleID), 
            Hostrem_percentage_total_seqs_removed~gen_material,
            paired = T) #s. p = 7.11e-15
#cDNA and DNA n = 48

###CB vs Feces faceted by cDNA and DNA#####
####Reads -N.S.#####
host_free_reads_WvF.cDNAandDNA<- ggplot(hostrem.samples, aes(x = sample_type, y= HostFree_Reads, color = sample_type, fill = sample_type)) +
  theme_bw() +
  labs(y= "Reads", color = "Sample Type", fill = "Sample Type") +
  facet_wrap(~gen_material, labeller = as_labeller(c("DNA" = "DNA", "cDNA" = "RNA(cDNA)"))) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2)  +
  scale_fill_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch basin"))+
  scale_color_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch basin"))+
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
# geom_pwc(method = "wilcox_test",
#          label = "p = {p.adj}",
#          hide.ns = TRUE,
#          step.increase = 0.08,
#          label.size = 6,
#          tip.length = 0.02)
host_free_reads_WvF.cDNAandDNA


#Checking Stats
cDNA.hostrem.samples <- subset(hostrem.samples, gen_material == "cDNA")
wilcox_test(cDNA.hostrem.samples, HostFree_Reads~ sample_type) #ns, p = 0.237
#Feces n=48,Water  n=12 
DNA.hostrem.samples <- subset(hostrem.samples, gen_material == "DNA")
wilcox_test(DNA.hostrem.samples, HostFree_Reads~ sample_type) #ns p = 0.82


####Percentage - SIGNIFICANT FOR DNA #####
host_free_percentage_WvF.cDNAandDNA<- ggplot(hostrem.samples, aes(x = sample_type, 
                                                                  y= Hostrem_percentage_total_seqs_removed,
                                                                  color = sample_type, fill = sample_type)) +
  theme_bw() +
  labs(y= "Percentage (%)\n of reads", color = "Sample Type", fill = "Sample Type") +
  facet_wrap(~gen_material, labeller = as_labeller(c("DNA" = "DNA", "cDNA" = "RNA(cDNA)"))) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  scale_fill_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch basin"))+
  scale_color_manual(values = sample.type.palette, labels = c("Feces" = "Feces", "Water" = "Catch basin"))+
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
host_free_percentage_WvF.cDNAandDNA


#Checking Stats
cDNA.hostrem.samples <- subset(hostrem.samples, gen_material == "cDNA")
wilcox_test(cDNA.hostrem.samples, Hostrem_percentage_total_seqs_removed~ sample_type) #ns, p = 0.935
#Feces n=48,Water  n=12 
DNA.hostrem.samples <- subset(hostrem.samples, gen_material == "DNA")
wilcox_test(DNA.hostrem.samples, Hostrem_percentage_total_seqs_removed~ sample_type) #S p = 0.0044


#####SUPPLEMENTARY FIGURE 2##########
sfigure2 <- cowplot::plot_grid(
  host_free_percentage_WvF.cDNAandDNA + theme(axis.title.y = element_text(size = 20)),
  host_free_percentage_cDNAvDNA.WandF + theme(axis.title.y = element_text(size = 20)), 
  ncol =1,
  labels = "AUTO", label_size = 22)+
  theme(
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
sfigure2 
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryFigure2.png", 
       plot = sfigure2,  
       device = "png", 
       dpi = 600, 
       width = 10, height =10)


###Feedlot - Not significant ####
host_free_reads_feedlot <- ggplot(hostrem.samples, aes(x = feedlot, y= HostFree_Reads, color = factor(feedlot), fill = factor(feedlot))) +
  theme_bw() +
  labs(title = "HOST-FREE READS ", y= "Reads", color = "Feedlot", fill = "Feedlot" ) +
  geom_boxplot(alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
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
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))
# geom_pwc(method = "wilcox_test", p.adjust.method = "BH",label = "p = {p.adj}",
#          hide.ns = TRUE,
#          step.increase = 0.08,
#          label.size = 5,
#          tip.length = 0.02)
host_free_reads_feedlot

#Checking Stats
wilcox_test(hostrem.samples, 
            HostFree_Reads~feedlot, 
            p.adjust.method = "BH") #Not significant

###Checking for batch effects - reseq vs first run - Not significant ####
host_free_reads_batches <- ggplot(hostrem.samples, aes(x = re_sequenced, y= HostFree_Reads, color = re_sequenced, fill = re_sequenced)) +
  theme_bw() +
  labs(title = "HOST-FREE READS ", y= "Reads",color = "Re-Sequenced", fill = "Re-Sequenced") +
  geom_boxplot(alpha = 0.1) +
  #geom_boxplot(aes(group = factor(feedlot)), alpha = 0.1) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
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
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))
host_free_reads_batches

#Checking Stats
wilcox_test(hostrem.samples, HostFree_Reads~re_sequenced) #NS (p = 0.271)
#Not resequenced n=114, resequenced  n=6  
