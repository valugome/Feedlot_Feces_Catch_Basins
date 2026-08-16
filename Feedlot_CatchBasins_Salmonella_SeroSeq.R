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
  "tidyverse", "dplyr", "stringr", "paletteer", "vegan", "cowplot",
  "ggplot2",
  "ggdendro", "randomcoloR", "ggpubr", "ggsignif", 
  "ggtext", "ggnewscale", "rstatix", "ggrepel", "ggh4x",
  "writexl", 
  "Polychrome", "colorspace"
)

load_packages(cran_pkgs)

#Load data####
##SEROVARS PROPORTION#####
seroseq_serovars <- read_csv('Data/SeroSeq_Data.csv')

##Pivot to longer format 
seroseq_serovars_long <- pivot_longer(seroseq_serovars,
                                      !c("Sample", "SampleID",
                                         "Matrix", "Feedlot", 
                                         "WL Reads", "invA", "Notes",
                                         "CSS Reads"),
                                      names_to = "Serovar",
                                      values_to = "Proportion")%>%
  ##Fixing matrix names, and feedlot names
  mutate(Matrix = case_when(Matrix == "feces" ~ "FECES",
                            Matrix == "water" ~ "CB"))


### Set labels for x axis (instead of SampleID, want matrix)
xlabels_matrix <- setNames(seroseq_serovars_long$Matrix, seroseq_serovars_long$SampleID)
xlabels_matrix ##OK now!

##Color palette
serovars_palette <- as.character(paletteer_d("colorBlindness::paletteMartin"))
serovars_palette <- serovars_palette[-1] ##Don't want black (first color)

##PLOT SEROSEQ - SEROVAR DATA
#Plot
seroseq_serovars.plot <- ggplot(seroseq_serovars_long, 
                                aes(x=SampleID, y= Proportion, fill = Serovar)) +
  theme_minimal() +
  facet_grid(~Feedlot, 
             scales = "free", 
             space = "free_x")+
  labs(y= "Relative Proportion",
       #title = "Salmonella Serovars per Feedlot"
       ) +
  geom_bar(stat = "summary", colour = "black", na.rm = TRUE)+
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0), labels = xlabels_matrix) +
  scale_fill_manual(values = serovars_palette)+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 25),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black", 
                                   angle = 90, 
                                   vjust = 0.5,
                                   hjust = 0.5),
        strip.text=element_text(size=18, color = "black", face = "bold"),
        # strip.background=element_rect(fill='black'
        #                               # , color='white'),
        strip.placement="outside", panel.spacing.y=unit(0.5, "lines"),
        panel.spacing = unit(0.8, "lines"), 
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = 1.5))
seroseq_serovars.plot

##SEROVARS PER SAMPLE TYPE#####
# Summarize number of serovars per sample
serovars_per_sample_df <- seroseq_serovars_long %>%
  filter(Proportion > 0) %>%  # Only proportions > 0
  group_by(SampleID, Matrix) %>%# Keep Matrix info per sample
  summarise(serovars_per_sample = n_distinct(Serovar), .groups = "drop") # Count unique serovars
serovars_per_sample_df

##AVERAGE NUMBER OF SEROVARS PER SAMPLE TYPE 
serovars_per_sample_df %>%
  group_by(Matrix)%>%
  summarise(mean_serovars_per_sample_type = mean(serovars_per_sample), 
            sd_serovars_per_sample = sd(serovars_per_sample)) #Feces = 2.16, CB 1   

# Count how many samples have each number of serovars by Matrix
serovar_counts <- serovars_per_sample_df %>%
  group_by(serovars_per_sample, Matrix) %>%
  summarise(num_samples = n(), .groups = "drop")

serovar_counts$serovars_per_sample <- as.factor(serovar_counts$serovars_per_sample) # Convert to factor for discrete x-axis

# Plot the stacked bar chart
serovars_per_sample <- ggplot(serovar_counts, aes(x = serovars_per_sample, y = num_samples, fill = Matrix)) +
  geom_bar(stat = "identity") +
  labs(x = "No. of serovars per sample", y = "No. of Samples", fill = "Sample Type") +
  scale_y_continuous(breaks = seq(1, 12, by = 2), expand = c(0, 0, 0, 1))+
  scale_fill_manual(values = c("#5D729D", "brown"),
                    labels = c("FECES" = "FECES", "CB" = "CATCH BASINS")) +
  theme_minimal() +
  theme(legend.position  = "inside",
        legend.position.inside = c(0.8, 0.9),
        legend.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        panel.border = element_rect(colour = "black", linewidth= 1, fill = NA),
        axis.title = element_text(size = 20, colour = "black"),
        axis.text = element_text(size = 16, colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.5))
serovars_per_sample

##FIGURE 1 ########
figure1 <- plot_grid(seroseq_serovars.plot,
                     serovars_per_sample, 
                     labels = "AUTO",
                     label_size = 22,
                     ncol = 1,
                     rel_heights = c(0.75, 0.25))
figure1



##BETADIV OF SALMONELLA SEROVARS####
#OTU table###
##Pivot to wide format for OTU-like table
seroseq_serovars_wide <- seroseq_serovars_long %>%
  pivot_wider(
    id_cols = Serovar,
    names_from = SampleID,
    values_from = Proportion
  )%>%
  column_to_rownames(var="Serovar")%>% ##OTU IDs
  replace(is.na(.), 0) #replace NA with 0


any(colSums(seroseq_serovars_wide)== 0) ## no samples with 0 OTUs

#Make into phyloseq
##Metadata####
metadata <- read.csv('Data/Metadata_Feedlot_CatchBasins.csv', 
                     check.names = F,
                     row.names = "sampleID")

metadata$SampleID<- rownames(metadata) #Making a SampleID column 

##PHYLOSEQ
OTU_seroseq <-phyloseq::otu_table(seroseq_serovars_wide, taxa_are_rows = TRUE)
phyloseq_seroseq <- phyloseq(OTU_seroseq, sample_data(metadata))

##Calculate jaccards
seroseq.jac <- vegdist(t(phyloseq_seroseq@otu_table), method = "jaccard") 

#Metadata for seroseq samples only
seroseq.sampledata.df <- data.frame(phyloseq_seroseq@sam_data)%>%
  mutate(feedlot = factor(feedlot, levels = c("1", "2", "3", "4", "5")))

#PERMANOVA
set.seed(98)
seroseq_adonis  <- adonis2(seroseq.jac ~ sample_type + feedlot, 
                           seroseq.sampledata.df, 
                           by = "margin", 
                           permutations = 9999)
seroseq_adonis  #30.8% of variation is due to feedlot p = 6e-04

#Interaction - NS
#PERMANOVA
set.seed(98)
seroseq_adonis_interaction  <- adonis2(seroseq.jac ~ sample_type * feedlot, 
                           seroseq.sampledata.df, 
                           by = "margin", 
                           permutations = 9999)
seroseq_adonis_interaction  #2.8% of variation is due to the interaction between sample type and feedlot p = 0.9837

######Supp table 5.0#######
stable5.0 <- data.frame(seroseq_adonis, check.names = F)%>%
  mutate(Dataset = "Salmonella Seroseq Serovars")%>%
  rownames_to_column(var = "Fixed Effect")%>%
  mutate(`Fixed Effect` = str_replace(`Fixed Effect`, "sample_type", "Sample Type"))%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))
stable5.0

###SUPPLEMENTARY TABLE 5.0 SALMONELLA SEROSEQ SECTION#######
stable5_salmonella_seroseq <- stable5.0%>%
  mutate(SumOfSqs = round(SumOfSqs, 2), 
         `R2` = round(`R2`,3), 
         `F` = round(`F`,2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3))%>%
  select(Dataset, `Fixed Effect`, `Df`, SumOfSqs, R2 , `F`, `Pr(>F)`)
stable5_salmonella_seroseq


##PERMDISP
#Run betadisp average distance to centroid
betadisp_seroseq_feedlot <- betadisper(seroseq.jac, seroseq.sampledata.df$feedlot)
##Then test by permuting
set.seed(98)
seroseq_feedlot.permdisp <- permutest(betadisp_seroseq_feedlot, permutations = 9999)
seroseq_feedlot.permdisp ##NS, p = 0.0682

