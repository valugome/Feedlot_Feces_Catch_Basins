#LOAD R PACKAGES ######
setwd('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Analysis_repository/Microbiome/Bacteria_archaea/R_Analysis')
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
  "tidyverse", "dplyr", "stringr", "paletteer", "vegan"
)

load_packages(cran_pkgs)

#Load data####
##Trimmomatic stats#####
trimmomatic <- read.csv('Data/Trimmomatic_Reads_Feedlot_CatchBasins.csv')

##SortMeRNA stats#####
sortmerna <- read.csv('Data/SortMeRNA_Reads_Feedlot_CatchBasins.csv')

##FLASH stats#####
flash <- read.csv('Data/FLASH_Reads_Feedlot_CatchBasins.csv')

##Host Removal stats#####
hostrem <- read.csv('Data/HostRem_Reads_Feedlot_CatchBasins.csv')

#Merge dataframes
final_feedlot_lagoon_stats <- left_join(trimmomatic, 
                                        sortmerna, 
                                        by = "SampleID")%>%
  left_join(flash, by = "SampleID")%>%
  left_join(hostrem, by = "SampleID")%>%
  filter(!grepl("F5W02", SampleID)) ##Not including F5W02
nrow(final_feedlot_lagoon_stats) #129 = 120 samples, 4 extraction blanks (2 EBs, 2 EBc), 3 NTC, 2 Zymo mock communities

#SEQUENCING EFFORT#####
#How many reads did we get from sequencing? Samples
final_feedlot_lagoon_stats%>%
  filter(!grepl("EB|NTC|Zymo", SampleID))%>%
  mutate(Trimm_NumberOfInputReads = (Trimm_NumberOfInputReads_Paired * 2))%>%
  mutate(gen_material = ifelse(grepl("c", SampleID), "RNA", "DNA"))%>%
  group_by(gen_material)%>%
  summarise(min_sequences_per_sample = min(Trimm_NumberOfInputReads),
            max_sequences_per_sample = max(Trimm_NumberOfInputReads),
            mean_sequences_per_sample = mean(Trimm_NumberOfInputReads),
            median_sequences_per_sample = median(Trimm_NumberOfInputReads))
#gen_material min_sequences_per_sample max_sequences_per_sample mean_sequences_per_sample median_sequences_per_sample
#DNA          56642584                115179420                 81702921.                   81437535
#RNA          42285580                104109092                 78857885.                   80742616

#How many reads did we get from sequencing? Zymo mock
final_feedlot_lagoon_stats%>%
  filter(grepl("Zymo", SampleID))%>%
  mutate(Trimm_NumberOfInputReads = (Trimm_NumberOfInputReads_Paired * 2))%>%
  summarise(min_sequences_per_sample = min(Trimm_NumberOfInputReads),
            max_sequences_per_sample = max(Trimm_NumberOfInputReads),
            mean_sequences_per_sample = mean(Trimm_NumberOfInputReads),
            median_sequences_per_sample = median(Trimm_NumberOfInputReads))
#min_sequences_per_sample max_sequences_per_sample mean_sequences_per_sample median_sequences_per_sample
#92868846                 95956714                  94412780                    94412780

#How many reads did we get from sequencing? Blanks and non template controls
final_feedlot_lagoon_stats%>%
  filter(grepl("NTC|EB", SampleID))%>%
  mutate(Trimm_NumberOfInputReads = (Trimm_NumberOfInputReads_Paired * 2))%>%
  summarise(min_sequences_per_sample = min(Trimm_NumberOfInputReads),
            max_sequences_per_sample = max(Trimm_NumberOfInputReads),
            mean_sequences_per_sample = mean(Trimm_NumberOfInputReads),
            median_sequences_per_sample = median(Trimm_NumberOfInputReads)) 
#Mean per sample (7 samples - 4 EBs and 3 NTC) 5441.714, median per sample 510. 

#QC SEQUENCES - TRIMMOMATIC ####
#What percentage of reads did we filter with trimmomatic? Samples
final_feedlot_lagoon_stats%>%
  filter(!grepl("EB|NTC|Zymo", SampleID))%>%
  pull(Trimm_percentage_pairedreads_cut, SampleID)%>%
  summary()
 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.954   2.653   3.196   3.141   3.51   4.932 

#What percentage of reads did we filter with trimmomatic? Zymo mock communities
final_feedlot_lagoon_stats%>%
  filter(grepl("Zymo", SampleID))%>%
  pull(Trimm_percentage_pairedreads_cut, SampleID)%>%
  summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.586   2.775   2.965   2.965   3.154   3.343

#What percentage of reads did we filter with trimmomatic? Blanks and non template controls
final_feedlot_lagoon_stats%>%
  filter(grepl("NTC|EB", SampleID))%>%
  pull(Trimm_percentage_pairedreads_cut, SampleID)%>%
  summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#3.252  77.090  92.613  78.002  97.983 100.000 

#SortMeRNA rRNA removal####
##SAMPLES####
#What percentage of reads did we filter with sortmerna? 
final_feedlot_lagoon_stats%>%
  filter(!grepl("EB|NTC|Zymo", SampleID))%>%
  pull(SortMeRNA_percentage_removed_reads, SampleID)%>%
  summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.4156  0.9433  1.0195  1.8922  1.0621 40.2354      60 

##BLANKS and NTC####
#What percentage of reads did we filter with sortmerna? 
final_feedlot_lagoon_stats%>%
  filter(grepl("NTC|EB", SampleID))%>%
  pull(SortMeRNA_percentage_removed_reads, SampleID)
#EB2      NTC      EBc     EBc2       EB     NTC1     NTC2 
#NA       NA 1.789345 7.009346       NA       NA       NA 

#FLASH MERGING####
##PERCENTAGE MERGED READS#####
summary(final_feedlot_lagoon_stats$FLASH_percentage_merged_reads)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    14.9    20.4    26.4    26.8    31.5    57.6       1 

###SAMPLES#####
final_feedlot_lagoon_stats %>%
  filter(!grepl("EB|NTC|Zymo", SampleID))%>%
  pull(FLASH_percentage_merged_reads)%>%
  summary
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14.9    20.3    26.0    26.1    30.9    44.8 

###ZYMOS#####
final_feedlot_lagoon_stats %>%
  filter(grepl("Zymo", SampleID))%>%
  pull(FLASH_percentage_merged_reads)%>%
  summary
# Min.   1st Qu.  Median.  Mean   3rd Qu.  Max. 
# 27.21   28.73   30.25   30.25   31.77   33.29 


#HOST REMOVAL####
##SAMPLES####
###MERGED READS#####
final_feedlot_lagoon_stats%>%
  filter(!grepl("EB|NTC|Zymo", SampleID))%>%
  mutate(sample_type=ifelse(grepl("W", SampleID), "Water", "Feces"))%>%
  group_by(sample_type)%>%
  summarise(mean_percentage_host_mapping_seqs = mean(Hostrem_percentage_merged_seqs_removed),
          min_percentage_host_mapping_seqs = min(Hostrem_percentage_merged_seqs_removed),
          max_percentage_host_mapping_seqs = max(Hostrem_percentage_merged_seqs_removed))

# sample_type mean_percentage_host_mapping_seqs min_percentage_host_mapping_seqs max_percentage_host_mapping_seqs
# <chr>                                   <dbl>                            <dbl>                            <dbl>
#  Feces                                   0.135                           0.0172                             1.39
# Water                                   0.124                           0.0185                             1.01

###UNMERGED READS#####
final_feedlot_lagoon_stats%>%
  filter(!grepl("EB|NTC|Zymo", SampleID))%>%
  mutate(sample_type=ifelse(grepl("W", SampleID), "Water", "Feces"))%>%
  group_by(sample_type)%>%
  summarise(mean_percentage_host_mapping_seqs = mean(Hostrem_percentage_unmerged_seqs_removed),
          min_percentage_host_mapping_seqs = min(Hostrem_percentage_unmerged_seqs_removed),
          max_percentage_host_mapping_seqs = max(Hostrem_percentage_unmerged_seqs_removed))

# sample_type mean_percentage_host_mapping_seqs min_percentage_host_mapping_seqs max_percentage_host_mapping_seqs
# <chr>                                   <dbl>                            <dbl>                            <dbl>
#  Feces                                   0.133                           0.0582                            0.672
#  Water                                   0.148                           0.0533                            0.637


###TOTAL READS#####
final_feedlot_lagoon_stats%>%
  filter(!grepl("EB|NTC|Zymo", SampleID))%>%
  mutate(sample_type=ifelse(grepl("W", SampleID), "Water", "Feces"),
         gen_material=ifelse(grepl("c", SampleID), "cDNA", "DNA"))%>%
  group_by(sample_type, gen_material)%>%
  summarise(mean_percentage_host_mapping_seqs = mean(Hostrem_percentage_total_seqs_removed),
            min_percentage_host_mapping_seqs = min(Hostrem_percentage_total_seqs_removed),
            max_percentage_host_mapping_seqs = max(Hostrem_percentage_total_seqs_removed))
# sample_type gen_material  mean_percentage_host_mapping_seqs    min_percentage_host_mapping_seqs  max_percentage_host_mapping_…¹
#   1 Feces       DNA                                   0.187                            0.121                           0.763
# 2 Feces       cDNA                                    0.0781                           0.0542                          0.222
# 3 Water       DNA                                     0.180                            0.0837                          0.685
# 4 Water       cDNA                                    0.105                            0.0495                          0.336

##ZYMOS####
###MERGED READS#####
final_feedlot_lagoon_stats%>%
  filter(grepl("Zymo", SampleID))%>%
  summarise(mean_percentage_host_mapping_seqs = mean(Hostrem_percentage_merged_seqs_removed),
          min_percentage_host_mapping_seqs = min(Hostrem_percentage_merged_seqs_removed),
          max_percentage_host_mapping_seqs = max(Hostrem_percentage_merged_seqs_removed))
# mean_percentage_host_mapping_seqs min_percentage_host_mapping_seqs max_percentage_host_mapping_seqs
#                         0.2644885                          0.26219                        0.2667869

###UNMERGED READS#####
final_feedlot_lagoon_stats%>%
  filter(grepl("Zymo", SampleID))%>%
  summarise(mean_percentage_host_mapping_seqs = mean(Hostrem_percentage_unmerged_seqs_removed),
          min_percentage_host_mapping_seqs = min(Hostrem_percentage_unmerged_seqs_removed),
          max_percentage_host_mapping_seqs = max(Hostrem_percentage_unmerged_seqs_removed))
# mean_percentage_host_mapping_seqs min_percentage_host_mapping_seqs max_percentage_host_mapping_seqs
#                         0.3305684                        0.3234374                        0.3376994

###TOTAL READS#####
final_feedlot_lagoon_stats%>%
  filter(grepl("Zymo", SampleID))%>%
  summarise(mean_percentage_host_mapping_seqs = mean(Hostrem_percentage_total_seqs_removed),
            min_percentage_host_mapping_seqs = min(Hostrem_percentage_total_seqs_removed),
            max_percentage_host_mapping_seqs = max(Hostrem_percentage_total_seqs_removed))
# mean_percentage_host_mapping_seqs min_percentage_host_mapping_seqs max_percentage_host_mapping_seqs
#                        0.3185669                        0.3145155                        0.3226183

##BLANKS#####
###MERGED READS#####
final_feedlot_lagoon_stats%>%
  filter(grepl("EB|NTC", SampleID))%>%
  summarise(mean_percentage_host_mapping_seqs = mean(Hostrem_percentage_merged_seqs_removed),
          min_percentage_host_mapping_seqs = min(Hostrem_percentage_merged_seqs_removed),
          max_percentage_host_mapping_seqs = max(Hostrem_percentage_merged_seqs_removed))

###UNMERGED READS#####
final_feedlot_lagoon_stats%>%
  filter(grepl("EB|NTC", SampleID))%>%
  summarise(mean_percentage_host_mapping_seqs = mean(Hostrem_percentage_unmerged_seqs_removed),
          min_percentage_host_mapping_seqs = min(Hostrem_percentage_unmerged_seqs_removed),
          max_percentage_host_mapping_seqs = max(Hostrem_percentage_unmerged_seqs_removed))

###TOTAL READS#####
final_feedlot_lagoon_stats%>%
  filter(grepl("EB|NTC", SampleID))%>%
  summarise(mean_percentage_host_mapping_seqs = mean(Hostrem_percentage_total_seqs_removed),
            min_percentage_host_mapping_seqs = min(Hostrem_percentage_total_seqs_removed),
            max_percentage_host_mapping_seqs = max(Hostrem_percentage_total_seqs_removed))


##SUPPLEMENTARY TABLE 2#######
stable2 <- final_feedlot_lagoon_stats %>%
  mutate(`Library Type` = ifelse(grepl("c", SampleID), "Metatranscriptomic (RNA (cDNA))", "Metagenomic (DNA)"),
         `Sample Type` = ifelse(grepl("W", SampleID), "Catch Basin", "Feces"))%>%
  arrange(SampleID)%>%
  mutate(`Trimmmomatic Number of Input Reads (Paired)` = Trimm_NumberOfInputReads_Paired, 
         `Trimmmomatic Percentage of Reads Dropped` = round(Trimm_percentage_pairedreads_cut,2),
         `Trimmmomatic Number of Output Reads (Paired)` = Trimm_total_pairedreads_left,
         `SortMeRNA Number of Input Reads (Total)` = SortMeRNA_input_number_reads,
         `SortMeRNA Number of Output Reads (Total)` = SortmeRNA_post_number_reads,            
         `SortMeRNA Percentage of Reads Removed` = round(SortMeRNA_percentage_removed_reads,2),    
         `FLASH Number of Input Reads (Total)` = FLASH_input_number_reads,
         `FLASH Number of Reads Merged (Total)` = FLASH_merged_num_seqs,
         `FLASH Percentage of Reads Merged (%)` = round(FLASH_percentage_merged_reads,2), 
         `Min Length of Merged Reads` = FLASH_merged_min_len,
         `Max Length of Merged Reads` = FLASH_merged_max_len,
         `Average Length of Merged Reads` = round(FLASH_merged_avg_len,2), 
         `FLASH Number of Reads Unmerged (Total)` = FLASH_unmerged_num_seqs,
         `FLASH Percentage of Reads Unmerged (%)` = round(FLASH_percentage_unmerged_reads,2),
         `Min Length of Unmerged Reads` = FLASH_unmerged_min_len,
         `Max Length of Unmerged Reads` = FLASH_unmerged_max_len,
         `Average Length of Unmerged Reads` = round(FLASH_unmerged_avg_len,2),
         `Host Removal Number of Input Merged + Unmerged Reads (Total)` = Hostrem_input_total_num_seqs,
         `Host Removal Percentage of Merged + Unmerged Reads Removed (%)` = round(Hostrem_percentage_total_seqs_removed,2),
         `Host Removal Number of Output Merged + Unmerged Reads (Total)` = Hostrem_output_total_num_seqs,
         `Host Removal Number of Input Merged Reads (Total)` = Hostrem_input_merged_num_seqs,
         `Host Removal Percentage of Merged Reads Removed (%)` = round(Hostrem_percentage_merged_seqs_removed,2),
         `Host Removal Number of Output Merged Reads (Total)` = Hostrem_output_merged_num_seqs,
         `Host Removal Number of Input Unmerged Reads (Total)` = Hostrem_input_unmerged_num_seqs,
         `Host Removal Percentage of Unmerged Reads Removed (%)` = round(Hostrem_percentage_unmerged_seqs_removed,2),
         `Host Removal Number of Output Unmerged Reads (Total)` = Hostrem_output_unmerged_num_seqs,
         )%>%
  select(
    SampleID, 
    `Sample Type`, 
    `Library Type`,
    `Trimmmomatic Number of Input Reads (Paired)`, 
    `Trimmmomatic Percentage of Reads Dropped`,
    `Trimmmomatic Number of Output Reads (Paired)`,
    `SortMeRNA Number of Input Reads (Total)`,
    `SortMeRNA Number of Output Reads (Total)`,            
    `SortMeRNA Percentage of Reads Removed`,   
    `FLASH Number of Input Reads (Total)`,
    `FLASH Number of Reads Merged (Total)`,
    `FLASH Percentage of Reads Merged (%)`,
    `Min Length of Merged Reads`,
    `Max Length of Merged Reads`,
    `Average Length of Merged Reads`,
    `FLASH Number of Reads Unmerged (Total)`,
    `FLASH Percentage of Reads Unmerged (%)`,
    `Min Length of Unmerged Reads`,
    `Max Length of Unmerged Reads`,
    `Average Length of Unmerged Reads`,
    `Host Removal Number of Input Merged + Unmerged Reads (Total)`,
    `Host Removal Percentage of Merged + Unmerged Reads Removed (%)`,
    `Host Removal Number of Output Merged + Unmerged Reads (Total)`,
    `Host Removal Number of Input Merged Reads (Total)`,
    `Host Removal Percentage of Merged Reads Removed (%)`,
    `Host Removal Number of Output Merged Reads (Total)`,
    `Host Removal Number of Input Unmerged Reads (Total)`,
    `Host Removal Percentage of Unmerged Reads Removed (%)`,
    `Host Removal Number of Output Unmerged Reads (Total)` 
  )
stable2
  
write_xlsx(stable2,
          "/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Writing/Paper_figures/SupplementaryTable2.xlsx")


