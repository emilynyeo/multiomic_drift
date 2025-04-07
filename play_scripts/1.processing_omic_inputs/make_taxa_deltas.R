# Make taxa Deltas

# Process Picrust CHANGE data 

library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
#install.packages("phangorn")
#install.packages("caper",repos="https://cloud.r-project.org",quiet=TRUE)
#install.packages("diversitree",repos="https://cloud.r-project.org",quiet=TRUE)
#install.packages("geiger",repos="https://cloud.r-project.org",quiet=TRUE)
#install.packages("nlme",repos="https://cloud.r-project.org",quiet=TRUE)
#install.packages("OUwie",repos="https://cloud.r-project.org",quiet=TRUE)
#install.packages("phangorn",repos="https://cloud.r-project.org",quiet=TRUE)
#install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
#devtools::install_github(repo = "aplantin/pldist")
library(pldist)

make_unique_names <- function(names) {
  name_counts <- table(names)  # Create a table of name counts
  occurrence <- integer(length(names)) # Vector to track occurrences
  
  # Loop through each name and assign a unique suffix if necessary
  for (i in seq_along(names)) {
    if (name_counts[names[i]] > 1) {
      occurrence[i] <- sum(names[1:i] == names[i])  # Count occurrences up to current index
      names[i] <- paste0(names[i], "_", occurrence[i])  # Append the occurrence number
    }
  }
  
  # Ensure unique names
  for (i in seq_along(names)) {
    while (sum(names[i] == names) > 1) {
      occurrence[i] <- occurrence[i] + 1
      names[i] <- sub("_\\d+$", "", names[i])  # Remove existing suffix
      names[i] <- paste0(names[i], "_", occurrence[i])  # Append new occurrence number
    }
  }
  return(names)
}
# Function to process the names ALL
process_names_all <- function(names) {
  cleaned_numbers_all <- character(length(names)) #vector to store cleaned numbers
  all_names <- grep("\\..*$", names, value = TRUE) #names matching the pattern
  for (i in seq_along(names)) { # Extract numbers and remove leading zeros
    if (names[i] %in% all_names) {
      extracted_number <- sub(".*-(\\d+)\\..*$", "\\1", names[i])
      cleaned_number <- sub("^0+", "", extracted_number)
      cleaned_numbers_all[i] <- cleaned_number
    } else {
      cleaned_numbers_all[i] <- names[i] # If name don't match, keep the original
    }
  }
  return(cleaned_numbers_all)  # Return the cleaned numbers
}
print_taxa_summary <- function(before, after, step) {
  removed <- before - after
  print(paste0("# Taxa removed in step ", step, ": ", removed))
  print(paste0("# Taxa remaining after step ", step, ": ", after))
}
## TAXA ####################################################################################################
load("~/projects/research/Stanislawski/BMI_risk_scores/microbiome_rs/data/PhyloseqObj.RData")
load("/Users/emily/projects/research/Stanislawski/BMI_risk_scores/microbiome_rs/data/Genus_Sp_tables.RData")
print_ps(drift.phy.count)

# Step 1: Remove taxa not seen more than 3 times in at least 10% of the samples
initial_taxa_count <- ntaxa(drift.phy.count)
GP_count <- filter_taxa(drift.phy.count, function(x) sum(x > 3) > (0.1*length(x)), TRUE)
print_taxa_summary(initial_taxa_count, ntaxa(GP_count), "1")

# Step 2: Apply coefficient of variation cutoff
initial_taxa_count_gp <- ntaxa(GP_count)
gpsf_count <- filter_taxa(GP_count, function(x) sd(x)/mean(x) > 1.0, TRUE)
print_taxa_summary(initial_taxa_count_gp, ntaxa(gpsf_count), "2")

# Get genus --- double check MATTS email 
tax_genus <- tax_glom(gpsf_count, "Genus")
genus_filtered_count <- otu_table(tax_genus) 
genus_filtered_count <- t(genus_filtered_count) %>% as.data.frame()
colnames(genus_filtered_count) <- as.data.frame(tax_table(tax_genus))$Genus
rm(drift.phy.clr,drift.phy.count,GP_count,gpsf_count,sp.clr, tax_genus,sp.count,
   drift.phy.count.r21116, drift.phy.ra, genus.clr, genus.count, genus.ra, sp.ra)

# Function to make column names unique
colnames(genus_filtered_count) <- make_unique_names(colnames(genus_filtered_count)) # rename the columns
genus_count <- genus_filtered_count %>% rownames_to_column(var = "subject_id") 

# Relative abundance conversion
#Genus_relative_abundance <- genus_filtered_count %>%
#  rownames_to_column(var = "subject_id") %>%  dplyr::rowwise() %>%
#  mutate(total = sum(c_across(g__Parabacteroides_B_862066:`g__Massilistercora`), na.rm = TRUE)) %>%
#  mutate(across(g__Parabacteroides_B_862066:`g__Massilistercora`, ~ .x / total)) %>%
#  dplyr::select(-total) %>% column_to_rownames(var = "subject_id")

##  removes any column where the proportion of zeros is greater than or equal to 80%.
threshold <- 0.80
zero_percentage <- colSums(genus_count == 0, na.rm = TRUE) / nrow(genus_count) # % zeros / column
genus_count_cleaned <- genus_count[, zero_percentage < threshold] # only colz < 20% zeros
dim(genus_count) - dim(genus_count_cleaned) #

heatmap(cor(genus_count_cleaned[, 2:ncol(genus_count_cleaned)]))
genus_count_cleaned$all_samples <- process_names_all(genus_count_cleaned$subject_id)

# make time column 
genus_count_cleaned <- genus_count_cleaned %>%
  mutate(time = case_when(grepl("BL", subject_id) ~ 0, grepl("6m", subject_id) ~ 6, 
                          grepl("12m", subject_id) ~ 12, TRUE ~ NA_real_)) %>% 
                          filter(!grepl("\\.3m$", subject_id)) %>% 
                          filter(!grepl("\\.18m$", subject_id)) 

### First make 0-6 m change
tax_0_6 <- genus_count_cleaned %>% dplyr::filter(time !=12)
tax_0_6$subjID <- gsub("\\.(BL|6m|12m)$", "", tax_0_6$subject_id)
tax_0_6 <- tax_0_6 %>% filter(subjID %in% names(table(subjID)[table(subjID) == 2]))

meta_tax <- tax_0_6 %>% dplyr::select(c(subjID, subject_id, time))
meta_tax <- meta_tax %>% filter(subjID %in% names(table(subjID)[table(subjID) == 2]))

# Filter meta data and arrange by subject_id
meta.pldist <- meta_tax %>%
  dplyr::rename(sampID = subject_id) %>%
  dplyr::filter(sampID %in% tax_0_6$subject_id, time %in% c(0, 6)) %>%
  #dplyr::filter(subjID %in% n.2) %>% 
  arrange(sampID)

# Prepare OTU data and check matching samples
otu.pldist <- tax_0_6 %>%
  dplyr::filter(subject_id %in% meta.pldist$sampID) %>%
  arrange(subject_id) %>% column_to_rownames("subject_id") %>% 
  dplyr::select(-c("all_samples", "time", "subjID"))

# Prepare OTU data for analysis
otu.data <- data_prep(otu.pldist, meta.pldist, paired = TRUE)
otu.data$otu.props[1:3, 1:3]  # OTU proportions
otu.data$otu.clr[1:3, 1:3]    # CLR-transformed proportions

# Perform transformation
res <- pltransform(otu.data, paired = TRUE, norm = TRUE)

# Assign results
quant.prop <- res$dat.quant.prop
quant.clr <- res$dat.quant.clr

save(quant.prop, file='~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/tax_0_6_pldist_RA.RData')
save(quant.clr , file='~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/tax_0_6_pldist_CLR.RData')

##########################################################

### First make 6-12 m change
tax_6_12 <- genus_count_cleaned %>% dplyr::filter(time !=0)
tax_6_12$subjID <- gsub("\\.(BL|6m|12m)$", "", tax_6_12$subject_id)
tax_6_12 <- tax_6_12 %>% filter(subjID %in% names(table(subjID)[table(subjID) == 2]))

meta_tax <- tax_6_12 %>% dplyr::select(c(subjID,subject_id, time))
meta_tax <- meta_tax %>% filter(subjID %in% names(table(subjID)[table(subjID) == 2]))

# Filter meta data and arrange by subject_id
meta.pldist <- meta_tax %>%
  dplyr::rename(sampID = subject_id) %>%
  dplyr::filter(sampID %in% tax_6_12$subject_id, time %in% c(6, 12)) %>%
  arrange(sampID)

# Prepare OTU data and check matching samples
otu.pldist <- tax_6_12 %>%
  dplyr::filter(subject_id %in% meta.pldist$sampID) %>%
  arrange(subject_id) %>% column_to_rownames("subject_id") %>% 
  dplyr::select(-c("all_samples", "time", "subjID"))

# Prepare OTU data for analysis
otu.data <- data_prep(otu.pldist, meta.pldist, paired = TRUE)
otu.data$otu.props[1:3, 1:3]  # OTU proportions
otu.data$otu.clr[1:3, 1:3]    # CLR-transformed proportions

# Perform transformation
res <- pltransform(otu.data, paired = TRUE, norm = TRUE)

# Assign results
quant.prop <- res$dat.quant.prop
quant.clr <- res$dat.quant.clr

save(quant.prop, file='~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/tax_6_12_pldist_RA.RData')
save(quant.clr , file='~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/tax_6_12_pldist_CLR.RData')
