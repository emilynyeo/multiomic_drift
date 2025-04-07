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

## Functional ####################################################################################################
pathways <- fread("/Users/emily/projects/research/Stanislawski/BMI_risk_scores/picrust2/june7/pathways_out/path_abun_unstrat_descrip.tsv") %>% .[, -1] %>% 
  t() %>% as.data.frame() %>% row_to_names(1) %>% rownames_to_column("SampleID") %>% mutate(across(-1, as.numeric))

# Remove variables with low presence thresholds 
threshold <- 0.80
zero_percentage <- colSums(pathways == 0, na.rm = TRUE) / nrow(pathways) # % zeros / column
path_all_time_cleaned <- pathways[, zero_percentage < threshold] # only colz < 20% zeros
dim(pathways) - dim(path_all_time_cleaned) 
any(path_all_time_cleaned < 0, na.rm = TRUE)

# MAYBE CLR TRANSFORMATION for PATHWAYs
#path_clr_transformed <- cbind(path_all_time_cleaned[, 1, drop = FALSE], 
#                              apply(path_all_time_cleaned[, -1], 2, clr) %>% as.data.frame())

# remove low varianve and high correlation 
preProcValues <- preProcess(path_all_time_cleaned[,2:312], 
                            method = c("nzv", "corr"), thresh = 0.95, fudge = 0.2,
                            na.remove = TRUE, numUnique = 15, verbose = TRUE, 
                            freqCut = 95/5, uniqueCut = 10, cutoff = 0.95)
preProcValues
path_count_cs <- predict(preProcValues, path_all_time_cleaned)
heatmap(cor(path_count_cs[, 2:147]))

# Process samples 
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
path_count_cs$all_samples <- process_names_all(path_count_cs$SampleID)
colnames(path_count_cs) <- make.names(colnames(path_count_cs), unique = TRUE)
# make time column 
path_count_cs <- path_count_cs %>%
  mutate(time = case_when(grepl("BL", SampleID) ~ 0, grepl("6m", SampleID) ~ 6, 
                          grepl("12m", SampleID) ~ 12, TRUE ~ NA_real_)) %>% filter(!grepl("\\.3m$", SampleID)) 

### First make 0-6 m change
path_0_6 <- path_count_cs %>% dplyr::filter(time !=12)
path_0_6$subjID <- gsub("\\.(BL|6m|12m)$", "", path_0_6$SampleID)
path_0_6 <- path_0_6 %>% filter(subjID %in% names(table(subjID)[table(subjID) == 2]))

meta_path <- path_0_6 %>% dplyr::select(c(subjID,SampleID, time))
meta_path <- meta_path %>% filter(subjID %in% names(table(subjID)[table(subjID) == 2]))

# Filter meta data and arrange by SampleID
meta.pldist <- meta_path %>%
  dplyr::rename(sampID = SampleID) %>%
  dplyr::filter(sampID %in% path_0_6$SampleID, time %in% c(0, 6)) %>%
  #dplyr::filter(subjID %in% n.2) %>% 
  arrange(sampID)

# Prepare OTU data and check matching samples
otu.pldist <- path_0_6 %>%
  dplyr::filter(SampleID %in% meta.pldist$sampID) %>%
  arrange(SampleID) %>% column_to_rownames("SampleID") %>% 
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

save(quant.prop, file='~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/pathway_0_6_pldist_RA.RData')
save(quant.clr , file='~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/pathway_0_6_pldist_CLR.RData')

##########################################################

### First make 6-12 m change
path_6_12 <- path_count_cs %>% dplyr::filter(time !=0)
path_6_12$subjID <- gsub("\\.(BL|6m|12m)$", "", path_6_12$SampleID)
path_6_12 <- path_6_12 %>% filter(subjID %in% names(table(subjID)[table(subjID) == 2]))

meta_path <- path_6_12 %>% dplyr::select(c(subjID,SampleID, time))
meta_path <- meta_path %>% filter(subjID %in% names(table(subjID)[table(subjID) == 2]))

# Filter meta data and arrange by SampleID
meta.pldist <- meta_path %>%
  dplyr::rename(sampID = SampleID) %>%
  dplyr::filter(sampID %in% path_6_12$SampleID, time %in% c(6, 12)) %>%
  arrange(sampID)

# Prepare OTU data and check matching samples
otu.pldist <- path_6_12 %>%
  dplyr::filter(SampleID %in% meta.pldist$sampID) %>%
  arrange(SampleID) %>% column_to_rownames("SampleID") %>% 
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

save(quant.prop, file='~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/pathway_6_12_pldist_RA.RData')
save(quant.clr , file='~/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa_change/Change Data/pathway_6_12_pldist_CLR.RData')
