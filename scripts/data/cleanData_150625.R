# =============================================================================
# CLEAN DATA TO PRODUCE DATA SET WITH PARTICIPANT METADATA AND PROTEOMICS
# Author: Alice Piller
# Date: May 2025
# =============================================================================

# =============================================================================
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script prepares a cleaned dataset integrating proteomics (SomaLogic),
# antibody levels, and participant metadata for a COVID-19 neutralization study.
# Key steps include:
#  - Loading and pre-processing raw data sources (SomaLogic proteomics, antibody 
#    binding summary data at two time points, and participant metadata)
#  - Excluding samples based on pre-defined quality and relevance criteria
#  - Correcting SuppO2 classifications for specific participants
#  - Merging with metadata and calculating derived variables (e.g., age binarisation)
#  - Classifying neutralization levels based on median FRNT50
#  - Averaging gene replicate measurements across probes
#  - Creating binary comorbidity indicators
#  - Renaming specific ambiguous protein names
#  - Merging and binarizing antibody values, and computing neutralization potency
#  - Exporting the final cleaned dataset to CSV
# 
# Output: cleanedData_ab_<date>.csv in the `cleaned` folder
# =============================================================================


# =============================================================================
# SETUP
# =============================================================================

## Clear environment
rm(list = ls())

## Install packages (commented out for journal reproducibility)
#install.packages("dplyr")
#install.packages("tidyr")

## Load necessary libraries
library(dplyr)
library(tidyr)

## Set current date tag for output file naming
date <- "150625"

# =============================================================================
# DIRECTORIES
# =============================================================================

## Set working directory to raw data folder inside the project
folder <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(folder, "/sourceDataAndFigures_Jun25/data/raw"))


# =============================================================================
# READ DATA
# =============================================================================

## Load SomaLogic proteomics data
somalogicMaster <- read.csv("somalogic_master_18_05_2023.csv")
matrixCleaned <- somalogicMaster

## Load antibody data from Penny (timepoint 1 and 2)
abTp1 <- read.csv("20250211_OMICS_Sigal_samples_binding_summary_ec50_tp1.csv")
abTp2 <- read.csv("20250211_OMICS_Sigal_samples_binding_summary_ec50_later_tp.csv")

## Load metadata
metadata <- read.csv("subsetfullMetadataWithNeut.csv")


# =============================================================================
# CLEAN DATA
# =============================================================================

# -----------------------------------------------------------------------------
# 1) Exclude unwanted samples based on various quality/eligibility filters
# -----------------------------------------------------------------------------

## Remove PID 039-02-0027 due to immunosuppression
sidforPid27 <- "EDT3-03910171"

## Exclude samples sampled too late - more than 10 days
timepointSomaSid <- c("EDT3-03910012", "EDT3-03910051", "EDT3-03910063", "EDT3-03910066", "EDT3-03910070", "EDT3-03910097", "EDT3-03910135", "EDT3-03910032", "EDT3-03910039", "EDT3-03910038", "EDT3-03910057", "EDT3-03910053", "EDT3-03910055", "EDT3-03910062", "EDT3-03910065", "EDT3-03910064", "EDT3-03910079", "EDT3-03910080", "EDT3-03910084", "EDT3-03910086", "EDT3-03910090", "EDT3-03910103", "EDT3-03910105", "EDT3-03910109", "EDT3-03910112", "EDT3-03910118", "EDT3-03910117", "EDT3-03910119", "EDT3-03910122", "EDT3-03910124", "EDT3-03910115", "EDT3-03910126", "EDT3-03910127", "EDT3-03910128", "EDT3-03910130", "EDT3-03910133", "EDT3-03910134", "EDT3-03910136", "EDT3-03910137", "EDT3-03910141", "EDT3-03910150", "EDT3-03910157", "EDT3-03910159", "EDT3-03910142", "EDT3-03910151", "EDT3-03910155", "EDT3-03910160", "EDT3-03910167", "EDT3-03910144")

## Exclude neuts done at time point 1 (too early) or time point 2 but less than 10 days from confirmed diagnostic swab
timepoint1NeutSid <- c("EDT3-03910015", "EDT3-03910012", "EDT3-03910020", "EDT3-03910034", "EDT3-03910039", "EDT3-03910057", "EDT3-03910053", "EDT3-03910055", "EDT3-03910065", "EDT3-03910064", "EDT3-03910088", "EDT3-03910090", "EDT3-03910091", "EDT3-03910092", "EDT3-03910116", "EDT3-03910114", "EDT3-03910123", "EDT3-03910127", "EDT3-03910138", "EDT3-03910143", "EDT3-03910150", "EDT3-03910158", "EDT3-03910166", "EDT3-03910170")

## Exclude neuts done more than 10 days post-swab from time point 2
timepoint2NeutSid <- c("EDT3-03920026", "EDT3-03920104")

## Combine all sample IDs to exclude
sidsToExclude <- c(sidforPid27, timepointSomaSid, timepoint1NeutSid)

# -----------------------------------------------------------------------------
# 2) Reclassify false SuppO2 cases based on delayed deaths
# -----------------------------------------------------------------------------

falseSuppO2sid <- c("EDT3-03910043", "EDT3-03910161")

# -----------------------------------------------------------------------------
# 3) Apply exclusions and SuppO2 corrections
# -----------------------------------------------------------------------------

matrixCleaned <- matrixCleaned %>%
  filter(!SID %in% sidsToExclude) %>%
  mutate(SuppO2_Binom = ifelse(SID %in% falseSuppO2sid, "0", SuppO2_Binom),
         SuppO2 = ifelse(SID %in% falseSuppO2sid, "N/A", SuppO2))


# =============================================================================
# MERGE WITH METADATA
# =============================================================================

## Merge to add Age at Enrollment and derive binary age group
metadata <- metadata %>%
  select(PID, Age_At_Enrolment.x)

matrixCleaned <- matrixCleaned %>%
  select(-Age_Group) %>%
  left_join(metadata, by = "PID") %>%
  relocate(Age_At_Enrolment.x, .before = AgeGroup_Binom) %>%
  mutate(AgeGroup_Binom = if_else(Age_At_Enrolment.x > 49, 1, 0))


# =============================================================================
# CLASSIFY NEUTRALIZATION LEVELS
# =============================================================================

## Use median FRNT50 to classify neutralization response as high/low
threshFrnt50 <- median(matrixCleaned$D614GNeutralization)

matrixCleaned <- matrixCleaned %>%
  mutate(Neuts_Binom = if_else(D614GNeutralization > threshFrnt50, "1", "0"))


# =============================================================================
# COLLAPSE GENE REPLICATES
# =============================================================================

## Convert wide gene columns to long, average replicate probes, convert back
genes <- matrixCleaned %>%
  select(PID, CRYBB2:last_col()) # Make sure CRYBB2 is first protein column

longData <- genes %>%
  pivot_longer(cols = -PID, names_to = "gene", values_to = "expression") %>%
  mutate(gene_base = gsub("_\\d+$", "", gene))

averagedData <- longData %>%
  group_by(PID, gene_base) %>%
  summarise(mean_expression = mean(expression), .groups = 'drop')

wideData <- averagedData %>%
  pivot_wider(names_from = gene_base, values_from = mean_expression)

matrixCleaned_no_genes <- matrixCleaned %>% select(D614GNeutralization:SubjectID)
matrixCleaned <- merge(matrixCleaned_no_genes, wideData, by = "PID")


# =============================================================================
# ADD COMORBIDITY BINARIZATION
# =============================================================================

## Create Comorbid_Binom where 1 = Hypertension or Diabetes
matrixCleaned <- matrixCleaned %>%
  mutate(Comorbid_Binom = if_else(Hyp_Binom == 1 | Diabetes_Binom == 1, 1, 0), .before = Diabetes_)


# =============================================================================
# CORRECT COMPOUND PROTEIN NAME
# =============================================================================

matrixCleaned <- matrixCleaned %>%
  rename(FTH1 = FTH1FTL)


# =============================================================================
# ADD ANTIBODY DATA
# =============================================================================

## Reformat antibody data and add binary classification
abTp1 <- abTp1 %>%
  mutate(SID = paste0("EDT3-0", SID)) %>%
  rename_with(~ paste0(., "_tp1"), c("Spike.D614G", "RBD.D614G", "Nucleocapsid"))
abTp2 <- abTp2 %>%
  mutate(SID = paste0("EDT3-0", SID)) %>%
  rename_with(~ paste0(., "_tp2"), c("Spike.D614G", "RBD.D614G", "Nucleocapsid"))

## Merge and compute Potency = FRNT50 / Spike binding
matrixCleaned <- matrixCleaned %>%
  left_join(abTp1, by = c("PID", "SID")) %>%
  left_join(abTp2, by = c("PID")) %>%
  mutate(Spike.D614G_tp1_Binom = if_else(Spike.D614G_tp1 > median(Spike.D614G_tp1, na.rm = T), 1, 0),
         RBD.D614G_tp1_Binom = if_else(RBD.D614G_tp1 > median(RBD.D614G_tp1, na.rm = T), 1, 0),
         Nucleocapsid_tp1_Binom = if_else(Nucleocapsid_tp1 > median(Nucleocapsid_tp1, na.rm = T), 1, 0),
         Spike.D614G_tp2_Binom = if_else(Spike.D614G_tp2 > quantile(Spike.D614G_tp2, 0.4, na.rm = T), 1, 0), # 0.4 quantile gave most significantly differential proteins
         RBD.D614G_tp2_Binom = if_else(RBD.D614G_tp2 > quantile(RBD.D614G_tp2, 0.4), 1, 0), # 0.4 quantile gave most significantly differential proteins
         Nucleocapsid_tp2_Binom = if_else(Nucleocapsid_tp2 > median(Nucleocapsid_tp2, na.rm = T), 1, 0)) %>%
  relocate(c("Spike.D614G_tp1", "Spike.D614G_tp1_Binom", 
             "RBD.D614G_tp1", "RBD.D614G_tp1_Binom", 
             "Nucleocapsid_tp1", "Nucleocapsid_tp1_Binom",
             "Spike.D614G_tp2", "Spike.D614G_tp2_Binom",
             "RBD.D614G_tp2", "RBD.D614G_tp2_Binom",
             "Nucleocapsid_tp2", "Nucleocapsid_tp2_Binom"), 
           .after = 31) %>%
  relocate(SID.y, .after = SID.x) %>%
  mutate(Potency = D614GNeutralization/Spike.D614G_tp2,
         Potency_Binom = if_else(Potency > quantile(Potency, 0.6), 1, 0), # 0.6 quantile gave most significantly differential proteins
         .after = Nucleocapsid_tp2)


# =============================================================================
# WRITE CSV
# =============================================================================

filename <- paste0("../cleaned/cleanedData_", date, ".csv")
write.csv(matrixCleaned, filename, row.names = FALSE)
