# =============================================================================
# DIFFERENTIAL PROTEIN EXPRESSION
# Author: Alice Piller
# Date: 14 May 2025
# =============================================================================

# =============================================================================
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script performs differential protein expression analysis using SomaLogic
# proteomics data. Specifically, it:
#
#  - Splits participants into two comparison groups based on a binary variable 
#    (e.g., high vs low neutralization capacity (Neuts_Binom))
#  - Conducts two-sample t-tests for each protein after log-transformation
#  - Adjusts p-values for multiple testing using the Storey FDR method
#  - Calculates log-transformed fold changes between groups
#  - Outputs two CSV files:
#      (1) All proteins tested with p-values, FDR, and logFC
#      (2) A filtered set of significantly differentially expressed proteins 
#          (FDR < 0.05)
#
# Output: CSV files saved to ../proteins/allProteins/<condition>/
# =============================================================================


# =============================================================================
# SETUP
# =============================================================================

## Clear environment
rm(list = ls())

## Install packages (commented for reproducibility)
# install.packages("dplyr")
# install.packages("qvalue")

## Load libraries
library(dplyr)
library(qvalue)


# =============================================================================
# DIRECTORIES
# =============================================================================

## Set working directory to folder with cleaned data
folder <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(folder, "/sourceDataAndFigures_Jun25/data/cleaned"))


# =============================================================================
# PARAMETERS
# =============================================================================

## 1) Define condition for group comparison
cond <- "Spike.D614G_tp2_Binom"

## 2) Define FDR threshold
threshFDR <- 0.05

## 3) Define range of protein columns in dataset
protein_start_col <- 47
protein_end_col <- 4779

## 4) Date tag for outputs
date <- "170625"


# =============================================================================
# READ DATA
# =============================================================================

## Load cleaned data
data <- read.csv("cleanedData_150625.csv")


# =============================================================================
# SPLIT DATA INTO COMPARISON GROUPS
# =============================================================================

## 1) Subset rows by binary condition
group1 <- data %>% filter(.data[[cond]] == 1)
group2 <- data %>% filter(.data[[cond]] == 0)

## 2) Select protein expression columns only
group1 <- group1 %>% select(c(protein_start_col:protein_end_col))
group2 <- group2 %>% select(c(protein_start_col:protein_end_col))


# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

## 1) Perform two-sample t-tests on log-transformed expression values
conductTTest <- function(protein) {
  a <- log(group1[[protein]])
  b <- log(group2[[protein]])
  t.test(a, b, var.equal = TRUE)
}

## List of protein column names
proteins <- colnames(group1)

## Run t-tests for each protein and store results
tTestResults <- sapply(proteins, conductTTest)


# =============================================================================
# CALCULATE SIGNIFICANCE AND FOLD CHANGE
# =============================================================================

## Extract p-values
pValues <- unlist(tTestResults["p.value", ])

## Adjust p-values using Storey's FDR method
pAdjusted <- qvalue(pValues)$qvalues

## Convert to -log10 scale
log10FDR <- -log10(pAdjusted)

## Calculate group means for each protein
meanGroup1 <- apply(group1, 2, mean)
meanGroup2 <- apply(group2, 2, mean)

## Compute log fold change
logFoldChange <- log(meanGroup1 / meanGroup2)

## Combine results
resultsMatrix <- data.frame(
  ProteinID = proteins,
  p.value = pValues,
  FDR = pAdjusted,
  log10FDR = log10FDR,
  meanGroup1 = meanGroup1,
  meanGroup2 = meanGroup2,
  logFoldChange = logFoldChange
)

## Filter by FDR threshold
resultsMatrixThresh <- resultsMatrix %>%
  filter(FDR <= threshFDR)


# =============================================================================
# WRITE RESULTS TO CSV
# =============================================================================

## Define output filenames
file1 <- paste0("../proteins/allProteins/", cond, "/fullResults", cond, "_", date, ".csv")
file2 <- paste0("../proteins/threshProteins/", cond, "/threshResults", cond, "_FDR", threshFDR, "_", date, ".csv")

## Export all results
write.csv(resultsMatrix, file1, row.names = FALSE)

## Export filtered results
write.csv(resultsMatrixThresh, file2, row.names = FALSE)
