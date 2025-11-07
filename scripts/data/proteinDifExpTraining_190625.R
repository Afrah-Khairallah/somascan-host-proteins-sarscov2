# =============================================================================
# DIFFERENTIAL PROTEIN EXPRESSION
# Alice Piller
# May 2024
# =============================================================================


# =============================================================================
# DESCRIPTION
# =============================================================================

# This script performs differential protein expression analysis on protein expression
# data using a t-test. It compares two groups (e.g., neutralization, or and 
# returns both raw and FDR-adjusted significance values, as well as log fold 
# changes. Results are exported as .csv files for downstream use.


# =============================================================================
# SETUP
# =============================================================================

## Clear environment
rm(list = ls())

## Load packages
library(dplyr)
library(qvalue)


# =============================================================================
# PARAMETERS
# =============================================================================

## Working directory
setwd("C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/sourceDataAndFigures_Jun25/data/predictiveModel")

## Significance threshold
threshFDR <- 0.05
threshFC <- 1.4

## Column range containing protein expression values
proteinStartCol <- 47
proteinEndCol <- 4779

## Set date
date <- "190625"


# =============================================================================
# READ DATA
# =============================================================================

data <- read.csv("training_0.6_190625.csv")


# =============================================================================
# SPLIT INTO COMPARISON GROUPS
# =============================================================================

## Assign group1 and group2 based on condition
group1 <- filter(data, Neuts_Binom == "1")
group2 <- filter(data, Neuts_Binom == "0")

## Select only protein columns
group1 <- group1 %>% select(proteinStartCol:proteinEndCol)
group2 <- group2 %>% select(proteinStartCol:proteinEndCol)


# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

## Define function to conduct t-tests on log-transformed protein values
conductTTest <- function(protein) {
  a <- log(group1[[protein]])
  b <- log(group2[[protein]])
  t.test(a, b, var.equal = TRUE)
}

## Apply t-test across all proteins
proteins <- colnames(group1)
tTestResults <- sapply(proteins, conductTTest)

## Extract p-values and adjust for multiple testing
pValues <- unlist(tTestResults["p.value", ])
pAdjusted <- qvalue(pValues)$qvalues
log10FDR <- -log10(pAdjusted)

## Compute log fold change between group means
meanGroup1 <- apply(group1, 2, mean)
meanGroup2 <- apply(group2, 2, mean)
logFoldChange <- log(meanGroup1 / meanGroup2)

## Create results table
resultsMatrix <- data.frame(
  ProteinID = proteins,
  p.value = pValues,
  FDR = pAdjusted,
  log10FDR = log10FDR,
  meanGroup1 = meanGroup1,
  meanGroup2 = meanGroup2,
  logFoldChange = logFoldChange
)

## Filter by FDR
resultsMatrixThresh <- resultsMatrix %>%
  filter(FDR <= threshFDR,
         abs(logFoldChange) > log(1.4))


# =============================================================================
# WRITE RESULTS TO CSV
# =============================================================================

## Export full results and filtered results
file1 <- paste0("proteins/fullResultsTraining_", date, ".csv")
write.csv(resultsMatrix, file1, row.names = FALSE)

file2 <- paste0("proteins/threshResultsTraining_", date, ".csv")
write.csv(resultsMatrixThresh, file2, row.names = FALSE)
