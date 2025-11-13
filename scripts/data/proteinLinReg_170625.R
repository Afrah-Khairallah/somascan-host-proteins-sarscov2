# =============================================================================
# LINEAR REGRESSION - PROTEIN ASSOCIATIONS WITH NUMERIC FRNT50
# Author: Alice Piller
# Date: May 2025
# =============================================================================

# =============================================================================
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script identifies proteins associated with continuous measures 
# (e.g., D614GNeutralization or Spike/binding titers) using univariate linear 
# models.
#
# Specifically, it:
#  - Reads cleaned protein expression data and log-transforms protein values
#  - Fits individual linear models for each protein against the chosen response
#  - Extracts coefficients and p-values
#  - Applies Storey FDR correction to p-values
#  - Outputs full and filtered result sets as CSV files
#
# Output: CSV files saved to ../genes/allGenes/ and ../genes/threshGenes/
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

## Set working directory to cleaned data folder
folder <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(folder, "/sourceDataAndFigures_Jun25/data/cleaned"))


# =============================================================================
# PARAMETERS
# =============================================================================

## 1) Response variable for modeling (e.g., "D614GNeutralization")
responseVar <- "D614GNeutralization"

## 2) FDR threshold
fdrThreshold <- 0.05

## 3) Protein column indices
proteinStartCol <- 47
proteinEndCol <- 4779

## 4) Date tag for output files
dateTag <- "170625"


# =============================================================================
# READ DATA
# =============================================================================

## Load data
data <- read.csv("cleanedData_150625.csv")

## Log-transform protein expression values
proteinsLog <- data %>%
  select(c(proteinStartCol:proteinEndCol)) %>%
  log()

## Combine response variable with protein data
modelData <- data %>%
  select({{ responseVar }}) %>%
  cbind(proteinsLog)


# =============================================================================
# UNIVARIATE LINEAR MODELS
# =============================================================================

## Extract protein names
proteins <- colnames(proteinsLog)

## Initialize containers for model outputs
univarModels <- list()
coefficients <- c()
pValues <- c()

## Loop through each protein and fit linear model
for (protein in proteins) {
  formula <- as.formula(paste0(responseVar, " ~ ", protein))
  univarModels[[protein]] <- lm(formula, data = modelData) %>% summary()
  coefficients <- c(coefficients, univarModels[[protein]]$coefficients[2, "Estimate"])
  pValues <- c(pValues, univarModels[[protein]]$coefficients[2, "Pr(>|t|)"])
}

## Apply FDR correction using Storey's method
pAdjusted <- qvalue(pValues)$qvalues

## Combine results into data frame
results <- data.frame(
  proteinId = proteins,
  coefficient = coefficients,
  pValue = pValues,
  fdr = pAdjusted
)

## Filter significant results
resultsFiltered <- results %>%
  filter(fdr < fdrThreshold)


# =============================================================================
# WRITE RESULTS TO CSV
# =============================================================================

## Define output paths
fileAll <- paste0("../proteins/allProteins/", responseVar, "/fullResults", responseVar, "_", dateTag, ".csv")
fileFiltered <- paste0("../proteins/threshProteins/", responseVar, "/threshResults", responseVar, "FDR_", fdrThreshold, "_", dateTag, ".csv")

## Save full results
write.csv(results, fileAll, row.names = FALSE)

## Save filtered results
write.csv(resultsFiltered, fileFiltered, row.names = FALSE)
