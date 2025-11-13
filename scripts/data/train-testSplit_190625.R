# =============================================================================
# TRAIN-TEST SPLIT
# Alice Piller
# July 2024
# =============================================================================


# =============================================================================
# DESCRIPTION
# =============================================================================

# This script reads in the cleaned dataset and generates stratified training and
# testing datasets based on the binary neutralization variable (Neuts_Binom).
# The output is saved to CSV for downstream predictive modeling.


# =============================================================================
# SETUP
# =============================================================================

## Clear environment
rm(list = ls())

## Set working directory
setwd("C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/sourceDataAndFigures_Jun25/data")

## Load packages
library(rsample)

## Ensure dplyr functions take precedence
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# =============================================================================
# PARAMETERS
# =============================================================================

## Proportion of data to allocate to training set
proportion <- 0.6

## Output date tag
date <- "190625"


# =============================================================================
# READ DATA
# =============================================================================

## Load cleaned dataset
data <- read.csv("./cleaned/cleanedData_150625.csv")


# =============================================================================
# TRAIN-TEST SPLIT
# =============================================================================

## Set seed for reproducibility
set.seed(123)

## Perform stratified split on Neuts_Binom variable
split <- initial_split(data, prop = proportion, strata = Neuts_Binom)
train <- training(split)
test <- testing(split)


# =============================================================================
# SAVE OUTPUT FILES
# =============================================================================

## Write training and testing datasets to CSV
write.csv(train, paste0("./predictiveModel/training_", proportion, "_", date, ".csv"), row.names = FALSE)
write.csv(test, paste0("./predictiveModel/testing_", proportion, "_", date, ".csv"), row.names = FALSE)
