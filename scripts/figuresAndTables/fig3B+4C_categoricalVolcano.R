# =============================================================================
# VOLCANO PLOT FOR SIGNIFICANT PROTEINS DERIVED FROM LINEAR REGRESSION
# Author: Alice Piller
# Date: June 2024
# =============================================================================

# =============================================================================
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates volcano plots of differentially expressed proteins 
# associated with a binary user-specified condition (e.g., SuppO2_Binom).
#
# It supports highlighting overlap with other protein-level results from up to
# two specified comparison conditions.
#
# Input:
# - Full and filtered differential expression result CSVs
# - Condition of interest and overlap conditions
#
# Output:
# - Publication-ready volcano plot PNG
# =============================================================================


# =============================================================================
# SETUP
# =============================================================================

## Clear environment
rm(list = ls())

## Install packages (commented out for reproducibility)
# install.packages("pacman")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("ggthemes")

## Load packages
pacman::p_load(dplyr, ggplot2, ggrepel, ggthemes)

## Set working directory
folder <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(folder, "/sourceDataAndFigures_Jun25/data/proteins"))


# =============================================================================
# PARAMETERS
# =============================================================================

condition <- "SuppO2_Binom"
overlapConditions <- c("Neuts_Binom", "D614GNeutralization")
fdrThreshold <- 0.05
dateTag <- "170625"


# =============================================================================
# FUNCTIONS
# =============================================================================

themePublication <- function(base_size = 17, base_family = "helvetica") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1), hjust = 0.5),
      text = element_text(face = "bold"),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = 16),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(), 
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_line(colour = "#f0f0f0"),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "lines"),
      legend.margin = margin(0, 0, 0, 0, unit = "cm"),
      legend.title = element_blank(),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
}


# =============================================================================
# READ DATA
# =============================================================================

## Define filepaths for full and significant result sets
filePathsFull <- list(
  SuppO2_Binom = "allProteins/SuppO2_Binom/fullResultsSuppO2_Binom_170625.csv",
  Neuts_Binom = "allProteins/Neuts_Binom/fullResultsNeuts_Binom_170625.csv",
  Spike.D614G_tp2_Binom = "allProteins/Spike.D614G_tp2_Binom/fullResultsSpike.D614G_tp2_Binom_170625.csv"
)

filePathsThresh <- list(
  SuppO2_Binom = "threshProteins/SuppO2_Binom/threshResultsSuppO2_Binom_FDR0.05_170625.csv",
  Neuts_Binom = "threshProteins/Neuts_Binom/threshResultsNeuts_Binom_FDR0.05_170625.csv",
  Spike.D614G_tp2_Binom = "threshProteins/Spike.D614G_tp2_Binom/threshResultsSpike.D614G_tp2_Binom_FDR0.05_170625.csv",
  D614GNeutralization = "threshProteins/D614GNeutralization/threshResultsD614GNeutralizationFDR_0.05_170625.csv"
)

## Load primary condition datasets
dfAll <- read.csv(filePathsFull[[condition]])
dfThresh <- read.csv(filePathsThresh[[condition]])

## Load overlap sets
overlapLists <- lapply(overlapConditions, function(cond) read.csv(filePathsThresh[[cond]])$ProteinID)


# =============================================================================
# LABEL DIFFERENTIALLY EXPRESSED PROTEINS
# =============================================================================

## Identify categories: Up, Down, Shared, or No change
dfAll <- dfAll %>%
  mutate(
    diffExpr = case_when(
      ProteinID %in% dfThresh$ProteinID & logFoldChange > 0 ~ "Up",
      ProteinID %in% dfThresh$ProteinID & logFoldChange < 0 ~ "Down",
      TRUE ~ "No change"
    )
  )

## Overwrite with shared labels if protein is shared
for (i in seq_along(overlapLists)) {
  sharedProteins <- overlapLists[[i]]
  sharedLabel <- "Shared (neut.)"
  dfAll$diffExpr[dfAll$ProteinID %in% sharedProteins & dfAll$ProteinID %in% dfThresh$ProteinID] <- sharedLabel
}

## Label proteins for plotting
labelProteins <- c("Up", "Down", "Shared (neut.)")
dfAll$label <- ifelse(dfAll$diffExpr %in% labelProteins, dfAll$ProteinID, NA)
dfAll$diffExpr <- factor(dfAll$diffExpr, levels = unique(c("Down", "Up", "Shared (neut.)", "No change")))
dfAll$log10Fdr <- -log10(dfAll$FDR)


# =============================================================================
# PLOT
# =============================================================================

plotColors <- c("#006EAE", "#C5373D", "#3D892E", "#8E99AB")
names(plotColors) <- levels(dfAll$diffExpr)

p <- ggplot(dfAll, aes(x = logFoldChange, y = log10Fdr, color = diffExpr, label = label)) +
  geom_point(size = 2, alpha = 0.65) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.73, alpha = 0.7) +
  geom_hline(yintercept = -log10(fdrThreshold), color = "black", linewidth = 0.73, alpha = 0.7) +
  geom_text_repel(size = 4, seed = 1, fontface = "bold") +
  labs(x = "log(fold change)",
       y = "-log10FDR") +
  scale_color_manual(values = plotColors) +
  themePublication()

print(p)


# =============================================================================
# SAVE PLOT
# =============================================================================

outputName <- paste0(condition, "+", paste(overlapConditions, collapse = "+"),
                     "_FDR", fdrThreshold, "_", dateTag, ".png")
outputPath <- paste0("../../figuresAndTables/fig5", outputName)

ggsave(outputPath, p, width = 7, height = 7, dpi = 600)
