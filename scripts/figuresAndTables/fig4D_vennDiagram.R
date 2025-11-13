# =============================================================================
# VENN DIAGRAM
# Author: Alice Piller
# Date: June 2024
# =============================================================================

# =============================================================================
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates a Venn diagram of the top differentially expressed
# proteins between two user-specified conditions (e.g., neutralization and severity).
# It uses color-coded overlapping sets and customized label positions.
#
# Output: A high-resolution PNG image of the Venn diagram.
# =============================================================================


# =============================================================================
# SETUP
# =============================================================================

## Clear environment
rm(list = ls())

## Load libraries
pacman::p_load(VennDiagram, ggplot2, dplyr, colorspace)

## Set working directory
folder <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(folder, "/sourceDataAndFigures_Jun25/data/proteins/threshProteins"))


# =============================================================================
# PARAMETERS
# =============================================================================

## Conditions to compare
conds <- c("Spike.D614G_tp2_Binom", "Neuts_Binom", "D614GNeutralization")

## File paths for each condition
filePaths <- list(
  Neuts_Binom = "Neuts_Binom/threshResultsNeuts_Binom_FDR0.05_170625.csv",
  D614GNeutralization = "D614GNeutralization/threshResultsD614GNeutralizationFDR_0.05_170625.csv",
  SuppO2_Binom = "SuppO2_Binom/threshResultsSuppO2_Binom_FDR0.05_170625.csv",
  Spike.D614G_tp2_Binom = "Spike.D614G_tp2_Binom/threshResultsSpike.D614G_tp2_Binom_FDR0.05_170625.csv"
)

## Plot labels for each condition
labelLookup <- list(
  SuppO2_Binom = expression(bold("Severity")),
  Neuts_Binom = expression(bold("High/Low FRNT"["50"])),
  D614GNeutralization = expression(bold("Numeric FRNT"["50"])),
  Spike.D614G_tp2_Binom = expression(bold("Anti-Spike Ab"))
)

## Fill colors for each condition
colorLookup <- list(
  SuppO2_Binom = "#E96900",
  Neuts_Binom = "#7e45a4",
  D614GNeutralization = "#be456e",
  Spike.D614G_tp2_Binom = "#3D892E"
)

## Date tag for file naming
dateTag <- "190625"


# =============================================================================
# READ DATA
# =============================================================================

## Load significant proteins (FDR < 0.05) for both conditions
data1 <- read.csv(filePaths[[conds[1]]])
data2 <- read.csv(filePaths[[conds[2]]])
data3 <- read.csv(filePaths[[conds[3]]])

## Extract protein IDs into named list
proteinLists <- list(
  a = data1$proteinId,
  b = data2$proteinId,
  c = data3$proteinId
)


# =============================================================================
# VENN DIAGRAM
# =============================================================================

## Color and label setup
fillColors <- unlist(colorLookup[conds])
categoryLabels <- c(labelLookup[[conds[1]]], labelLookup[[conds[2]]], labelLookup[[conds[3]]])

## Customize label position and distance between conditions
labelDist <- list(Neuts_Binom = 0.07,
                  D614GNeutralization = 0.1,
                  SuppO2_Binom = 0.06,
                  Spike.D614G_tp2_Binom = 0.065)

labelPos <- list(Neuts_Binom = 150,
                 D614GNeutralization = 0,
                 SuppO2_Binom = 45,
                 Spike.D614G_tp2_Binom = 200)

## Extract distances and positions based on conditions
dist <- c(labelDist[[conds[1]]], labelDist[[conds[2]]], labelDist[[conds[3]]])
position <- c(labelPos[[conds[1]]], labelPos[[conds[2]]], labelPos[[conds[3]]])

## Generate Venn diagram
vennPlot <- venn.diagram(
  x = proteinLists,
  category.names = categoryLabels,
  filename = NULL,
  fill = fillColors,
  cat.col = darken(fillColors),
  alpha = 0.5,
  lwd = 0,
  cat.fontfamily = "arial",
  cat.fontface = "bold",
  fontface = "italic",
  cat.cex = 2.8,
  cex = 2.8,
  cat.pos = position,
  cat.dist = dist,
  margin = 0.19,
  disable.logging = TRUE
)

## Clear page and draw
grid.newpage()
grid.draw(vennPlot)


# =============================================================================
# SAVE PLOT
# =============================================================================

## Define filename and save
fileName <- paste0("../../../figuresAndTables/fig4/fig4D_venn_", conds[1], "_", conds[2], "_", dateTag, ".png")
ggsave(fileName, vennPlot, width = 10, height = 10, dpi = 600)
