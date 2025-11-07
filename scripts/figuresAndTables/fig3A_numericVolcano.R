# =============================================================================
# VOLCANO PLOT FOR SIGNIFICANT PROTEINS DERIVED FROM LINEAR REGRESSION
# Author: Alice Piller
# Date: 14 May 2024
# =============================================================================

# =============================================================================
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates a volcano plot showing protein associations with a
# numeric condition (e.g., neutralization titers) based on linear regression
# coefficients and FDR-adjusted p-values.
#
# Specifically, it:
#  - Loads a full set of differential results and a filtered significant subset
#  - Classifies proteins as "Up", "Down", or "No change" based on coefficient
#    direction and significance
#  - Filters extreme effect size outliers
#  - Adds protein labels for significant features
#  - Creates and saves a volcano plot as a high-resolution PNG
#
# Output: volcano plot saved to ../../figuresAndTables/fig3/
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

## Load libraries
pacman::p_load(dplyr, ggplot2, ggrepel, ggthemes)

## Set working directory
folder <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(folder, "/sourceDataAndFigures_Jun25/data/proteins"))


# =============================================================================
# PARAMETERS
# =============================================================================

## 1) Condition used for modeling
condition <- "D614GNeutralization"

## 2) FDR threshold for significance
fdrThreshold <- 0.05

## 3) Absolute effect size thresholds for removing outliers
effectThresholdUpper <- 4400
effectThresholdLower <- -5000

## 4) Date tag for output
dateTag <- "170625"


# =============================================================================
# FUNCTIONS
# =============================================================================

## Custom ggplot2 theme for publication-quality plots
themePublication <- function(baseSize = 17, baseFamily = "helvetica") {
  theme_foundation(base_size = baseSize, base_family = baseFamily) +
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

## Load full and filtered protein regression results
proteinsAll <- read.csv("allProteins/D614GNeutralization/fullResultsD614GNeutralization_170625.csv")
proteinsFiltered <- read.csv("threshProteins/D614GNeutralization/threshResultsD614GNeutralizationFDR_0.05_170625.csv")


# =============================================================================
# LABEL DIFFERENTIALLY EXPRESSED PROTEINS
# =============================================================================

## Identify significant proteins
sigProteins <- proteinsFiltered$proteinId

## Label and filter data
proteinsAll <- proteinsAll %>%
  filter(coefficient < effectThresholdUpper & coefficient > effectThresholdLower) %>%
  mutate(diffExpr = if_else(proteinId %in% sigProteins,
                            true = if_else(coefficient > 0, "Up", "Down"),
                            false = "No change"),
         .after = fdr)

## Count by label
classCount <- proteinsAll %>% count(diffExpr)

## Label proteins for display
proteinsAll$deLabel <- NA
highlight <- which(proteinsAll$diffExpr %in% c("Up", "Down"))
proteinsAll$deLabel[highlight] <- proteinsAll$proteinId[highlight]

## Reorder factor for plotting
proteinsAll$diffExpr <- factor(proteinsAll$diffExpr, levels = c("Down", "Up", "No change"))

## Add -log10(fdr)
proteinsAll$negLog10Fdr <- -log10(proteinsAll$fdr)


# =============================================================================
# PLOT
# =============================================================================

## Define color scheme
plotColors <- c("#006EAE", "#C5373D", "#8E99AB")
names(plotColors) <- c("Down", "Up", "No change")

## Volcano plot
volcanoPlot <- ggplot(data = proteinsAll, 
                      aes(x = coefficient, y = negLog10Fdr, color = diffExpr, label = deLabel)) +
  geom_point(size = 2, alpha = 0.65) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.73, alpha = 0.7) +
  geom_hline(yintercept = -log10(fdrThreshold), color = "black", linewidth = 0.73, alpha = 0.7) +
  geom_text_repel(size = 4, seed = 1, fontface = "bold") +
  labs(x = "Effect size", y = "-log10(FDR)") +
  scale_color_manual(values = plotColors) +
  themePublication()

volcanoPlot


# =============================================================================
# SAVE PLOT
# =============================================================================

## Define file name and save
fileName <- paste0("../../figuresAndTables/fig3/", condition, "Volcano_out", 
                   effectThresholdUpper, "_", effectThresholdLower, "_", dateTag, ".png")

ggsave(fileName, volcanoPlot, width = 7, height = 7, dpi = 600)
