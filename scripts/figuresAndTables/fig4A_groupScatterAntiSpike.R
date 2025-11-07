# =============================================================================
# GROUPED SCATTER PLOT - ANTI-SPIKE ANTIBODY
# Author: Alice Piller
# Date: June 2024
# =============================================================================

# =============================================================================
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates a grouped scatter plot of Anti-Spike antibody levels
# stratified by response category ("High" vs "Low") using the 40th percentile
# as the threshold.
#
# The plot includes:
#  - Shaded areas showing the binarization threshold
#  - Scatter points for each participant
#  - Geometric mean titers (GMT) and standard deviation bars
#  - A side density plot
#
# Output: A high-resolution PNG saved to /figuresAndTables/figure4_antiSpikeAntibody/A/
# =============================================================================


# =============================================================================
# SETUP
# =============================================================================

## Clear environment
rm(list = ls())

## Load libraries
pacman::p_load(dplyr, ggplot2, ggsignif, cowplot, stringr, ggthemes, colorspace, ggbeeswarm)

## Set working directory
folder <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(folder, "/sourceDataAndFigures_Jun25/data/cleaned"))


# =============================================================================
# PARAMETERS
# =============================================================================

## Define color palette
colorHigh <- "#e96900"
colorLow <- "#a54891"


# =============================================================================
# FUNCTIONS
# =============================================================================

## Custom publication-quality theme
themePublication <- function(baseSize = 17, baseFamily = "helvetica") {
  theme_foundation(base_size = baseSize, base_family = baseFamily) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
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
      legend.position = "none",
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "lines"),
      legend.margin = margin(0, 0, 0, 0, unit = "cm"),
      legend.title = element_blank(),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
}

## Format scientific notation for axis labels
formatScientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e\\+", "e", l)
  l <- gsub("e", "10^", l)
  parse(text = l)
}


# =============================================================================
# READ DATA
# =============================================================================

data <- read.csv("cleanedData_150625.csv")


# =============================================================================
# DATA PREPARATION
# =============================================================================

## Convert binary response to label
responseVar <- "Spike.D614G_tp2_Binom"
data[[responseVar]] <- ifelse(data[[responseVar]] == 1, "High", "Low")


# =============================================================================
# PLOT SETUP
# =============================================================================

## Compute geometric means and bounds
geomMeans <- data %>%
  group_by(Spike.D614G_tp2_Binom) %>%
  summarise(
    geometricMean = exp(mean(log(Spike.D614G_tp2))),
    ymin = exp(mean(log(Spike.D614G_tp2)) - sd(log(Spike.D614G_tp2))),
    ymax = exp(mean(log(Spike.D614G_tp2)) + sd(log(Spike.D614G_tp2)))
  )
geomMeans$labels <- paste("GMT =", round(geomMeans$geometricMean, 0))
geomMeans$yPos <- exp(8.8) + geomMeans$ymax

## Define color map
colors <- c("High" = colorHigh, "Low" = colorLow)

## Calculate FRNT50 threshold (40th percentile)
threshFRNT50 <- quantile(data$Spike.D614G_tp2, 0.4)


# =============================================================================
# SCATTER PLOT
# =============================================================================

p1 <- ggplot(data, aes(x = Spike.D614G_tp2_Binom, y = Spike.D614G_tp2, fill = Spike.D614G_tp2_Binom)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = threshFRNT50, ymax = Inf, alpha = .2, fill = lighten(colorHigh)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 190, ymax = threshFRNT50, alpha = .2, fill = lighten(colorLow)) +
  geom_quasirandom(shape = 21, size = 4, stroke = 0.5, alpha = 0.75, na.rm = TRUE) +
  labs(x = "Anti-spike Ab category", y = bquote(bold("Level"))) +
  scale_fill_manual(values = colors) +
  scale_colour_manual(values = "black") +
  scale_y_continuous(trans = "log10", labels = formatScientific) +
  geom_errorbar(
    data = geomMeans,
    aes(x = Spike.D614G_tp2_Binom, y = geometricMean, ymin = ymin, ymax = ymax),
    position = position_dodge(width = 0.2),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_segment(
    data = geomMeans,
    aes(x = c(0.95, 1.95), xend = c(1.05, 2.05), y = geometricMean, yend = geometricMean),
    color = "black", linewidth = 1.2
  ) +
  geom_text(aes(x = 1, y = 12000, label = geomMeans$labels[1]), size = 4) +
  geom_text(aes(x = 2, y = 900, label = geomMeans$labels[2]), size = 4) +
  themePublication()


# =============================================================================
# DENSITY PLOT
# =============================================================================

p2 <- ggplot(data, aes(x = Spike.D614G_tp2, fill = Spike.D614G_tp2_Binom)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(trans = "log10") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )

## Align density plot
p2Final <- ggdraw() + draw_plot(p2, -0.13, -0.02, 1, 1)


# =============================================================================
# FINAL PLOT
# =============================================================================

finalPlot <- plot_grid(p1, p2Final, labels = NULL, ncol = 2, rel_widths = c(4, 1), scale = c(1, 1), align = "h")
finalPlot


# =============================================================================
# SAVE PLOT
# =============================================================================

fileName <- paste0("../../figuresAndTables/fig4/fig4A_groupScatterAntiSpike.png")
widthPlot <- 6

ggsave(fileName, finalPlot, width = widthPlot, height = 0.8 * widthPlot, bg = "white", dpi = 600)
