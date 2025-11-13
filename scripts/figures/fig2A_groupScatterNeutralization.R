# =============================================================================
# GROUPED SCATTER PLOT - NEUTRALIZATION
# Author: Alice Piller
# Date: June 2025
# =============================================================================

# =============================================================================
# DESCRIPTION
# -----------------------------------------------------------------------------
# This script generates a grouped scatter plot of SARS-CoV-2 D614G FRNT50 
# neutralization values, stratified by binary neutralization response ("High" 
# vs "Low") as defined by the median FRNT50 value.
#
# The plot includes:
#  - Shaded areas indicating high vs low neutralization thresholds
#  - Scatter points for individual FRNT50 values
#  - Overlayed geometric mean titers (GMTs) and standard deviation intervals
#  - A horizontal dashed line for the limit of detection (FRNT50 = 25)
#  - An adjacent density plot to show distribution of values
# 
# The script uses custom publication-quality themes and color schemes. 
# Output: A high-resolution PNG image of the combined plot is saved to the 
# 'figuresAndTables' directory.
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
# install.packages("ggsignif")
# install.packages("cowplot")
# install.packages("stringr")
# install.packages("ggthemes")
# install.packages("colorspace")
# install.packages("ggbeeswarm")

## Load required packages
pacman::p_load(dplyr, ggplot2, ggsignif, cowplot, stringr, ggthemes, colorspace, ggbeeswarm)

## Set working directory to cleaned data folder
folder <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(folder, "/sourceDataAndFigures_Jun25/data/cleaned"))


# =============================================================================
# PARAMETERS
# =============================================================================

## Define colour palette for neutralization groups
col_1 <- "#CA9B23"  # High neutralization (gold)
col_0 <- "#006EAE"  # Low neutralization (blue)


# =============================================================================
# FUNCTIONS
# =============================================================================

## Custom publication-style theme for ggplot2
theme_Publication <- function(base_size = 17, base_family = "helvetica") {
  theme_foundation(base_size = base_size, base_family = base_family) +
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

## Format axis labels as base-10 scientific notation
scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e\\+", "e", l)
  l <- gsub("e", "10^", l)
  parse(text = l)
}


# =============================================================================
# READ DATA
# =============================================================================

## Load cleaned metadata and protein data
df <- read.csv("cleanedData_150625.csv")


# =============================================================================
# DATA PREPARATION
# =============================================================================

## Convert neutralization group from numeric to label
df$Neuts_Binom <- ifelse(df$Neuts_Binom == 1, "High", "Low")


# =============================================================================
# PLOT SETUP
# =============================================================================

## Compute geometric mean titers (GMT) and standard deviation intervals
geom_means <- df %>%
  group_by(Neuts_Binom) %>%
  summarise(
    geometric_mean = exp(mean(log(D614GNeutralization))),
    ymin = exp(mean(log(D614GNeutralization)) - sd(log(D614GNeutralization))),
    ymax = exp(mean(log(D614GNeutralization)) + sd(log(D614GNeutralization)))
  )

## Prepare GMT labels and y-axis placement
geom_means$labels <- paste("GMT =", round(geom_means$geometric_mean, 0))
geom_means$y_pos <- exp(8.8) + geom_means$ymax

## Define colour mapping
colours <- c("High" = col_1, "Low" = col_0)

## Compute median FRNT50 value
frnt50_median <- median(df$D614GNeutralization)


# =============================================================================
# SCATTER PLOT
# =============================================================================

p1 <- ggplot(df, aes(x = Neuts_Binom, y = D614GNeutralization, fill = Neuts_Binom)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = frnt50_median, ymax = Inf, alpha = .2, fill = lighten(col_1)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 1, ymax = frnt50_median, alpha = .2, fill = lighten(col_0)) +
  geom_quasirandom(shape = 21, size = 4, stroke = 0.5, alpha = 0.75, na.rm = TRUE) +
  labs(x = "Neutralization category",
       y = bquote(bold("FRNT"[50]*""))) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = "black") +
  scale_y_continuous(
    trans = "log10",
    labels = scientific,
    limits = c(1, 10000)
  ) +
  geom_errorbar(
    data = geom_means,
    aes(x = Neuts_Binom, y = geometric_mean, ymin = ymin, ymax = ymax),
    position = position_dodge(width = 0.2),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_segment(
    data = geom_means,
    aes(x = c(0.95, 1.95), xend = c(1.05, 2.05), y = geometric_mean, yend = geometric_mean),
    color = "black", linewidth = 1.2
  ) +
  geom_text(aes(x = 1, y = 10000, label = geom_means$labels[1]), size = 4) +
  geom_text(aes(x = 2, y = 700, label = geom_means$labels[2]), size = 4) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "red") +
  theme_Publication()

p1


# =============================================================================
# DENSITY PLOT
# =============================================================================

p2 <- ggplot(df, aes(x = D614GNeutralization, fill = Neuts_Binom)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = colours) +
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

## Align the density plot for combination with the main plot
p2_final <- ggdraw() + draw_plot(p2, -0.13, -0.02, 1, 1)

p2_final


# =============================================================================
# FINAL PLOT
# =============================================================================

## Combine scatter and density plots into a final figure
p_final <- plot_grid(p1, p2_final, labels = NULL, ncol = 2, rel_widths = c(4, 1), scale = c(1, 1), align = "h")

p_final


# =============================================================================
# SAVE PLOT
# =============================================================================

## Set output file name and save final image
filename <- paste0("../../figuresAndTables/fig2/fig2A_groupScatterNeutralization.png")
w <- 6
ggsave(filename, p_final, width = w, height = 0.8 * w, bg = "white", dpi = 600)
