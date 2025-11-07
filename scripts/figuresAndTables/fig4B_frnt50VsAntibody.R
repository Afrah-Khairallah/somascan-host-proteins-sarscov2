# =============================================================================
# CORRELATE FRNT50 AND ANTIBODY EC50
# =============================================================================
# Alice Piller
# June 2025


# =============================================================================
# DESCRIPTION
# =============================================================================

# This script generates a scatter plot and log-log linear regression line
# to visualize the relationship between FRNT50 neutralization titers and 
# antibody levels (e.g., anti-Spike, anti-RBD, or anti-Nucleocapsid). 
# The user can optionally remove values at the lower limit of detection (LOD).


# =============================================================================
# SETUP
# =============================================================================

## Clear environment
rm(list = ls())

## Set working directory
folder <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(folder, "/sourceDataAndFigures_Jun25/data/cleaned"))

## Load packages
pacman::p_load(dplyr, ggplot2, nls2, ggthemes, boot)


# =============================================================================
# PARAMETERS
# =============================================================================

condition <- "Spike.D614G_tp2"  # Antibody column to correlate with FRNT50

date <- "170625"

removeLodObs <- 0  # Set to 1 to remove LOD values

labelLookup <- list(
  RBD.D614G_tp2 = "Anti-RBD antibody",
  Spike.D614G_tp2 = "Anti-spike antibody",
  Nucleocapsid_tp2 = "Nucleocapsid antibody"
)

colorLookup <- list(
  RBD.D614G_tp2 = "#018F99",
  Spike.D614G_tp2 = "#3D892E",
  Nucleocapsid_tp2 = "#96a00a"
)


# =============================================================================
# FUNCTIONS
# =============================================================================

## Plotting theme for publication-quality visuals
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
      axis.text = element_text(face = "bold"), 
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_line(colour = "#f0f0f0"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.key.size = unit(1.1, "lines"),
      legend.margin = margin(0, 0, 0, 0, unit = "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
}

## Format scientific notation for axis ticks
scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e[+]?(-?\\d+)$", "10^\\2", l)
  parse(text = l)
}


# =============================================================================
# READ DATA
# =============================================================================

data <- read.csv("cleanedData_150625.csv")


# =============================================================================
# MODEL PREPARATION
# =============================================================================

modelDf <- data %>%
  mutate(abValue = .data[[condition]]) %>%
  select(D614GNeutralization, abValue) %>%
  mutate(
    logNeutralization = log(D614GNeutralization),
    logAbValue = log(abValue)
  ) %>%
  na.omit()

## Optionally remove observations at limit of detection
if (removeLodObs == 1) {
  modelDf <- modelDf %>% filter(abValue != 200)
}


# =============================================================================
# PLOT
# =============================================================================

scatterPlot <- ggplot(modelDf, aes(x = abValue, y = D614GNeutralization)) +
  geom_point(fill = colorLookup[[condition]], shape = 21, size = 3, stroke = 0.5, alpha = 0.75) +
  labs(x = labelLookup[[condition]], y = bquote(bold("FRNT"[50]*""))) +
  scale_colour_manual(values = "black") +
  scale_y_continuous(trans = "log10", labels = scientific, limits = c(1, 10000)) +
  scale_x_continuous(trans = "log10") +
  themePublication()


# =============================================================================
# LINEAR MODEL
# =============================================================================

linMod <- lm(logNeutralization ~ logAbValue, data = modelDf)

summaryStats <- summary(linMod)
adjR2 <- summaryStats$adj.r.squared
pVal <- signif(summaryStats$coefficients[2, 4], 2)

xNew <- data.frame(logAbValue = seq(min(modelDf$logAbValue), max(modelDf$logAbValue), length.out = 100))
predictions <- exp(predict(linMod, newdata = xNew))
fitLine <- data.frame(x = exp(xNew$logAbValue), y = predictions)

finalPlot <- scatterPlot +
  geom_line(data = fitLine, aes(x = x, y = y), color = "red", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("R² = ", round(adjR2, 3), "; p = ", pVal),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black")

finalPlot


# =============================================================================
# SAVE PLOT
# =============================================================================

## Define file name and save
fileName <- if (removeLodObs == 1) {
  paste0("../../figuresAndTables/fig4/fig4B_frnt50Vs", condition, "_no_lods_", date, ".png")
} else {
  paste0("../../figuresAndTables/fig4/fig4B_frnt50Vs", condition, "_", date, ".png")
}

w <- 5

ggsave(fileName, finalPlot, width = w, height = 0.88 * w, bg = "white", dpi = 600)
