# =============================================================================
# GSEA SIGNIFICANT PATHWAYS PLOT
# Alice Piller
# May 2024
# =============================================================================


# =============================================================================
# DESCRIPTION
# =============================================================================

# This script generates a publication-quality plot of significantly upregulated
# and downregulated hallmark gene sets (FDR < 0.1) for a selected set of variables
# including neutralization titers, binned neutralization categories, and anti-spike
# antibody response. Shared and condition-specific pathways are annotated.


# =============================================================================
# SETUP
# =============================================================================

## Clear workspace
rm(list = ls())

## Load required libraries
pacman::p_load(ggplot2, dplyr, ggthemes, readxl, stringr, RColorBrewer, tidyr, ggpubr)


# =============================================================================
# PARAMETERS
# =============================================================================

## Conditions for GSEA comparison
conds <- c("D614GNeutralization", "Neuts_Binom", "Spike.D614G_tp2_Binom")

## Human-readable labels for each condition
labelsLookup <- list(
  Neuts_Binom = "Neut. cat.", 
  D614GNeutralization = "Neut. num.",
  Spike.D614G_tp2_Binom = "Anti-spike Ab"
)

## Custom colors for each condition
colorsLookup <- list(
  Neuts_Binom = "#7e45a4", 
  D614GNeutralization = "#be456e",
  Spike.D614G_tp2_Binom = "#018F99"
)


# =============================================================================
# FUNCTIONS
# =============================================================================

## Custom theme for consistent publication-quality plots
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
      legend.key.size = unit(0.2, "lines"),
      legend.margin = margin(0, 0, 0, 0, unit = "cm"),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
}


# =============================================================================
# READ DATA
# =============================================================================

## Set working directory to the location of the GSEA results
path <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
resultsPath <- paste0(path, "/sourceDataAndFigures_Jun25/data/gsea")
setwd(resultsPath)

## Load all relevant GSEA result files for selected conditions
files <- c()
for (cond in conds){
  files <- append(files, list.files(path = paste0("./", cond), pattern = paste0("_", cond, "[._]"), full.names = TRUE))
}

## Read upregulated and downregulated pathway sheets
pathways <- list()
for (file in files){
  up <- read_excel(file, sheet = "Upregulated")
  dn <- read_excel(file, sheet = "Downregulated")
  group <- str_match(basename(file), "gseaResults_(.*?)\\.xlsx$")[2]
  up <- up %>% mutate(Group = group)
  dn <- dn %>% mutate(Group = group)
  pathways[[length(pathways) + 1]] <- up
  pathways[[length(pathways) + 1]] <- dn
}

## Combine and preprocess pathway data
pathways <- bind_rows(pathways)
pathways$FDR.q.val <- ifelse(pathways$FDR.q.val == 0, 1e-10, pathways$FDR.q.val)

## Filter significant results and format pathway names
pathways <- pathways %>%
  filter(FDR.q.val < 0.1) %>%
  mutate(`-log10(FDR)` = -log10(FDR.q.val)) %>%
  rename(PATHWAY = NAME) %>%
  mutate(PATHWAY = str_replace_all(PATHWAY, "HALLMARK_", "")) %>%
  mutate(PATHWAY = str_replace_all(PATHWAY, "_", " ")) %>%
  mutate(PATHWAY = str_to_sentence(str_to_lower(PATHWAY))) %>%
  mutate(PATHWAY = case_when(
    PATHWAY == "Mtorc1 signaling" ~ "mTORC1 signaling",
    PATHWAY == "Pi3k akt mtor signaling" ~ "PI3k/AKT/mTOR signaling",
    PATHWAY == "Uv response up" ~ "UV response up",
    PATHWAY == "Dna repair" ~ "DNA repair",
    TRUE ~ PATHWAY
  )) %>%
  mutate(Direction = if_else(NES > 0, "Upregulated", "Downregulated")) %>%
  mutate(Direction = factor(Direction, levels = c("Upregulated", "Downregulated")))

## Determine pathway order by shared/condition-specific grouping
shared <- pathways %>%
  group_by(PATHWAY) %>%
  filter(n() > 1) %>%
  summarise(NES_average = mean(NES)) %>%
  arrange(desc(NES_average)) %>%
  pull(PATHWAY)

cond1 <- setdiff(pathways$PATHWAY[pathways$Group == conds[1]], shared)
cond2 <- setdiff(pathways$PATHWAY[pathways$Group == conds[2]], shared)
cond3 <- setdiff(pathways$PATHWAY[pathways$Group == conds[3]], shared)

pathways$PATHWAY <- factor(pathways$PATHWAY, levels = c(shared, cond1, cond2, cond3))

## Apply readable labels to group variable
pathways <- pathways %>%
  mutate(Group = case_when(
    Group == conds[1] ~ labelsLookup[[conds[1]]],
    Group == conds[2] ~ labelsLookup[[conds[2]]],
    Group == conds[3] ~ labelsLookup[[conds[3]]]
  ))


# =============================================================================
# PLOT
# =============================================================================

## Set up plot colors
colours <- unlist(colorsLookup[conds])
names(colours) <- unlist(labelsLookup[conds])

## Define bracket positions for pathway group annotation
x1a <- 0.6
x1b <- length(shared) + 0.4
x2a <- length(shared) + 0.6
x2b <- length(c(shared, cond1)) + 0.4
x3a <- length(c(shared, cond1)) + 0.6
x3b <- length(c(shared, cond1, cond2)) + 0.4

## Build plot
p <- ggplot(pathways, aes(x = PATHWAY, y = NES, size = `-log10(FDR)`)) +
  geom_segment(aes(x = PATHWAY, xend = PATHWAY, y = 0, yend = NES, alpha = `-log10(FDR)`), size = 0.8) +
  geom_point(aes(fill = Group), shape = 21, alpha = 0.5) +
  geom_point(aes(colour = Group), shape = 21, stroke = 1.5) +
  geom_bracket(xmin = x1a, xmax = x1b, y.position = 2.7, label = "Shared", tip.length = 0.01, coord.flip = TRUE, fontface = "bold") +
  geom_bracket(xmin = x2a, xmax = x2b, y.position = 2.7, label = labelsLookup[[conds[1]]], tip.length = 0.01, coord.flip = TRUE, fontface = "bold") +
  geom_bracket(xmin = x3a, xmax = x3b, y.position = 2.7, label = labelsLookup[[conds[2]]], tip.length = 0.01, coord.flip = TRUE, fontface = "bold") +
  labs(x = "Pathway",
       size = expression(bold(-log[10](FDR))),
       alpha = expression(bold(-log[10](FDR)))) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  scale_size_continuous(range = c(3, 12)) +
  coord_flip() +
  themePublication() +
  guides(fill = guide_legend(override.aes = list(size = 2.5)))

p


# =============================================================================
# SAVE PLOT
# =============================================================================

## Save final plot to file
fileName <- paste0("../../figuresAndTables/fig6/fig6_gsea", conds[1], "_", conds[2], "_", conds[3], ".png")
ggsave(fileName, plot = p, width = 10.5, height = 10.5, dpi = 600, bg = "white")
