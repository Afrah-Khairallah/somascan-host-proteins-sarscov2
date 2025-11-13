library(dplyr)

# Load existing data
df <- read.csv("/data/cleaned/Merged_Metadata_with_Outpatient_Days.csv")
min_max_df <- read.csv("/data/cleaned/min_max_tooutpatient.csv")

# Ensure PID columns are character
df$PID <- as.character(df$PID)
min_max_df$PID_ <- as.character(min_max_df$PID)

# View result
head(data)
######################################################
data <- df %>%
  left_join(min_max_df, by = c("PID" = "PID")) %>%
  filter(SuppO2_Binom == 0)  # Keep only rows where SuppO2_Binom is 0
##########################################################


threshold <- 6  # Set threshold for binary classification
data$NLR3_Binom <- ifelse(data$NLR_max > threshold, 1, 0)
table(data$NLR3_Binom)


median(data$CD4LymphocytesAbs__min, na.rm = TRUE)

#Low CD4 =  1
threshold2 <- 350  # Set threshold for binary classification
data$CD4LymphocytesAbs__min <- ifelse(data$CD4LymphocytesAbs__min > threshold2, 0, 1)
table(data$CD4LymphocytesAbs__min)

#Older than 50 = 1
threshold2 <- 50  # Set threshold for binary classification
data$Age_At_Enrolment.x <- ifelse(data$Age_At_Enrolment.x > threshold2, 1, 0)
table(data$Age_At_Enrolment.x)

#Low Lymphocyte = 1
threshold2 <- 1.1  # Set threshold for binary classification
data$Lymph_abs_min <- ifelse(data$Lymph_abs_min > threshold2, 0, 1)
table(data$Lymph_abs_min)

data$Comorbid_Binom <- ifelse(data$Diabetes_Binom == 1 | data$Hyp_Binom == 1, 1, 0)

> table(data$NLR3_Binom, data$Neuts_Binom)

> table(data$CD4LymphocytesAbs__min, data$Neuts_Binom)

> table(data$Age_At_Enrolment.x, data$Neuts_Binom)

> table(data$Lymph_abs_min, data$Neuts_Binom)

> table(data$Comorbid_Binom, data$Neuts_Binom)

> table(data$Sex_Bino1, data$Neuts_Binom)

> table(data$HIV_Binom, data$Neuts_Binom)



## Select relevant variables and omit rows with NA (NLR_Binom has one NA)
df <- data %>%
  select(Neuts_Binom, Age_At_Enrolment.x,Sex_Bino1, NLR3_Binom,Comorbid_Binom, HIV_Binom,CD4LymphocytesAbs__min,Lymph_abs_min) %>%
  na.omit()

## Import packages
pacman::p_load(dplyr, caret, ggplot2, ggthemes, conflicted)

## Specifying to use select from dplyr because of package clashes
conflict_prefer("select", "dplyr")


# FUNCTIONS ---------------------------------------------------------------
theme_Publication <- function(base_size = 19, base_family = "helvetica") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1), hjust = 0.5),
      text = element_text(face = "bold"),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = 18),
      axis.title.y = element_blank(),
      #axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(face = "bold"), 
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom", # Change this to "right" if preferred
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "cm"),
      legend.title = element_blank(),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
}


# Converts scientific to notation from 1e02 format to 10^2
scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e\\+", "e", l)
  # turn the 'e+' into plotmatHigh format
  l <- gsub("e", "10^", l)
  # return this as an expression
  parse(text=l)
}

## Extracts the odds ratio, confidence interval and p-value for a given predictor from the multivariate model
get_summary_values <- function(model, predictor){
  # Coefficients from input model
  coeff <- summary(model)$coefficients
  # Odds ratio
  estimate <- coeff[, "Estimate"][predictor]
  OR <- exp(estimate)
  # Confidence intervals
  log_CI <- confint(model)[predictor,]
  CI <- exp(log_CI)
  # P-value
  P_val <- coeff[, "Pr(>|z|)"][predictor]
  return(list(model = model,
              estimate = estimate,
              OR = OR,
              CI = CI,
              P_val = P_val))
}

## Creates a summary row for a given predictor formatted as such: OR (95% CI)   P-value
## Takes output of get_summary_values() and formats for export as a Word table
summary_row_table <- function(predictor, results){
  OR <- signif(results[[predictor]]$OR, 2)
  CI <- signif(results[[predictor]]$CI, 2)
  P_val <- signif(results[[predictor]]$P_val, 2)
  # OR and CI String
  OR_CI_str <- paste0(OR, " (", CI[1], "; ", CI[2], ")")
  row <- data.frame(OR_95CI = OR_CI_str,
                    p_value = P_val)
  return(row)
}

## Creates a summary row for a given predictor to be used in the odds ratio plot
## Columns are Condition, Odds ratio, CI 2.5, CI 97.5, and the p-value of the estimate
## Summary rows of all univariates will be compiled into single data frame used for plotting
summary_row_plot <- function(predictor, results){
  OR <- results[[predictor]]$OR
  CI <- results[[predictor]]$CI
  P_val <- results[[predictor]]$P_val
  row <- data.frame(Condition = predictor, OR = OR, CI2.5 = CI[1], CI97.5 = CI[2], P_value = P_val)
  return(row)
}


# CREATE UNIVARIATE MODELS -------------------------------------------------

## List of explanatory variables
predictors <- setdiff(names(df), "Neuts_Binom")

## Create empty list and data frames to store results from univariate analyses
univariate_results <- list()
univariate_results_table <- data.frame(OR95CI = c(), P_value = c())
univariate_results_plot_df <- data.frame(Condition = c(), OR = c(), CI2.5 = c(), CI97.5 = c(), P_value = c())

## Conduct univariate analyses and store results
for(predictor in predictors){
  # Train univariate model
  formula <- as.formula(paste("Neuts_Binom ~", predictor))
  model <- glm(formula, data = df, family = binomial)
  # Get odds ratio, confidence interval and p-value of predictor
  univariate_results[[predictor]] <- get_summary_values(model, predictor)
  # Table: summarise results in table to be exported as Word table
  univariate_results_table <- rbind(univariate_results_table, summary_row_table(predictor, univariate_results))
  # Data frame: summarise results in data frame for plotting
  univariate_results_plot_df <- rbind(univariate_results_plot_df, summary_row_plot(predictor, univariate_results))
}


# MULTIVARIATE MODEL ------------------------------------------------------

## Construct multivariate formula with all factors
formula <- as.formula(paste("Neuts_Binom ~", paste(predictors, collapse = " + ")))

## Multivariate model
multi_model <- glm(formula, data = df, family = binomial)

## Create empty list and data frames to store results from multivariate analysis
multivariate_results <- list()
multivariate_results_table <- data.frame(`OR (95% CI)` = c(), `P-value` = c())
multivariate_results_plot_df <- data.frame(Condition = c(), OR = c(), CI2.5 = c(), CI97.5 = c(), P_value = c())

## Conduct multivariate analysis and store results
for(predictor in predictors){
  # Get odds ratio, confidence interval and p-value of predictor
  multivariate_results[[predictor]] <- get_summary_values(multi_model, predictor)
  # Table: summarise results in table to be exported as Word table
  multivariate_results_table <- rbind(multivariate_results_table, summary_row_table(predictor, multivariate_results))
  # Data frame: summarise results in data frame for plotting
  multivariate_results_plot_df <- rbind(multivariate_results_plot_df, summary_row_plot(predictor, multivariate_results))
}


# WRITE SUMMARY CSV TABLES ---------------------------------------------------------

## Write to csv
# write.csv(univariate_results_table, file = "../../figuresAndTables/supplementary/univariateResults.csv")
# write.csv(multivariate_results_table, file = "../../figuresAndTables/supplementary/multivariateResults.csv")
# 

# PLOT SETUP --------------------------------------------------------------------

## Updated list of variable names and their plot labels (unchanged)
variable_to_label <- list(
  Age_At_Enrolment.x = "Age > 50", 
  Sex_Bino1 = "Male",
  NLR3_Binom = "High NLR", 
  Comorbid_Binom ="Comorbidities",
  HIV_Binom = "PLWH",
  CD4LymphocytesAbs__min = "Low CD4",
  Lymph_abs_min = "Lymphopenia"
)

## Replace variable names with display labels
univariate_results_plot_df$Condition <- sapply(univariate_results_plot_df$Condition, function(x) variable_to_label[[x]])
multivariate_results_plot_df$Condition <- sapply(multivariate_results_plot_df$Condition, function(x) variable_to_label[[x]])

## Reverse the order of labels
ordered_labels <- rev(unname(unlist(variable_to_label)))

univariate_results_plot_df$Condition <- factor(univariate_results_plot_df$Condition, levels = ordered_labels)
multivariate_results_plot_df$Condition <- factor(multivariate_results_plot_df$Condition, levels = ordered_labels)

## Custom labels in reverse order for plotting (bold)
custom_labels <- rev(c(
  bquote(bold("Age > 50")),
  bquote(bold("Male")),
  bquote(bold("High NLR")),
  bquote(bold("Comorbidities")),
  bquote(bold("PLWH")),
  bquote(bold("Low CD4")),
  bquote(bold("Lymphopenia"))
))

## Assign colors according to reversed order
colours <- c(
  "#C5373D", "#A54891", "#006EAE", "#0096A0", "#429130", "#96A00A", "#E9C54E", "#F29742"
)
names(colours) <- ordered_labels



# PLOT --------------------------------------------------------------------


## Plot univariate odds ratios
p_uni <- ggplot(univariate_results_plot_df, aes(x = Condition, y = OR)) +
  # Title
  ggtitle("Univariate") +
  # Confidence interval bars
  geom_errorbar(
    aes(x = Condition, y = OR, ymin = CI2.5, ymax = CI97.5),
    width = 0.75,
    linewidth = 0.8
  ) +
  # Point outline
  geom_point(aes(colour = Condition), size = 5, shape = 21, stroke = 1.5) +
  # Point fill
  geom_point(aes(fill = Condition, colour = Condition), size = 5, shape = 21, alpha = 0.5) +
  # p-values
  geom_text(aes(label = paste("p =", signif(P_value, 2))), vjust = -1.1) + 
  # Horizontal line at OR = 1
  geom_hline(yintercept = 1, colour = "#96A0B3", size = 1, linetype = "dashed") +
  # Adding custom labels
  scale_x_discrete(labels = custom_labels) +
  # Log-scaling y-axis
  scale_y_continuous(transform = "log10") +
  # Colours
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  scale_shape_binned(16) +
  # Flip x and y axes
  coord_flip() +
  # Theme
  guides(fill = guide_legend(override.aes = list(size = 2.5))) +
  theme_Publication()

p_uni

## Multivariate plot
p_multi <- ggplot(multivariate_results_plot_df, aes(x = Condition, y = OR)) +
  # Title
  ggtitle("Multivariate") +
  # Confidence interval bars
  geom_errorbar(
    aes(x = Condition, y = OR, ymin = CI2.5, ymax = CI97.5),
    width = 0.75,
    linewidth = 0.8
  ) +
  # Point outline
  geom_point(aes(colour = Condition), size = 5, shape = 21, stroke = 1.5) +
  # Point fill
  geom_point(aes(fill = Condition, colour = Condition), size = 5, shape = 21, alpha = 0.5) +
  # p-values
  geom_text(aes(label = paste("p =", signif(P_value, 2))), vjust = -1.1) +
  # Horizontal line at OR = 1
  geom_hline(yintercept = 1, colour = "#96A0B3", size = 1, linetype = "dashed") +
  # Adding custom labels
  scale_x_discrete(labels = custom_labels) +
  # Log-scaling y-axis
   scale_y_continuous(
    trans = "log10",
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  # Colours
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  scale_shape_binned(16) +
  # Flip x and y axes
  coord_flip() +
  # Theme
  guides(fill = guide_legend(override.aes = list(size = 2.5))) +
  theme_Publication()

p_multi


## Save plot
ggsave("uniOddsRatios_categorical_non_severe.png", plot = p_uni, width = 8.1, height = 7, dpi = 600)
ggsave("multiOddsRatioscategorical_non_severe_non_Sci.png", plot = p_multi, width = 8.1, height = 7, dpi = 600)
