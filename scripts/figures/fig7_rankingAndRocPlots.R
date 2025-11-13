# LOGISTIC REGRESSION WITH STEPWISE REGRESSION
## ALICE PILLER
## JULY 2024



# DESCRIPTION -------------------------------------------------------------

# This script constructs a logistic regression model to predict neutralization 
# response.

# Prerequisites:
# 1) Data must be split into training and testing data using train-testSplit.R 
# script.
# 2) Significant differentially regulated proteins must be determined using 
# training data and the differentialGeneExpressionForTraining.R script.

# The following is carried out on the training data:
# Neutralization response is iteratively fitted, using bootstrapping, against 
# the top significant proteins, as determined in 2), with forward and backward 
# step wise regression only selecting the best predictors that optimize Akaike 
# information criterion (AIC). The selected variables are added to a running 
# tally which serves to rank the proteins upon completion of the bootstrap 
# iterations.

# The top 3 proteins from the ranked list are used to construct a multivariate 
# logistic regression model, and univariate models are constructed from the 
# single proteins. The performance of the multivariate and 3 univariate models
# is assessed on the test data using ROC curves and AUC values.



# SETUP -------------------------------------------------------------------

## Remove variables in global environment
rm(list = ls())

## Install packages
# install.packages("MASS")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("pROC")
# install.packages("caret")
# install.packages("conflicted")
# install.packages("boot")

## Import packages
pacman::p_load(MASS, dplyr, ggplot2, pROC, caret, conflicted, boot, ggthemes)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

## Set working directory
<<<<<<< HEAD:scripts/figuresAndTables/fig7_rankingAndRocPlots.R
### Provide path to sourceDataAndFigures folder
path <- "C:/Users/Alice.Piller/OneDrive - AHRI/Documents/COVID-19 Neutralisation/"
setwd(paste0(path, "/sourceDataAndFigures_Jun25/data/predictiveModel"))

=======
###sourceDataAndFigures folder
setwd(paste0(path, "/data/predictiveModel"))
>>>>>>> e4eb264 (Update local repo: delete old files, add new ones):scripts/figures/fig7_rankingAndRocPlots.R


# PARAMETERS --------------------------------------------------------------

## Transformations and scaling
logTrans <- TRUE
scale <- TRUE

## Choose the number of resampling iterations for bootstrap stepwise regularisation
n_iterations_bs <- 1000



# FUNCTIONS ---------------------------------------------------------------

## Theme to generate Publication quality figures 
theme_Publication <- function(base_size = 17, base_family = "helvetica") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
      text = element_text(face = "bold"),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = 16),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(face = "bold", size = 15), 
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_line(colour = "#f0f0f0"),
      legend.key = element_rect(colour = "black", linewidth = 0.8),
      legend.position = "top",
      #legend.position = "none",
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "lines"),
      legend.margin = margin(0, 0, 0, 0, unit = "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      #legend.title = element_text(face = "italic"),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
}


## Produces a ROC curve for the input model for a given predictor
rocPlotUni <- function(model, predictor, col = "black"){
  # Predict response probabilities from test data using input model
  pred_probs <- predict.glm(model, type = "response", newdata = test)
  # Actual responses
  observed <- test$Neuts_Binom
  # Create ROC object
  roc_obj <- roc(observed, pred_probs) 
  # ROC curve
  roc_plot <- ggroc(roc_obj, size = 1.5, alpha = 1, legacy.axes = TRUE, colour = col) +
    geom_abline(intercept = 0, linetype = "dashed", color = "#6E788D", linewidth = 1) +
    labs(
      title = paste0(predictor, " (AUC = ", round(roc_obj$auc, 2), ")"),
    ) +
    theme_Publication()
  return(roc_plot)
}

## Conducts forward and backward stepwise regression iteratively over bootstrap samples and keeps a running tally of protein predictors selected in each model
bootstrapSelection <- function(data, n_iter = 1000) {
  # Empty list to store selected protein predictors
  selected_vars <- list()
  for (i in 1:n_iter) {
    # Sample 100x with replacement
    sample_index <- sample(1:nrow(data), size = 100, replace = TRUE)
    sample_data <- data[sample_index, ]
    # Conduct step wise regression on sample
    model <- stepAIC(glm(Neuts_Binom ~ ., data = sample_data, family = binomial()), direction = "both", trace = FALSE)
    # Add selected variables to list
    selected_vars[[i]] <- names(coef(model))[-1]
  }
  return(table(unlist(selected_vars)))
}



# DATA ---------------------------------------------------------------

## Read in training and testing data
train_full <- read.csv("training_0.6_190625.csv")
test_full <- read.csv("testing_0.6_190625.csv")

## Top differentially regulated proteins in training data
topProteins <- read.csv("./proteins/threshResultsTraining_190625.csv") %>%
  arrange(FDR)

## Create data frames for training and testing comprising the response (Neuts_Binom) and the predictors (topProteins)
train_subset <- train_full %>%
  select(Neuts_Binom, all_of(topProteins$ProteinID))
test_subset <- test_full %>%
  select(Neuts_Binom, all_of(topProteins$ProteinID))

## Log-transform genes
if (logTrans == TRUE){
  train_log <- train_subset %>%
    mutate(across(-1, ~ log(.))) # Log-transform all genes
  test_log <- test_subset %>%
    mutate(across(-1, ~ log(.)))
}

## Scale genes
### Scaling parameters for the test set are acquired from the training set
if(scale == TRUE){
  ## Scaling training set genes
  train_genes <- train_log[, -1] # minus Neuts_Binom
  train_genes_scaled <- scale(train_genes)
  # Extract mean and standard deviation from training data to scale testing data
  scale_centre <- attr(train_genes_scaled, "scaled:center")
  scale_scale <- attr(train_genes_scaled, "scaled:scale")
  # Add Neuts_Binom back
  train_scaled <- train_genes_scaled %>%
    as.data.frame() %>%
    mutate(Neuts_Binom = train_log$Neuts_Binom, .before = 1)
  
  ## Scaling testing set genes
  test_genes <- test_log[, -1] # minus Neuts_Binom
  test_genes_scaled <- scale(test_genes, center = scale_centre, scale = scale_scale)
  test_scaled <- test_genes_scaled %>%
    as.data.frame() %>%
    mutate(Neuts_Binom = test_log$Neuts_Binom, .before = 1)
  train <- train_scaled
  test <- test_scaled
}



# STEPWISE REGRESSION WITH BOOTSTRAPPING ----------------------------------

## Bootstrapping
## 1) Training data is sampled 100x with replacement
## 2) Logistic model is trained with stepwise regression, which optimises AIC
## 3) Remaining (top) predictors are tallied
## 4) Process is repeated
## 5) Output is a list of the 12 genes ranked by predictive power (according to training data)
set.seed(123)
bootstrapResults <- bootstrapSelection(train, n_iter = n_iterations_bs)
# 
## Rank genes from those that were most included in models from bootstrapping
rankedGenes <- sort(bootstrapResults, decreasing = T)
rankedGenes

## Results (so that you don't need to run the bootstrap each time)
# rankedGenes <- list(
#   MLN = 871,
#   FAP = 786,
#   HSPA8 = 674,
#   RTP4 = 555,
#   HP = 272,
#   RRBP1 = 251,
#   VWF = 137,
#   GBP1 = 120,
#   PLAT = 14,
#   VOPP1 = 13,
#   CDON = 10,
#   ADSSL1 = 8
# )



# TESTING -----------------------------------------------------------------

## ROC Curve

### Multivariate with top 3 proteins
modelMulti <- glm(Neuts_Binom ~ MLN + FAP + HSPA8, data = train, family = binomial, maxit = 50)
rocCurveMulti <- rocPlotUni(modelMulti, predictor = "Multivariate", col = "#A54891")
rocCurveMulti

## HSPA8
modelHSPA8 <- glm(Neuts_Binom ~ HSPA8, data = train, family = binomial, maxit = 50)
rocCurveHSPA8 <- rocPlotUni(modelHSPA8, predictor = "HSPA8", col = "#006EAE")
rocCurveHSPA8

## FAP
modelFAP <- glm(Neuts_Binom ~ FAP, data = train, family = binomial, maxit = 50)
rocCurveFAP <- rocPlotUni(modelFAP, predictor = "FAP", col = "#0096A0")
rocCurveFAP

## FAP
modelMLN <- glm(Neuts_Binom ~ MLN, data = train, family = binomial, maxit = 50)
rocCurveMLN <- rocPlotUni(modelMLN, predictor = "MLN", col = "#429130")
rocCurveMLN



# SAVE PLOTS ---------------------------------------------------------------

w = 6.5
ggsave("../../figuresAndTables/fig7/fig7A_rocCurveMulti.png", rocCurveMulti, width = w, height = w, bg = "white", dpi = 600)
ggsave("../../figuresAndTables/fig7/fig7B_rocCurveHSPA8.png", rocCurveHSPA8, width = w, height = w, bg = "white", dpi = 600)
ggsave("../../figuresAndTables/fig7/fig7C_rocCurveFAP.png", rocCurveFAP, width = w, height = w, bg = "white", dpi = 600)
ggsave("../../figuresAndTables/fig7/fig7D_rocCurveMLN.png", rocCurveMLN, width = w, height = w, bg = "white", dpi = 600)

