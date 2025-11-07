library(dplyr)

# Load existing data
df <- read.csv("/home/afrahkhairallah/Documents/Alex/Paper/06052025/Merged_Metadata_with_Outpatient_Days.csv")
min_max_df <- read.csv("/home/afrahkhairallah/Documents/Alex/Paper/06052025/14052025/min_max_tooutpatient.csv")

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

table(data$NLR3_Binom, data$Neuts_Binom)

table(data$CD4LymphocytesAbs__min, data$Neuts_Binom)

table(data$Age_At_Enrolment.x, data$Neuts_Binom)

table(data$Lymph_abs_min, data$Neuts_Binom)

table(data$Comorbid_Binom, data$Neuts_Binom)

table(data$Sex_Bino1, data$Neuts_Binom)

table(data$HIV_Binom, data$Neuts_Binom)


## Select relevant variables and omit rows with NA (NLR_Binom has one NA)
df <- data %>%
  select(Neuts_Binom, Age_At_Enrolment.x,Sex_Bino1, NLR3_Binom,Comorbid_Binom, HIV_Binom,CD4LymphocytesAbs__min,Lymph_abs_min) %>%
  na.omit()
write.csv(df, "df.csv")  

  
library(ggplot2)
library(dplyr)
library(scales)

# Color palette from the provided image
pub_cols <- c(
  "High NLR + ≥1 Risk"           = "#E76F51",  # Dark teal
  "High NLR only"                = "#1497eb",  # Orange
  "Risk factors only"            = "#264653",  # Coral
  "Neither High NLR nor Risk"    = "#E9C46A"   # Yellow
)


df_low <- df %>% filter(Neuts_Binom == 0)
#df_high <- df %>% filter(Neuts_Binom == 1)

# Pie chart with black text labels
pie_nlr_risk <- df_low %>%
  mutate(
    HighNLR = NLR3_Binom == 1,
    Risk    = Comorbid_Binom == 1 | Age_At_Enrolment.x == 1 | Sex_Bino1 == 1,
    Group   = case_when(
      HighNLR &  Risk ~ "High NLR + ≥1 Risk",
      HighNLR & !Risk ~ "High NLR only",
      !HighNLR & Risk ~ "Risk factors only",
      TRUE            ~ "Neither High NLR nor Risk"
    )
  ) %>%
  count(Group) %>%
  mutate(
    prop = n / sum(n),
    label = paste0(n, " (", percent(prop, accuracy = 1), ")")
  ) %>%
  ggplot(aes(x = "", y = n, fill = Group)) +
    geom_col(width = 1, color = "black", size = 1.1) +  # Black borders added
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      color = "white",        # black text instead of white
      size = 4.5,
      fontface = "bold",
      lineheight = 0.9
    ) +
    scale_fill_manual(values = pub_cols) +
    labs(
      title = "High NLR & Risk-Factor Status",
      fill = "Group"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

# View the chart
pie_nlr_risk


## ------------------------------------------------------------------
## 3.  Save as high-resolution PNG (600 dpi) -------------------------
ggsave(
  filename = "HighNLR_RiskFactors_pie_Low_Neuts.png",
  plot     = pie_nlr_risk,
  width    = 6,          # inches
  height   = 6,
  dpi      = 600,
  bg       = "white"     # keeps background clean in PNGs
)

