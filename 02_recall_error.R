#---------------- About ------------------
# Project: Recall error in the Patient-Oriented Eczema Measure (POEM) 

# This script compares d-POEM (derived from daily diaries) and r-POEM (recalled) scores, and calculates
# recall bias and recall noise in the score.


# Initialisation ----------------------------------------------------------

library(tidyverse)
library(gmodels)
library(reshape2)
library(Metrics)
library(nlme)
library(cowplot)
library(gridExtra)
library(forcats)

source("functions.R")

data_dir <- '/homedir/'
dir <- here::here()

q_labels <- c("bleeding", "cracked", "dry", "flaky", "itchy", "oozing", "sleep")

missing_values_analysis <- TRUE # Determines if daily_df is the original dataframe or daily_imputed_df with missing daily pres/abs imputed.



# Load in dataset ---------------------------------------------------------
# ---- A. Daily dataframe
if(missing_values_analysis){
  load(file = paste0(dir,"/daily_imputed_df.RData"))
  daily_df <- daily_imputed_df
  rm(daily_imputed_df)
} else{
  daily_df <- load_daily(remove_NA = TRUE)
}

# ---- B. POEM dataframe
POEM_full_df <- load_POEM()

# For the purpose of our study, we cannot calculate a POEM score for POEM scores recorded at days < 7 of the study (as we do not have 7 days of daily scores to calculate a POEM score).
POEM_df <- POEM_full_df %>%
  filter(Analysis.Relative.Day > 6)





# Calculate POEM from daily symptom scores --------------------------------

# Calculate POEM score from the scores of the past 7 days' symptoms with a set threshold.
thr <- ifelse(missing_values_analysis, NA, 2) # daily scores in missing values analysis have already been thresholded
observed_calculated_POEM_df <- generate_observed_calculated(threshold = thr, daily_df = daily_df, POEM_df = POEM_df)




# Exploratory plots -------------------------------------------------------

# Visualizing frequency of daily symptoms in full dataset.
pl1 <- visualize_freq(df = daily_df, to_plot = "daily")


# Visualizing frequency of observed POEM (r-POEM) item scores in full dataset.
pl2 <- visualize_freq(df = POEM_full_df, to_plot = "POEM")


# Visualizing frequency of observed POEM (r-POEM) item scores in the 202 observations or 777 (for imputed dataset) where we have daily symptoms for the week before.
pl3 <- visualize_freq(df = observed_calculated_POEM_df, to_plot = "observed_POEM")


# Visualizing frequency of calculated POEM (d-POEM) item scores in the 202 observations or 777 (for imputed dataset) where we have daily symptoms for the week before.
pl4 <- visualize_freq(df = observed_calculated_POEM_df, to_plot = "calculated_POEM")


# Visualizing difference in observed (r-) and calculated (d-) POEM scores.
bw <- function(b, x) {b / bw.nrd0(x)}

df <- observed_calculated_POEM_df %>% 
  filter(Item == "total") %>%
  mutate(err = observed_POEM - calculated_POEM)

pl5 <- df %>%
  ggplot(aes(x = err)) + 
  geom_density(adjust = bw(1.8, df$err)) +
  geom_vline(aes(xintercept = 0), # plot x = 0
             color="black", size=1) +
  geom_vline(aes(xintercept = mean(err)), # plot bias
             color="blue", linetype="dashed", size=1)  +
  labs(x = "Difference between r-POEM and d-POEM", y = "Density") +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor.x = element_blank())




# Varying thresholds to calculate POEM from daily -------------------------

if(missing_values_analysis == FALSE){
  threshold_low = 1
  threshold_high = 10
  thresholds = seq(threshold_low, threshold_high)
  observed_calculated_POEM_thresholds <- generate_varying_thresholds(threshold_low = threshold_low, threshold_high = threshold_high, daily_df = daily_df, POEM_df = POEM_df)
  
  # Plot root mean square error for each threshold
  # A. For total scores
  observed_calculated_POEM_thresholds %>% filter(Item == "total") %>%
    group_by(threshold) %>%
    summarise(RMSE = rmse(observed_POEM, calculated_POEM)) %>%
    ungroup() %>%
    ggplot(aes(y = RMSE, x = threshold)) + 
    geom_point() + 
    scale_x_continuous(labels = as.character(thresholds), breaks = thresholds) + 
    scale_y_continuous(labels = as.character(seq(0,10,2)), breaks = seq(0,10,2), limits = c(0,10.5))+
    labs(x = "Threshold", y = "RMSE") +
    theme_bw(base_size = 15) +
    theme(panel.grid.minor.x = element_blank())
  
  # B. For individual items
  observed_calculated_POEM_thresholds %>% filter(Item != "total") %>%
    group_by(threshold, Item) %>%
    summarise(RMSE = rmse(observed_POEM, calculated_POEM)) %>%
    ungroup() %>%
    ggplot(aes(y = RMSE, x = threshold, color = fct_reorder2(Item, threshold, RMSE), group = fct_reorder2(Item, threshold, RMSE))) + 
    geom_point() + 
    geom_line() + 
    scale_x_continuous(labels = thresholds, breaks = thresholds) + 
    labs(x = "Threshold", y = "RMSE", color = "") +
    scale_colour_viridis_d() +
    scale_fill_viridis_d() +
    theme_bw(base_size = 15) +
    theme(panel.grid.minor.x = element_blank())
}




# Quantifying recall bias and noise in POEM measurements ------------------

# Bias is the average difference between r-POEM and d-POEM and noise is the SD.
recall_error_df <- observed_calculated_POEM_df %>%
  group_by(Item) %>%
  summarise(avg = mean(observed_POEM - calculated_POEM),
            sd = sd(observed_POEM - calculated_POEM)) %>%
  ungroup() 

# Rename POEM row.
recall_error_df[recall_error_df == "total"] <- "POEM"


# Plot recall bias
pl6 <- recall_error_df %>%
  ggplot(aes(y = Item, x = avg)) +
  geom_point(size = 4, color = "#03396c") +
  geom_point(size = 3.3, color = "#d1e1ec") +
  geom_vline(xintercept = 0, colour = "black") + #colour = "darkgreen"
  labs(x = "Recall bias (score)", y = "", colour = "") +
  theme_bw(base_size = 15) + scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +
  scale_y_discrete(limits = rev(c("POEM", q_labels))) +
  annotate("text",  x = -Inf, y = -Inf, label = "report lower score", vjust = -0.3, hjust = -0.1, colour = "black") + #colour = "darkgreen"
  annotate("text",  x = Inf, y = -Inf, label = "report higher score", vjust = -0.3, hjust = 1.1, colour = "black") #colour = "darkgreen"


# Plot recall noise
pl7 <- recall_error_df %>%
  ggplot(aes(y = Item, x = sd)) +
  geom_point(size = 4, color = "#03396c") +
  geom_point(size = 3.3, color = "#d1e1ec") +
  labs(x = "Recall noise (score)", y = "", colour = "") + # b
  theme_bw(base_size = 15) + scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
  scale_y_discrete(limits = rev(c("POEM", q_labels)))
