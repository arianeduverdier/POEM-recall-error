#---------------- About ------------------
# Project: Recall error in the Patient-Oriented Eczema Measure (POEM) 

# This script fits the recalled days model on the data and plots the results of the model.


# Initialisation ----------------------------------------------------------

library(tidyverse)
library(here)
library(rstan)
library(bayesplot)
library(Metrics)
options(mc.cores = parallel::detectCores()) # Parallel computing
rstan_options(auto_write=TRUE) 

set.seed(1744834965) # Reproducibility

source("functions.R")

plot_dir <- '/homedir/'
mdl_name <- "recalled_days_model"
stan_mdl <- paste0(here(),"/Model/", mdl_name, ".stan")
fit_dir <- here("results/Recall error model")

n_chains <- 4
n_it <- 4000

q_labels <- c("bleeding", "cracked", "dry", "flaky", "itchy", "oozing", "sleep")

missing_values_analysis <- TRUE # Determines if daily_df is the original dataframe or daily_imputed_df with missing daily pres/abs imputed.




# Load in dataset ---------------------------------------------------------

# ---- A. Daily dataframe
if(missing_values_analysis){
  load(file = here("daily_imputed_df.RData"))
  daily_df <- daily_imputed_df
  rm(daily_imputed_df)
} else{
  daily_df <- load_daily(remove_NA = TRUE)
}

# ---- B. POEM dataframe
# For the purpose of our study, we cannot calculate a POEM score for POEM scores recorded at days < 7 of the study (as we do not have 7 days of daily scores to 
# calculate a POEM score even with imputation) so we can remove those POEM scores.
POEM_df <- load_POEM() %>%
  filter(Analysis.Relative.Day > 6)




# Calculate number of days in week ----------------------------------------

threshold <- ifelse(missing_values_analysis, NA, 2)

# Calculate for each severity item (Q1 to Q7), for all instances we can, the number of days
# in the week represented by POEM with present symptoms from the daily symptom scores.

# 1. Calculate POEM score from the scores of the past 7 days' symptoms with a set threshold.
observed_calculated_POEM_df <- generate_observed_calculated(threshold = threshold, daily_df = daily_df, POEM_df = POEM_df) %>%
  filter(Item != "total")

# 2. Calculate number of days using columns of the observed_calculated_POEM_df.
if(missing_values_analysis){
  # Columns of observed_calculated_POEM_df have day1 to day7_observed pres/abs scores.
  observed_calculated_POEM_df$NumDays <- rowSums(observed_calculated_POEM_df %>%
                                                   subset(select = paste0("day", seq(1,7), "_observed"))) 
} else{
  # Columns of the observed_calculated_POEM_df have day1 to day7_observed daily scores.
  observed_calculated_POEM_df$NumDays <- rowSums(observed_calculated_POEM_df %>%
                                                   subset(select = all_of(paste0("day", seq(1,7), "_observed"))) >= threshold) 
}

# Rename patients
observed_calculated_POEM_df <- observed_calculated_POEM_df %>%
  group_by(Patient) %>%
  dplyr::mutate(Patient = cur_group_id()) %>%
  ungroup()




# Select 3 example patients --------------------------------------------------

# Select 3 random patients from the dataset with mild, moderate, and severe symptoms and model patients according to their intensity patterns. 
# Severity bandings according to: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3920642/
# 0–2 (clear/almost clear); 3–7 (mild); 8–16 (moderate); 17–24 (severe); 25–28 (very severe)
total_calculated_POEM_df  <- observed_calculated_POEM_df %>% group_by(Patient, Day) %>%
  summarise(total_d_POEM = sum(calculated_POEM)) %>%
  ungroup() %>%
  mutate(severity_POEM = ifelse(total_d_POEM %in% seq(3, 7), "mild", ifelse(total_d_POEM %in% seq(8, 16), "moderate", ifelse(total_d_POEM %in% seq(17, 24), "severe", "NA"))))
  
# Randomly select a mild patient
random_mild_pt <- total_calculated_POEM_df %>% filter(severity_POEM == "mild") %>%
  dplyr::sample_n(1)
  
# Randomly select a moderate patient
random_moderate_pt <- total_calculated_POEM_df %>% filter(severity_POEM == "moderate") %>%
  dplyr::sample_n(1)
  
# Randomly select a severe patient
random_severe_pt <- total_calculated_POEM_df %>% filter(severity_POEM == "severe") %>%
  dplyr::sample_n(1)
  
mild_pt <- observed_calculated_POEM_df %>% filter(Patient == random_mild_pt$Patient, Day == random_mild_pt$Day) 
moderate_pt <- observed_calculated_POEM_df %>% filter(Patient == random_moderate_pt$Patient, Day == random_moderate_pt$Day) 
severe_pt <- observed_calculated_POEM_df %>% filter(Patient == random_severe_pt$Patient, Day == random_severe_pt$Day) 




# FOR EACH ITEM --------------------------------------------------------------------------------------
# Fit model for one severity item at a time.

item <- q_labels[1] 

# Prepare stan data -------------------------------------------------------

df <- observed_calculated_POEM_df %>% filter(Item == item) 
# observed_calculated_POEM_df$observed_POEM is the recalled POEM score from 0-4, call this Qi where i = 1 to 7 for the 7 items
# observed_calculated_POEM_df$NumDays is the number of days in week with symptoms (calculated from daily scores)

mild_pt_days <- (mild_pt %>% filter(Item == item) %>% select(NumDays))[[1]]
moderate_pt_days <- (moderate_pt %>% filter(Item == item) %>% select(NumDays))[[1]]
severe_pt_days <- (severe_pt %>% filter(Item == item) %>% select(NumDays))[[1]]
cal_days_fake = c(mild_pt_days, moderate_pt_days, severe_pt_days)

format_stan_data <- function(df) {
  with(df,
       list(N = length(observed_POEM),
            calculated_days = NumDays, # Number of days in week with symptoms (calculated from daily scores.
            observed_POEM = observed_POEM, # Recalled POEM score from 0-4.
            N_fake = length(cal_days_fake), # Number of example patients.
            day_fake = cal_days_fake, # Calculated days for each of the example patients.
            run = 1
       ))
}
data_stan <- format_stan_data(df)


# Run fit model -----------------------------------------------------------

fit <- stan(file = stan_mdl,
                    data = data_stan,
                    iter = n_it, 
                    chains = n_chains, 
                    seed = seed,
                    control = list(adapt_delta = 0.99))

saveRDS(fit, file = paste0(fit_dir, "/fit_", mdl_name, "_", item, ".rds"))


# Check fit model  --------------------------------------------------------

fit <- readRDS(file = paste0(fit_dir, "/fit_", mdl_name, "_", item, ".rds"))

# - Diagnostics
check_hmc_diagnostics(fit)

neff = summary(fit)$summary[, "n_eff"]
all((neff/4e3) > 0.01)

Rhat = summary(fit)$summary[, "Rhat"]
all(Rhat < 1.1)


# - Plots of parameter values
pairs(fit, pars = c("b", "sigma"))

plot(fit, pars = c('b', 'sigma'))

plot(fit, pars = c('b', 'sigma'), plotfun = "trace")


# - Posterior predictive checking
recalled_days <- rstan::extract(fit, pars = "recalled_days")[[1]]

# True calculated days vs. latent recalled days.
ppc_dens_overlay(df$NumDays, recalled_days[1:500, ])

# Convert recalled days to POEM score.
POEM_dict <- data.frame(days = seq(0, 7), scores=c(0, 1, 1, 2, 2, 3, 3, 4))
observed_POEM_rep <- plyr::mapvalues(recalled_days, from = POEM_dict$days, to = POEM_dict$scores)
# Observed POEM score 0-4 vs. POEM score 0-4 from replicated latent recalled days.
ppc_dens_overlay(df$observed_POEM, observed_POEM_rep[1:500, ])


# For each POEM score value from 0 to 4, the distribution of probabilities of replicated POEM score == value for each replication. 
# With vertical line that represents the probability of true POEM score == value in the dataset.
lapply(0:4, function(q) {HuraultMisc::post_pred_pval(yrep = observed_POEM_rep, 
                                                     y = df$observed_POEM, 
                                                     test_statistic = function(x) {mean(x == q)}, 
                                                     plot = T)
                        }$plot + labs(x = paste0("mean(score == ", q, ")"))) %>% 
  cowplot::plot_grid(plotlist = ., ncol = 2) 
# A5 landscape

# Can also get the corresponding p-values. The greater the p-value the better, because it means that the test statistic
# in the true dataset isn't significantly different to the test statistic of the replications.
lapply(0:4, function(q) {HuraultMisc::post_pred_pval(yrep = observed_POEM_rep, 
                                                     y = df$observed_POEM, 
                                                     test_statistic = function(x) {mean(x == q)}, 
                                                     plot = F)})


# RMSE of observed_POEM_rep vs. calculated_POEM for each of the replications.
# To see if the recalled days replications still capture the big difference that we saw in observed vs. calculated POEM) in the POEM investigation.
ggplot(data = data.frame(x = apply(observed_POEM_rep, 1, function(x) {rmse(df$calculated_POEM, x)}))) +
  geom_density(aes_string(x = "x"), fill = "#9ecae1", alpha = 0.8) +
  geom_vline(xintercept = rmse(df$observed_POEM, df$calculated_POEM), size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "RMSE") +
  theme_bw(base_size = 15)
# A6 landscape

# RMSE of observed_POEM_rep vs. observed_POEM for each of the replications.
ggplot(data = data.frame(x = apply(observed_POEM_rep, 1, function(x) {rmse(df$observed_POEM, x)}))) +
  geom_density(aes_string(x = "x"), fill = "#9ecae1", alpha = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "RMSE") +
  theme_bw(base_size = 15)


# - Posterior shrinkage (plot prior vs. posterior)
par <- HuraultMisc::summary_statistics(fit, pars = c("sigma","b")) %>% select(subset = -c("Index"))
par0 <- readRDS(paste0(fit_dir, "/par0_", mdl_name, ".rds"))

HuraultMisc::plot_prior_posterior(par0, par, pars = c("sigma","b")) 
# A6 landscape




# COMBINE ITEMS --------------------------------------------------------------------------------------

# - Plot posterior predictive distribution for 3 example patients.

rep <- data.frame(Patients = c(1, 2, 3))
for (item in q_labels){
  
  fit <- readRDS(file = paste0(fit_dir, "/fit_", mdl_name, "_", item, ".rds"))
  
  recalled_days_fake <- rstan::extract(fit, pars = "recalled_days_fake")[[1]]
  # Convert recalled days to POEM score.
  POEM_dict <- data.frame(days = seq(0, 7), scores=c(0, 1, 1, 2, 2, 3, 3, 4))
  observed_POEM_rep <- plyr::mapvalues(recalled_days_fake, from = POEM_dict$days, to = POEM_dict$scores)
  
  for (pt in c(1, 2, 3)){
    rep[pt, item][[1]] <- list(observed_POEM_rep[ , pt])
  }
}

# Aggregate symptoms to produce total POEM score from 0-28.
rep <- rep %>%
  mutate(POEM_rep = pmap(list(itchy, sleep, bleeding, oozing, cracked, flaky, dry),
                          function(...) {Reduce(`+`, list(...))}))


# Now plot POEM posterior predictive distribution for all 3 example patients.
mild_pt_dPOEM <- sum(mild_pt %>% select(calculated_POEM))
moderate_pt_dPOEM <- sum(moderate_pt %>% select(calculated_POEM))
severe_pt_dPOEM <- sum(severe_pt %>% select(calculated_POEM))
d_POEM <- c(mild_pt_dPOEM, moderate_pt_dPOEM, severe_pt_dPOEM)
  
mild_pt_symp <- mild_pt %>% select(Item, calculated_POEM) %>% rename(pt1 = calculated_POEM)
moderate_pt_symp <- moderate_pt %>% select(calculated_POEM) %>% rename(pt2 = calculated_POEM)
severe_pt_symp <- severe_pt %>% select(calculated_POEM) %>% rename(pt3 = calculated_POEM)
d_POEM_symp <- cbind(mild_pt_symp, moderate_pt_symp, severe_pt_symp)
row.names(d_POEM_symp) <- q_labels

for(pt in c(1, 2, 3)){
  # For total POEM score
  ssd <- HuraultMisc::extract_distribution(object = unlist(rep[pt, ]$POEM_rep), type = "discrete", support = seq(0,28,1)) # probability mass function is type "discrete"

  pl <- ggplot(data = ssd, aes(x = Value, y = Probability, 
                                      fill = factor(ifelse(Value == d_POEM[pt], "Highlighted", "Normal")))) +
    geom_col(color = "black", show.legend = FALSE) + #fill = "grey"
    scale_fill_manual(name = "Value", values = c("#8337a3", "grey")) +
    scale_x_continuous(labels = seq(0,28,4), breaks = seq(0,28,4), limits = c(-1, 29), expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0), labels = seq(0, 0.20, 0.02), breaks = seq(0, 0.20, 0.02), limits = c(0, 0.20)) + #labels = seq(0, max(ssd$Probability), 0.02)
    labs(x = "Predicted r-POEM", y = "Probability") +
    geom_vline(aes(xintercept = mean(unlist(rep[pt, ]$POEM_rep))), # plot mean r-POEM predicted by model
               color = "blue", linetype = "dashed", size = 1)  +
    theme_bw(base_size = 15) +
    theme(panel.grid.minor.x = element_blank())
  
  
  # For POEM score of each symptom.
  for (item in q_labels){
    ssd <- HuraultMisc::extract_distribution(object = unlist(rep[pt, ][[item]]), type = "discrete", support = seq(0,4,1)) # probability mass function is type "discrete"

    pt_col_names <- c("pt1", "pt2", "pt3")
    pl_symp <- ggplot(data = ssd, aes(x = Value, y = Probability, 
                                      fill = factor(ifelse(Value == d_POEM_symp[item, pt_col_names[pt]], "Highlighted", "Normal")))) +
      geom_col(color = "black", show.legend = FALSE) + #fill = "grey"
      scale_fill_manual(name = "Value", values = c("#8337a3", "grey")) +
      scale_x_continuous(labels = seq(0, 4, 1), breaks = seq(0, 4, 1), limits = c(-1, 5), expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), labels = seq(0, max(ssd$Probability), 0.1), breaks = seq(0, max(ssd$Probability), 0.1), limits = c(-0.01, max(ssd$Probability) + 0.01)) +
      labs(x = "Predicted r-POEM", y = "Probability") +
      geom_vline(aes(xintercept = mean(unlist(rep[pt, ][[item]]))), # plot mean r-POEM predicted by model
                 color = "blue", linetype = "dashed", size = 1)  +
      theme_bw(base_size = 15) +
      theme(panel.grid.minor.x = element_blank())
  }
}

# Average and 95% CI of predicted r-POEM for each of the example patients
ci_df <- data.frame()
for(pt in c(1, 2, 3)){
  ci <- bayesplot::mcmc_intervals_data(as.data.frame(unlist(rep[pt, ]$POEM_rep)), prob = 0.95, prob_outer = 0.95, point_est = "mean") %>% 
    rename(patient = parameter) %>% 
    mutate(patient = pt)
  ci_df <- rbind(ci_df, ci)
}



# - Plot recall noise and bias for all 7 items.

# For all 7 items, retrieve recall noise (sigma) and recall bias (b) from fit.
l_params <- lapply(q_labels,
                  function(q) {
                    fit <- readRDS(file = paste0(fit_dir, "/fit_", mdl_name, "_", q, ".rds"))
                    sigma <- rstan::extract(fit, pars = c("sigma"))[[1]]
                    b <- rstan::extract(fit, pars = c("b"))[[1]]
                    return(list(sigma = sigma, b = b))
                  })
l_params <- lapply(1:2, function(i) lapply(l_params, "[[", i))
sigma <- l_params[[1]]
b <- l_params[[2]]
names(sigma) <- q_labels
names(b) <- q_labels


# Plot 80 and 95% confidence intervals of recall noise and bias for each item
p1 <- bayesplot::mcmc_intervals(as.data.frame(sigma), prob = 0.8, prob_outer = 0.95) + labs(x = "Recall noise (pseudo-days)")+ xlim(c(0, 4)) + 
  theme_bw(base_size = 15)

p2 <- bayesplot::mcmc_intervals(as.data.frame(b), prob = 0.8, prob_outer = 0.95) + labs(x = "Recall bias (pseudo-days)") + # b
  theme_bw(base_size = 15) + scale_x_continuous(limits = c(-4.5, 4.5), breaks = seq(-4, 4, by = 1)) +
  geom_vline(xintercept = 0, colour = "black") + 
  annotate("text",  x=-Inf, y = -Inf, label = "recall fewer days", vjust=-0.3, hjust=-0.1, colour = "black") + 
  annotate("text",  x=Inf, y = -Inf, label = "recall more days",vjust = -0.3, hjust=1.1, colour = "black") 

