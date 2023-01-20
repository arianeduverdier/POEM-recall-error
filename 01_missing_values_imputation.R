#---------------- About ------------------
# Project: Recall error in the Patient-Oriented Eczema Measure (POEM) 

# This script imputes missing daily symptom presence/absence using a Markov chain model.


# Initialisation ----------------------------------------------------------

library(tidyverse)
library(EczemaPred)
library(HuraultMisc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstan)
library(Metrics)
options(mc.cores = parallel::detectCores()) # Parallel computing

set.seed(2021) # For reproducibility

source("functions.R")

plot_dir <- here::here("Plots/")
fits_dir <- here::here("results/Missing values")

q_labels <- c("bleeding", "cracked", "dry", "flaky", "itchy", "oozing", "sleep")

max_score <- 10 # daily scores range from 0 to 10
thr <- 2 # daily scores >= 2 are considered as present, absent otherwise

n_chains <- 4
n_it <- 2000

fit_cross_validation <- FALSE



# Load in dataset --------------------------------------------------------

daily_df <- load_daily(remove_NA = TRUE)
POEM_df <- load_POEM()

# Obtain the observed POEM scores where we do not have all 7 days present.
POEM_df <- generate_numdays_present(daily_df = daily_df, POEM_df = POEM_df)
POEM_missing_daily <- POEM_df %>% filter(num_days_daily_scores_present < 7)

# Obtain the indices (days) of the missing days in the daily_df needed to calculate POEM scores
# corresponding to those measured in POEM_df.
daily_missing_df <- missing_daily_days(daily_df = daily_df, POEM_df = POEM_df)

# Obtain the last previous daily symptom recording for each missing value.
find_last_occurence.factorial <- function(pt, day){
  stopifnot(day > 0)
  df <- daily_df %>% filter(Patient == pt, Analysis.Relative.Day == day)
  if(nrow(df) > 0){
    return(df %>% select(all_of(c(q_labels, "Analysis.Relative.Day"))))
  } else{
    find_last_occurence.factorial(pt = pt, day = (day - 1))
  }
}

daily_missing_pastinfo_df <- apply(daily_missing_df, 1, function(x){
  last_occ <- find_last_occurence.factorial(pt = x["Patient"], day = x["Analysis.Relative.Day"]) %>% 
    rename(Last.Occurence.Day = Analysis.Relative.Day) %>% 
    mutate(Patient = x["Patient"], 
           Analysis.Relative.Day = x["Analysis.Relative.Day"], 
           dt = Analysis.Relative.Day - Last.Occurence.Day)
}) %>% bind_rows()




# 1. Predicting missing daily scores- PPC and CV -----------------------

# Use the Markov chain model implemented in EczemaPred: https://ghurault.github.io/EczemaPred/articles/MC.html
K <- max_score + 1
# Requirements of EczemaPred MC model:
#- train and test should have columns y0 (for the current state), y1 (for the next state) and dt (for the time delay between states).
#- y0 and y1 should take integer values between 1 and K.
#- dt should take integer values greater than or equal to 1.
#- Missing values are not allowed.


# ---- A. Prior predictive check

t_max <- max(daily_df$Analysis.Relative.Day, POEM_df$Analysis.Relative.Day)
# y0 and y1 must take integer values between 1 and K. but our scores range from 0 to 10, so will need to shift them by 1
df0 <- data.frame(y0 = rep(1:K, each = t_max),
                  y1 = 1, # does not matter
                  dt = 1)

#  the categories are ordinal, we may assume that the probabilities of transitioning from state 5 to 1 is much smaller than the probability from transitioning from state 5 to 4. In that case, we could use an RBF-like Dirichlet prior for the transition probabilities
prior_RBF_MC <- function(K, alpha, l) {
  # RBF-like Dirichlet prior
  # Each transition pmf is normalised to a pseudo-count of K
  #
  # Args:
  # K : number of states
  # alpha: scaling factor, determining the concentration of the Dirichlet distribution
  # l: length scale
  #
  # Returns:
  # K*K matrix
  
  p <- matrix(NA, nrow = K, ncol = K)
  for (i in 1:K) {
    for (j in 1:K) {
      p[i, j] <- exp(-(i - j)^2 / 2 / l^2)
    }
    p[i, ] <- p[i, ] / sum(p[i, ])
  }
  p <- alpha * K * p
  return(list(p = p))
}

model <- EczemaModel("MC", K = K, prior = prior_RBF_MC(K, 1, 5))

fit_prior <- sample_prior(model, data = df0, chains = 1, refresh = 0)

pl <- plot_transition_MC(fit_prior) + scale_x_continuous(breaks = seq(1:11)) + scale_y_continuous(breaks = seq(1:11))


# ---- B. Cross-validation of predictions 

# Create 10 folds separating by patients (train on 9 folds and test on 1 fold) and obtain the predicted presence/absence.
# Estimate general performance by averaging results across all the test folds.
model <- EczemaModel("MC", K = K, prior = prior_RBF_MC(K, 1, 5))

# Assign patients randomly to 1 of 10 folds
k <- 10 # number of folds
seq_folds <- seq(1:k)
pts <- sample(unique(daily_df$Patient)) # Randomly shuffle patients
patient_folds_assignment <- split(pts, pts%%k)

# Create a dataframe to store all the results of the cross validation, particularly the metrics when evaluated on the test dataset
cv_results_df <- data.frame(Item = rep(q_labels, each = k), Fold = rep(seq_folds, times = length(q_labels)), AUC_prob = NA, Accuracy_dich = NA)

test_full_df <- data.frame()

for (item in q_labels){ 
  # Prepare full dataframe in format needed for MC model
  df_MC <- daily_df %>%
    rename(Score = all_of(item), Time = Analysis.Relative.Day) %>% 
    select(Patient, Time, Score) %>%
    drop_na() %>%
    group_by(Patient) %>% 
    mutate(y0 = Score + 1,
           y1 = lead(Score) + 1,
           dt = lead(Time) - Time) %>%
    ungroup()
  
  for (fold_num in 1:k){
    # Train on patients in folds 1 to k minus fold_num
    patients_train <- as.numeric(unlist(patient_folds_assignment[seq_folds[-fold_num]]))
    train_df <- df_MC %>%
      filter(Patient %in% patients_train) %>%
      select(y0, y1, dt) %>%
      drop_na()
    
    # Test on patients in fold_num
    patients_test <- as.numeric(unlist(patient_folds_assignment[fold_num]))
    test_df <- df_MC %>%
      filter(Patient %in% patients_test) %>%
      select(Patient, Time, y0, y1, dt) %>%
      drop_na()
    
    if(fit_cross_validation){
      # Fit the model
      fit <- EczemaFit(model,
                       train = train_df,
                       iter = n_it,
                       chains = n_chains,
                       pars = c("p", "P"),
                       init = 0)
      
      saveRDS(fit, file = paste0(fits_dir,"/fit_MC_cv_", fold_num, "_", item, ".rds"))
    } else{
      # Load in the fit
      fit <- readRDS(file = paste0(fits_dir,"/fit_MC_cv_", fold_num, "_", item, ".rds"))
    }
    
    test_df <- predict_daily_presabs(fit, test_df)
    
    test_df <- test_df %>%
      mutate(symp = as.numeric(y1 > thr), # we added + 1, so it is >thr rather than >= thr
             err = symp - prob_symp) %>%
      mutate(Fold = fold_num, Item = item)
    
    test_full_df <- rbind(test_full_df, test_df)
    
    # Obtain the AUC and the confusion matrix
    # true pres/abs of symptom vs. predicted prob of pres/abs of symptom (0 to 1)
    auc_prob <- as.numeric(pROC::auc(response = test_df$symp, predictor = test_df$prob_symp))
    
    # true pres/abs of symptom vs. predicted dichotomized prob of pres/abs of symptom (0 or 1)
    confusion_matrix <- caret::confusionMatrix(data = factor(test_df$impute, levels = 0:1), # df_MC$impute is the dichotomized version of the daily score predictions to be either 0 or 1
                                               reference = factor(test_df$symp, levels = 0:1))
    
    accuracy_dichotomized <- confusion_matrix$overall[["Accuracy"]]
    
    # Add metrics to the cv results dataframe.
    cv_results_df[(cv_results_df$Item == item) & (cv_results_df$Fold == fold_num), ] <- list(Item = item, Fold = fold_num, AUC_prob = auc_prob, Accuracy_dich = accuracy_dichotomized)
  }
}


# Plot cross validation results of AUROC for each item for horizon = 1
pl1 <- test_full_df %>%
  filter(dt == 1) %>%
  group_by(Item) %>%
  summarize(auroc = as.numeric(pROC::ci.auc(response = symp, predictor = prob_symp))[2],
            auroc_low_95 = as.numeric(pROC::ci.auc(response = symp, predictor = prob_symp))[1],
            auroc_high_95 = as.numeric(pROC::ci.auc(response = symp, predictor = prob_symp))[3]) %>%
  ggplot(aes(x = auroc, y = Item, xmin = auroc_low_95, xmax = auroc_high_95)) +
  geom_pointrange(position = position_dodge(width = .66), color = "#03396c") + #"#03396c" OR "#005b96"
  geom_point(size = 4, color = "#03396c") +
  geom_point(size = 3.3, color = "#d1e1ec") +
  theme_bw(base_size = 15) + 
  labs(y = "", x = "AUROC") +
  scale_x_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
  scale_y_discrete(limits=rev)



# Plot cross validation results of accuracy as a function of prediction horizon
pred <- tibble(dt = 1:30, Fold = 0) 

b <- lapply(q_labels,
            function(q) {
              
              perf_fit <- test_full_df %>%
                filter(Item == q) %>%
                mutate(Acc = 1 - abs(symp - impute)) %>%
                gamm4::gamm4(Acc ~ s(dt, bs = "cr"), # a Generative Additive Mixed effect model; s includes the terms that are smoothed, bs="cr" means penalized with cubic regression splines
                             random =~ (1 | Fold),
                             data = .)
              
              pred_fit <- predict(perf_fit$gam, newdata = pred, se.fit = TRUE)
              
              pred %>%
                mutate(Mean = pred_fit$fit,
                       SE = pred_fit$se.fit,
                       Item = q)
            }) %>%
  bind_rows()

prop_test <- test_full_df %>%
  group_by(dt) %>%
  summarise(N = n() / 7)

pl2 <- b %>%
  ggplot(aes(x = dt, y = Mean, colour = Item)) +
  geom_line(linewidth = 1.5) +
  coord_cartesian(ylim = c(0.5, 1), xlim = c(1, 14), expand = FALSE) +
  scale_x_continuous(breaks = 1:20,
                     sec.axis = dup_axis(breaks = prop_test$dt,
                                         labels = signif(prop_test$N),
                                         name = "Number of test observations")) +
  labs(x = "Prediction horizon (days)", y = "Accuracy", colour = "") +
  scale_colour_manual(values = HuraultMisc::cbbPalette) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor.x = element_blank())




# 2. Fit MC model to full daily dataset for each item ---------------------

# Fit model for each symptom separately.
for (item in q_labels){
  df <- daily_df %>%
    rename(Score = all_of(item), Time = Analysis.Relative.Day) %>% 
    select(Patient, Time, Score) %>%
    drop_na()
  
  df_MC <- df %>%
    group_by(Patient) %>%
    mutate(y0 = Score + 1,
           y1 = lead(Score) + 1,
           dt = lead(Time) - Time) %>%
    ungroup() %>%
    select(y0, y1, dt) %>%
    drop_na()
  
  
  fit <- EczemaFit(model,
                   train = df_MC,
                   iter = n_it,
                   chains = n_chains,
                   pars = c("p", "P"),
                   init = 0)
  
  saveRDS(fit, file = paste0(fits_dir,"/fit_MC_missing_values", item, ".rds"))
}

# Check fitted properly
check_hmc_diagnostics(fit)

neff = summary(fit)$summary[, "n_eff"]
all((neff/4e3) > 0.01)

Rhat = summary(fit)$summary[, "Rhat"]
all(Rhat < 1.1)

plot(fit, pars = "p")
plot_transition_MC(fit)

pairs(fit, pars = paste0("p[1,", 1:(max_score + 1), "]")) # transitions probabilities from state 1

plot(fit, pars = paste0("p[1,", 1:(max_score + 1), "]"), plotfun = "trace")

par0 <- HuraultMisc::summary_statistics(fit_prior, pars = "p")
par <- HuraultMisc::summary_statistics(fit, pars = "p")
HuraultMisc::plot_prior_posterior(par0, par, pars = "p", match_exact = FALSE)
HuraultMisc::plot_prior_influence(par0, par, pars = "p", match_exact = FALSE) 
plot_transition_MC(fit)

HuraultMisc::summary_statistics(fit, "p") %>%
  full_join(true_param, by = c("Variable" = "Parameter", "Index")) %>%
  rename(True = Value) %>%
  ggplot(aes(x = Variable)) +
  geom_pointrange(aes(y = Mean, ymin = `5%`, ymax = `95%`, colour = "Posterior")) +
  geom_point(aes(y = True, colour = "Truth"), size = 2) +
  coord_flip(ylim = c(0, 1)) +
  scale_colour_manual(values = c("Posterior" = "black", "Truth" = "#E69F00")) +
  labs(x = "", y = "Estimate", colour = "") +
  theme_bw(base_size = 20) +
  theme(legend.position = "top")

HuraultMisc::plot_coverage(do.call(cbind, rstan::extract(fit, pars = true_param[["Parameter"]])),
                           true_param[["Value"]])




# 3. Obtain predictions of presence/absence for missing values ------------

# Predict the missing daily presence/absence values based on the fit of each item on the full daily symptom dataset.
for (item in q_labels){
  fit <- readRDS(file = paste0(fits_dir,"/fit_MC_missing_values", item, ".rds"))
  
  predict_df_MC <- daily_missing_pastinfo_df %>%
    rename(y0 = all_of(item), Time = Analysis.Relative.Day) %>% 
    mutate(y0 = y0 + 1) %>%
    select(Patient, Time, y0, dt) %>%
    filter(dt <= 7) # prediction horizon <= 7 days
  
  predict_df_MC <- predict_daily_presabs(fit, predict_df_MC)
  
  # Impute with the dichotomized daily presence/absence probability (0 if probability < 0.5, 1 otherwise).
  daily_missing_df <- full_join(daily_missing_df, predict_df_MC %>% rename(Analysis.Relative.Day = Time) %>% select(-c("dt", "y0", "prob_symp")), by = c("Patient","Analysis.Relative.Day")) %>%
    drop_na() # Drop the instances we cannot impute.
  names(daily_missing_df)[names(daily_missing_df) == 'impute'] <- item
}




# 4. Impute daily presence score ------------------------------------------

# Dichotomize daily_df to presence/absence of symptoms if scores are >= threshold.
daily_thresholded_df <- data.frame(daily_df) # make a copy
daily_thresholded_df[q_labels] <- lapply(daily_thresholded_df[q_labels], function(x){as.numeric(x >= thr)})

# Merge the daily dataframe with the imputed missing values.
col_keep <- c("Patient","ID", "Analysis.Relative.Day", q_labels)
daily_imputed_df <- full_join(daily_thresholded_df %>% select(all_of(col_keep)), daily_missing_df, by = col_keep)

save(daily_imputed_df, file = paste0(here("daily_imputed_df.RData")))
