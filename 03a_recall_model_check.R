#---------------- About ------------------
# Project: Recall error in the Patient-Oriented Eczema Measure (POEM) 

# This script performs a prior predictive check for the recalled days model.


# Initialisation ----------------------------------------------------------

library(tidyverse)
library(here)
library(rstan)
options(mc.cores = parallel::detectCores()) # Parallel computing
rstan_options(auto_write = TRUE)

set.seed(649949093) # Reproducibility

mdl_name <- "recalled_days_model"
stan_mdl <- paste0(here(), "/Model/", mdl_name, ".stan")
fit_dir <- here("results/Recall error model")

run_prior <- FALSE
n_chains <- 4
n_it <- 4000

if (run_prior) {
  compiled_model <- stan_model(stan_mdl)
}


# Prior predictive check -------------------------------------------------

# See what the model (based only on the prior, no data) would predict recalled number of days for 10 random fake patients.
N_random <- 10
calculated_days_random <- sample(seq(1,7), N_random, replace=TRUE)
# Instead of using 10 random fake patients, could also use actual dataset (N and calculated_days only), and see what prior would predict as r-POEM.

if (run_prior) {
  
  data_prior <- list(
    N = N_random,
    calculated_days = calculated_days_random, 
    observed_POEM = vector(),
    N_fake = 0,
    day_fake = vector(),
    run = 0 # Don't evaluate the likelihood.
  )
  
  fit_prior <- sampling(compiled_model,
                        data = data_prior,
                        iter = n_it,
                        chains = n_chains)
  
  saveRDS(fit_prior, file = paste0(fit_dir, "/prior_", mdl_name, ".rds"))
  par0 <- HuraultMisc::summary_statistics(fit_prior, pars = c("sigma","b")) %>% select(subset = -c("Index"))
  saveRDS(par0, file = paste0(fit_dir, "/par0_", mdl_name, ".rds"))
  
} else {
  
  fit_prior <- readRDS(paste0(fit_dir, "/prior_", mdl_name, ".rds"))
  par0 <- readRDS(paste0(fit_dir, "/par0_", mdl_name, ".rds"))
  
}

check_hmc_diagnostics(fit_prior)

pairs(fit_prior, pars = c("b", "sigma"))
plot(fit_prior, pars = c('b', 'sigma'), plotfun = "trace")

# Visualize prior distribution
plot(fit_prior,
     pars = c("b", "sigma"),
     plotfun = "hist")


# Plot for one of the randomly generated patients, what their predicted r-POEM score would be by the model.
recalled_days <- rstan::extract(fit_prior, pars = "recalled_days")[[1]]
# Convert recalled days to (r-)POEM score.
POEM_dict <- data.frame(days = seq(0, 7), scores=c(0, 1, 1, 2, 2, 3, 3, 4))
observed_POEM_rep <- plyr::mapvalues(recalled_days, from = POEM_dict$days, to = POEM_dict$scores)
# Convert calculated days to (d-)POEM score.
d_POEM <- plyr::mapvalues(calculated_days_random, from = POEM_dict$days, to = POEM_dict$scores)

pt = 1
ssd <- HuraultMisc::extract_distribution(object = unlist(observed_POEM_rep[, pt]), type = "discrete", support = seq(0, 4, 1)) # probability mass function is type "discrete"

pl <- ggplot(data = ssd, aes(x = Value, y = Probability, 
                                  fill = factor(ifelse(Value == d_POEM[[1]], "Highlighted", "Normal")))) +
  geom_col(color = "black", show.legend = FALSE) +
  scale_fill_manual(name = "Value", values = c("#8337a3", "grey")) +
  scale_x_continuous(labels = seq(0, 4, 1), breaks = seq(0, 4, 1), limits = c(-1, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = seq(0, max(ssd$Probability), 0.1), breaks = seq(0, max(ssd$Probability), 0.1), limits = c(-0.01, max(ssd$Probability) + 0.01)) +
  labs(x = "Predicted r-POEM", y = "Probability") +
  geom_vline(aes(xintercept = mean(unlist(observed_POEM_rep[, pt]))), # Plot mean r-POEM predicted by model.
             color = "blue", linetype = "dashed", size = 1)  +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor.x = element_blank())

