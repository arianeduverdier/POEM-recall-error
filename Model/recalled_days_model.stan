functions{
  real ordered_probit_POEM_lpmf(int y, real eta, vector ct) {
    
    real out;
    
    if (y > 4) {
      reject("y must take discrete values between 0 and 4");
    }
    
    if (y == 0) {
      out = ordered_probit_lpmf(0 + 1 | eta, ct);
    } else if (y == 1) {
      out = log_sum_exp(ordered_probit_lpmf(1 + 1 | eta, ct), ordered_probit_lpmf(2 + 1 | eta, ct));
    } else if (y == 2) {
      out = log_sum_exp(ordered_probit_lpmf(3 + 1 | eta, ct), ordered_probit_lpmf(4 + 1 | eta, ct));
    } else if (y == 3) {
      out = log_sum_exp(ordered_probit_lpmf(5 + 1 | eta, ct), ordered_probit_lpmf(6 + 1 | eta, ct));
    } else {
      out = ordered_probit_lpmf(7 + 1 | eta, ct);
    }
    
    return(out);
    
  }
  
}

data {
  int<lower = 0, upper = 1> run; // Switch to evaluate the likelihood.
  int<lower = 0> N; // Number of observations.
  real<lower = 0, upper = 7> calculated_days[N]; // d-POEM symptom scores.
  int<lower = 0, upper = 4> observed_POEM[N * run]; // r-POEM symptom scores.
  int<lower = 0> N_fake; // Number of example patients.
  vector[N_fake] day_fake; // Generate example patients with [pt1, pt2, pt3] calculated_days.
}

transformed data {
  vector[7] ct;
  for (i in 1:7) ct[i] = i - 0.5; // Equally spaced cut-offs.
}

parameters {
  real<lower = 0> sigma; // Recall noise (standard deviation).
  real b; // Recall bias.
}

transformed parameters {
  real eta_recalled_days[N]; // Location of the recalled days distribution (in the latent space).
  real z_eta_recalled_days[N];
  vector[7] z_ct = ct / sigma;
  
  for (i in 1:N) {
    eta_recalled_days[i] = calculated_days[i] + b;
    z_eta_recalled_days[i] = eta_recalled_days[i] / sigma;
  }
  
}

model {
  sigma ~ normal(0, 2);
  b ~ normal(0, 2);
  
  if(run == 1){
    for (i in 1:N) {
      target += ordered_probit_POEM_lpmf(observed_POEM[i] | z_eta_recalled_days[i], z_ct); 
    }
  }
}

generated quantities {
  real recalled_days[N];
  real recalled_days_fake[N_fake];
  
  for (i in 1:N) {
    recalled_days[i] = ordered_probit_rng(z_eta_recalled_days[i], z_ct) - 1;
  }
  
  // Generate 3 example patients with calculated_days[3] equal to day_fake[3].
  for (i in 1:N_fake){
    recalled_days_fake[i] = ordered_probit_rng((b + day_fake[i]) / sigma, z_ct) - 1;
  }
  
}
