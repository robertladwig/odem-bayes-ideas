data {
  // Parameters of priors
  real<lower=0> lambda_mu_min;
  real<lower=0> lambda_mu_max;
  real<lower=0> lambda_sigma;

  // Error distributions
  real<lower=0> err_sigma;

  // Data dimensions
  int<lower=1> d; // number of dates
  int N_obs; // number of dates with observations
  int<lower=1, upper=d> ii_obs[N_obs];

  // Data
  real DO_tot_init; // need a starting point. this is one option.
  real DO_obs_tot[N_obs]; // how do we accept missing data? https://mc-stan.org/docs/2_23/stan-users-guide/sliced-missing-data.html
}
parameters {
  real<lower=0> lambda_mu;
  real<lower=0> lambda[d];
}
transformed parameters {
  real DO_tot[d];
  real dDOdt_tot[d];

  DO_tot[1] = DO_tot_init;
  dDOdt_tot[1] = 0;
  for(i in 2:d) {
    dDOdt_tot[i] = -lambda[i] * DO_tot[i-1];
    DO_tot[i] = DO_tot[i-1] + dDOdt_tot[i];
    //if(stratified) {
    //  DO_epi[i] = DO_epi[i-1] + dDOdt_epi[i];
    //  DO_hypo[i] = DO_hypo[i-1] + dDOdt_hypo[i];
    //} else {
    //  DO_epi[i] = DO_tot[i];
    //  DO_hypo[i] = DO_tot[i];
    //}
  }
}
model {
  lambda_mu ~ uniform(lambda_mu_min, lambda_mu_max);
  lambda ~ normal(lambda_mu, lambda_sigma);
  for(i in 1:N_obs) {
    DO_obs_tot[i] ~ normal(DO_tot[ii_obs[N_obs]], err_sigma);
  }
  //for(i in 1:N_obs_epi) {
  //  DO_obs_epi[i] ~ normal(DO_epi[ii_obs[N_obs]], err_sigma);
  //}
  //for(i in 1:N_obs_hypo) {
  //  DO_obs_hypo[i] ~ normal(DO_hypo[ii_obs[N_obs]], err_sigma);
  //}
}
