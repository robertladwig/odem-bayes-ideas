data {
  // Parameters of priors
  real<lower=0> NEP_mu_min;
  real<lower=0> NEP_mu_max;
  real<lower=0> NEP_sigma;
  real<lower=0> SED_mu_min;
  real<lower=0> SED_mu_max;
  real<lower=0> SED_sigma;

  // Error distributions
  real<lower=0> err_sigma;

  // Data dimensions
  int<lower=1> d; // number of dates
  int N_obs; // number of dates with observations
  int<lower=1, upper=d> ii_obs[N_obs];

  // Data
  real DO_epi_init; // need a starting point. this is one option.
  real khalf;
  real theta[d];
  real volume[d];
  real area[d];
  real k600[d];
  real o2sat[d];
  real tddepth[d];
  real DO_obs_epi[N_obs]; // how do we accept missing data? https://mc-stan.org/docs/2_23/stan-users-guide/sliced-missing-data.html
}
parameters {
  real<lower=0> NEP_mu;
  real<lower=0> NEP[d];
  real<lower=0> SED_mu;
  real<lower=0> SED[d];
}
transformed parameters {
  real DO_epi[d];
  real dDOdt_epi[d];

  DO_epi[1] = DO_epi_init;
  dDOdt_epi[1] = 0;
  for(i in 2:d) {
    dDOdt_epi[i] = NEP[i-1] * (DO_epi[i-1]/(khalf + DO_epi[i-1])) * theta[i] -
    SED[i-1] *  (DO_epi[i-1]/(khalf + DO_epi[i-1])) * theta[i] * area[i-1]/volume[i-1] +
    k600[i-1] * (o2sat[d-1] - DO_epi[i-1])/tddepth[i-1];
    DO_epi[i] =  (DO_epi[i-1] + dDOdt_epi[i]) * volume[i-1]/volume[i];
    //if(stratified) {
    //  DO_epi[i] = DO_epi[i-1] + dDOdt_epi[i];
    //  DO_hypo[i] = DO_hypo[i-1] + dDOdt_hypo[i];
    //} else {
    //  DO_epi[i] = DO_epi[i];
    //  DO_hypo[i] = DO_epi[i];
    //}
  }
}
model {
  NEP_mu ~ uniform(NEP_mu_min, NEP_mu_max);
  NEP ~ normal(NEP_mu, NEP_sigma);
  SED_mu ~ uniform(SED_mu_min, SED_mu_max);
  SED ~ normal(SED_mu, SED_sigma);
  for(i in 1:N_obs) {
   DO_obs_epi[i] ~ normal(DO_epi[ii_obs[N_obs]], err_sigma); //normal(DO_obs_epi[N_obs], err_sigma);//
   // DO_obs_epi[i] ~ normal(DO_obs_epi[N_obs], err_sigma);
  }
  //for(i in 1:N_obs_epi) {
  //  DO_obs_epi[i] ~ normal(DO_epi[ii_obs[N_obs]], err_sigma);
  //}
  //for(i in 1:N_obs_hypo) {
  //  DO_obs_hypo[i] ~ normal(DO_hypo[ii_obs[N_obs]], err_sigma);
  //}
}
