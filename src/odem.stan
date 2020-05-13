data {
  // Parameters of priors
  real<lower=0> NEP_mu_min;
  real<lower=0> NEP_mu_max;
  real<lower=0> NEP_sigma;
  real<lower=0> SED1_mu_min;
  real<lower=0> SED1_mu_max;
  real<lower=0> SED1_sigma;
  real<lower=0> MIN_mu_min;
  real<lower=0> MIN_mu_max;
  real<lower=0> MIN_sigma;
  real<lower=0> SED2_mu_min;
  real<lower=0> SED2_mu_max;
  real<lower=0> SED2_sigma;

  // Error distributions
  real<lower=0> err_sigma;

  // Data dimensions
  int<lower=1> d; // number of dates
  int N_obs; // number of dates with observations
  int N_obs_mix;
  int<lower=1, upper=d> ii_obs[N_obs];
  int<lower=1, upper=d> ii_obs_mix[N_obs_mix];

  // Data
  real DO_epi_init; // need a starting point. this is one option.
  real DO_hyp_init;
  real DO_tot_init;
  real khalf;
  real theta0[d];
  real theta1[d];
  real theta2[d];
  real volume_tot[d];
  real volume_epi[d];
  real area_epi[d];
  real volume_hyp[d];
  real area_hyp[d];
  real k600[d];
  real o2sat[d];
  real k600t[d];
  real o2satt[d];
  real tddepth[d];
  real DO_obs_epi[N_obs]; // how do we accept missing data? https://mc-stan.org/docs/2_23/stan-users-guide/sliced-missing-data.html
  real DO_obs_hyp[N_obs];
  real DO_obs_tot[N_obs_mix];
  real stratified[d];
}
parameters {
  real<lower=0> NEP_mu;
  real<lower=0, upper = 100> NEP[d];
  real<lower=0> SED1_mu;
  real<lower=0, upper=3200> SED1[d];
  real<lower=0> MIN_mu;
  real<lower=0, upper=1000> MIN[d];
  real<lower=0> SED2_mu;
  real<lower=0, upper=3200> SED2[d];
}
transformed parameters {
  real DO_epi[d];
  real dDOdt_epi[d];
  real DO_tot[d];
  real dDOdt_tot[d];
  real DO_hyp[d];
  real dDOdt_hyp[d];
  real flux_epi[d];
  real flux_hyp[d];
  real flux_tot[d];
  real delvol_epi[d];
  real delvol_hyp[d];
  real x_do1[d];
  real x_do2[d];

  DO_epi[1] = DO_epi_init;
  DO_hyp[1] = DO_hyp_init;
  DO_tot[1] = DO_tot_init;
  dDOdt_epi[1] = 0;
  dDOdt_hyp[1] = 0;
  flux_epi[1] = 0;
  dDOdt_tot[1] = 0;
  flux_tot[1] = 0;
  flux_hyp[1] = 0;
  x_do1[1] = 0; // new
  x_do2[1] = 0; // new
  delvol_epi[1] = 0; // new
  delvol_hyp[1] = 0; // new

  for(i in 2:d) {

  if (stratified[i] == 0) {

    dDOdt_tot[i] = NEP[i-1] * (DO_tot[i-1]/(khalf + DO_tot[i-1])) * theta0[i] -
    SED1[i-1] *  (DO_tot[i-1]/(khalf + DO_tot[i-1])) * theta0[i] * area_epi[i-1]/volume_tot[i-1] +
    k600t[i-1] * (o2satt[d-1] - DO_tot[i-1])/25. -
    MIN[i-1] * (DO_tot[i-1]/(khalf + DO_tot[i-1])) * theta0[i];

    if(fabs(dDOdt_tot[i])>fabs(DO_tot[i-1])){
      if(dDOdt_tot[i] < 0){
        flux_tot[i] = - DO_tot[i-1];
      }else{
         flux_tot[i] = dDOdt_tot[i];
      }
    }else{
      flux_tot[i] = DO_tot[i];
    }

    DO_tot[i] =  (DO_tot[i-1] + flux_tot[i]) * volume_tot[i-1]/volume_tot[i];

    dDOdt_epi[i] = dDOdt_tot[i];
    dDOdt_hyp[i] = dDOdt_tot[i];
    DO_hyp[i] = DO_tot[i];
    DO_epi[i] = DO_tot[i];
    flux_hyp[i] = flux_tot[i];
    flux_epi[i] = flux_tot[i];
    x_do1[i] = 0; // new
    x_do2[i] = 0; // new
    delvol_epi[i] = 0; // new
    delvol_hyp[i] = 0;

  } else if (stratified[i] == 1){

    delvol_epi[i] = (volume_epi[i] -  volume_epi[i-1])/volume_epi[i-1];
    if (delvol_epi[i] >= 0){
      x_do1[i] = DO_hyp[i-1];
    } else {
      x_do1[i] = DO_epi[i-1];
    }

    delvol_hyp[i] = (volume_hyp[i] -  volume_hyp[i-1])/volume_hyp[i-1];
    if (delvol_hyp[i] >= 0){
      x_do2[i] = DO_epi[i-1];
    } else {
      x_do2[i] = DO_hyp[i-1];
    }

    // FluxAtm multiplied by 0.1 to help fits
    dDOdt_epi[i] = NEP[i-1] * (DO_epi[i-1]/(khalf + DO_epi[i-1])) * theta1[i] -
    SED1[i-1] *  (DO_epi[i-1]/(khalf + DO_epi[i-1])) * theta1[i] * area_epi[i-1]/volume_epi[i-1] +
    k600[i-1] * (o2sat[d-1] - DO_epi[i-1])/tddepth[i-1] +
    delvol_epi[i] * x_do1[i]; // New was i-1

    if(fabs(dDOdt_epi[i])>fabs(DO_epi[i-1])){
      if(dDOdt_epi[i] < 0){
        flux_epi[i] = - DO_epi[i-1];
      }else{
         flux_epi[i] = dDOdt_epi[i];
      }
    }else{
      flux_epi[i] = dDOdt_epi[i];
    }

    DO_epi[i] =  (DO_epi[i-1] + flux_epi[i]) * volume_epi[i-1]/volume_epi[i];

    dDOdt_hyp[i] =  - MIN[i-1] * (DO_hyp[i-1]/(khalf + DO_hyp[i-1])) * theta2[i] -
    SED2[i-1] *  (DO_hyp[i-1]/(khalf + DO_hyp[i-1])) * theta2[i] * area_hyp[i-1]/volume_hyp[i-1] +
    delvol_hyp[i] * x_do2[i]; // New, was i-1

    if(fabs(dDOdt_hyp[i])>fabs(DO_hyp[i-1])){
      if(dDOdt_hyp[i] < 0){
        flux_hyp[i] = - DO_hyp[i-1];
      }else{
         flux_hyp[i] = dDOdt_hyp[i];
      }
    }else{
      flux_hyp[i] = dDOdt_hyp[i];
    }

    DO_hyp[i] =  (DO_hyp[i-1] + flux_hyp[i])* volume_hyp[i-1]/volume_hyp[i];

    dDOdt_tot[i] = dDOdt_epi[i];
    DO_tot[i] = (DO_epi[i-1] + DO_hyp[i-1] ) * volume_tot[i-1]/volume_tot[i];
    flux_tot[i] = flux_epi[i];

    }


  }
}
model {
  // Variables declared in the model block cannot be used elsewhere
  // NEP_mu ~  uniform(0,1);// uniform(NEP_mu_min, NEP_mu_max);
  // NEP ~ normal(0.5, 0.01);//normal(NEP_mu, NEP_sigma);
  // SED1_mu ~ uniform(SED1_mu_min, SED1_mu_max);
  // SED1 ~ normal(100, 1);
  // MIN_mu ~ uniform(MIN_mu_min, MIN_mu_max);
  // MIN ~ normal(2000, 1);
  // SED2_mu ~ uniform(SED2_mu_min, SED2_mu_max);
  // SED2 ~ normal(1000, 1);

  for(i in 2:d) {
    NEP[i] ~ normal(NEP[i-1],1);
    MIN[i] ~ normal(MIN[i-1],100);
    SED1[i] ~ normal(SED1[i-1],100);
    SED2[i] ~ normal(SED2[i-1],100);
  }

  for(i in 1:N_obs) {
   // NEP[ii_obs[i]] ~ normal(100,1);
   // MIN[ii_obs[i]] ~ normal(3000,100);
   // SED1[ii_obs[i]] ~ normal(3000,100);
   // SED2[ii_obs[i]] ~ normal(3000,100);
   DO_obs_epi[i] ~ normal(DO_epi[ii_obs[i]], 10); // error in "i", sigma of 1 way too low
   DO_obs_hyp[i] ~ normal(D_Ohyp[ii_obs[i]], 10); // error in "i", sigma of 1 way too low
  }

  for(i in 1:N_obs_mix) {
   // NEP[ii_obs_mix[i]] ~ normal(100,1);
   // MIN[ii_obs_mix[i]] ~ normal(3000,100);
   // SED1[ii_obs_mix[i]] ~ normal(3000,100);
   DO_obs_tot[i] ~ normal(DO_tot[ii_obs_mix[i]], 10); // error in "i", sigma of 1 way too low
  }

}
generated quantities {
  real Fnep[d];
  real Fsed1[d];
  real Fatm[d];
  real Fmin[d];
  real Fsed2[d];
  real Ftotepi[d];
  real Ftotepi2[d];
  real Fentrain_epi[d];
  real Fentrain_hyp[d];
  Ftotepi2[1] = DO_epi_init;
  Ftotepi[1] = 0;
  Fnep[1] = 0;
  Fsed1[1] = 0;
  Fatm[1] = 0;
  Fmin[1] = 0;
  Fsed2[1] = 0;
  Fentrain_epi[1] = 0;
  Fentrain_hyp[1] = 0;
  for (i in 2:d){
    Ftotepi2[i] = (DO_epi[i-1] + flux_epi[i]) * volume_epi[i-1]/volume_epi[i];
    Ftotepi[i] = NEP[i-1] * (DO_epi[i-1]/(khalf + DO_epi[i-1])) * theta1[i] -
    SED1[i-1] *  (DO_epi[i-1]/(khalf + DO_epi[i-1])) * theta1[i] * area_epi[i-1]/volume_epi[i-1] +
    k600[i-1] * (o2sat[d-1] - DO_epi[i-1])/tddepth[i-1];
    Fnep[i] = NEP[i-1] * (DO_epi[i-1]/(khalf + DO_epi[i-1])) * theta1[i];
    Fsed1[i] = -SED1[i-1] *  (DO_epi[i-1]/(khalf + DO_epi[i-1])) * theta1[i] * area_epi[i-1]/volume_epi[i-1];
    Fatm[i] = k600[i-1] * (o2sat[d-1] - DO_epi[i-1])/tddepth[i-1];
    Fmin[i] =  - MIN[i-1] * (DO_hyp[i-1]/(khalf + DO_hyp[i-1])) * theta2[i];
    Fsed2[i] = -SED2[i-1] *  (DO_hyp[i-1]/(khalf + DO_hyp[i-1])) * theta2[i] * area_hyp[i-1]/volume_hyp[i-1];
    Fentrain_epi[i] = delvol_epi[i] * x_do1[i];
    Fentrain_hyp[i] = delvol_hyp[i] * x_do2[i];
  }

}
