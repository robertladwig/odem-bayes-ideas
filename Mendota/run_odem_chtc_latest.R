cat('\f')
rm(list= ls())

library(tidyverse)

#oneyear <- 2008
#twoyear <- 2007:2008
fiveyear <- 1979:2019

input <- readr::read_csv(
  './input.txt',
  col_names=c('datetime', 'thermocline_depth', 'temperature_epi', 'temperature_hypo', 'temperature_total', 'volume_total', 'volume_epi', 'volume_hypo', 'area_thermocline', 'area_surface', 'upper_meta', 'lower_meta', 'year', 'day_of_year', 'wind'),
  col_types=cols(datetime=col_datetime(), year=col_integer(), day_of_year=col_integer(), .default=col_double()))
#in1yr <- filter(input, year %in% oneyear)
#in2yr <- filter(input, year %in% twoyear)
in5yr <- filter(input, year %in% fiveyear)


library(zoo)
#  time stamp in the first row, then area weighted average total oxygen in the second row, area weighted average epilimnion oxygen in the third row, and area weighted average hypolimnion oxygen in the fourth row
obs <- read.table(
  './observed.txt',
  header=FALSE,
  sep=' ',
  as.is=TRUE) %>%
  t() %>%
  as_tibble(.name_repair='minimal') %>%
  setNames(., nm=c('dateint', 'DO_tot', 'DO_epi', 'DO_hypo')) %>%
  mutate(date = zoo::as.Date(dateint, origin='1979-04-01')) %>% # just guessing at origin and therefore at dates
  select(date, everything())
#obs1yr <- filter(obs, lubridate::year(date) %in% oneyear)
#obs2yr <- filter(obs, lubridate::year(date) %in% twoyear)
obs5yr <- filter(obs, lubridate::year(date) %in% fiveyear)


in1yr= in5yr
obs1yr = obs5yr
#obs1yr = obs1yr[((round(nrow(obs5yr) * 1/3))):nrow(obs5yr),]

idx1 = which(!is.na(in1yr$thermocline_depth))[1]
idx2 = rev(which(!is.na(in1yr$thermocline_depth)))[1]
in1yr = in1yr #[idx1:idx2,]
in1yr$strat <- ifelse(is.na(in1yr$thermocline_depth),0,1)
strat.pos <- c()
for (ii in 1:length(in1yr$strat)){
  if (in1yr$strat[ii] == 1 && in1yr$strat[ii-1] == 0){
    strat.pos <- append(strat.pos, ii)
  }
}

idy = match(obs1yr$date,zoo::as.Date(in1yr$datetime))
# idx = idy[!is.na(idy)]
idx = idy[!is.na(obs1yr$DO_tot)]
idxx = idy[!is.na(obs1yr$DO_epi)]
# idz = which(!is.na(idx))
idz = match(idx, idy)
idzz =match(idxx, idy)
DO_obs_epi = rep(NA, length(in1yr$datetime))
DO_obs_epi[idxx] = obs1yr$DO_epi[idzz]
DO_obs_hyp = rep(NA, length(in1yr$datetime))
DO_obs_hyp[idxx] = obs1yr$DO_hypo[idzz]
DO_obs_tot = rep(NA, length(in1yr$datetime))
DO_obs_tot[idx] = obs1yr$DO_tot[idz]

library(lazyeval)
in1yr <- in1yr[,-c(1)] %>% mutate_each( funs_( interp( ~replace(., is.na(.),0) ) ) )

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(LakeMetabolizer)
# we started with 10 days but I've played with 200 to get more "data"

simdata <- tibble(
  DO_obs_epi = DO_obs_epi * 1000,
  DO_obs_hyp = DO_obs_hyp * 1000,
  DO_obs_tot = DO_obs_tot * 1000,
  day = seq(1, nrow(in1yr))
  # DO_obs_tot_true = 10 * exp(-0.04*day) + rnorm(length(day), 0, 0.0002),
  # have_obs = ifelse(day == 1, 0, round(runif(10))),
  # DO_obs_tot = ifelse(have_obs == 1, DO_obs_tot_true, NA)
)

WhenToEstimateParams = sort(c(idxx, idx, strat.pos)) # strat.pos

nParamEstimates = length(WhenToEstimateParams)
ParamIndex = rep(nParamEstimates,nrow(simdata))
ParamIndex[1:WhenToEstimateParams[1]-1] = 1 # use the first parameter value
for (i in 2:nParamEstimates){
  iCurrent = WhenToEstimateParams[i-1]:(WhenToEstimateParams[i]-1)
  ParamIndex[iCurrent] = i
}
# End new code
##################

dummyinput <- list(
  NEP_mu_min = 0,
  NEP_mu_max = 0.5,#10000,
  NEP_sigma = 1e-32, #1000,#0.000001,
  SED1_mu_min = 0,
  SED1_mu_max = 1,#1500,
  SED1_sigma = 1e-32, #5000,
  MIN_mu_min = 0,
  MIN_mu_max = 1,#10000,
  MIN_sigma = 1e-32, #1000,#0.000001,
  SED2_mu_min = 0,
  SED2_mu_max = 1,#1500,
  SED2_sigma = 1e-32, #5000,
  theta0 = 1.08^(in1yr$temperature_total - 20),
  theta1 = 1.08^(in1yr$temperature_epi - 20),
  theta2 = 1.08^(in1yr$temperature_hypo - 20),
  k600t = k600.2.kGAS.base(k.cole.base(in1yr$wind),temperature = in1yr$temperature_total, gas = "O2"),
  o2satt = o2.at.sat.base(temp = in1yr$temperature_total, altitude = 300) * 1000,
  k600 = k600.2.kGAS.base(k.cole.base(in1yr$wind),temperature = in1yr$temperature_epi, gas = "O2"),
  o2sat = o2.at.sat.base(temp = in1yr$temperature_epi, altitude = 300) * 1000,
  volume_epi = in1yr$volume_epi,
  volume_tot = in1yr$volume_total,
  area_epi = in1yr$area_surface,
  volume_hyp = in1yr$volume_hypo,
  area_hyp = in1yr$area_thermocline,
  tddepth = in1yr$thermocline_depth,
  ii_obs = idxx,
  ii_obs_mix = idx,
  wtr_epi = in1yr$temperature_epi,
  wtr_hyp = in1yr$temperature_hypo,
  wtr_tot = in1yr$temperature_total,
  khalf = 500, # New, was 3000
  err_sigma = 0.0003,
  d = nrow(simdata),
  DO_epi_init = 15 * 1000, #simdata$DO_obs[1],
  DO_hyp_init = 15 * 1000,
  DO_tot_init = 15 * 1000,
  stratified = in1yr$strat,
  strat_pos = strat.pos,
  len_strat_pos = length(strat.pos),
  d_strat_pos = length(strat.pos),
  ##################
  # New code by Paul
  i_Param = ParamIndex,  # n_ParamEst indeces within d
  n_ParamEst = nParamEstimates # number of times params are estimated
)

# simdata$DO_obs_epi = simdata$DO_obs_epi * 1.5
dummyinput$DO_obs_epi = simdata$DO_obs_epi[idxx]
dummyinput$DO_obs_hyp = simdata$DO_obs_hyp[idxx]
dummyinput$DO_obs_tot = simdata$DO_obs_tot[idx]
dummyinput$N_obs = length(dummyinput$ii_obs)
dummyinput$N_obs_mix = length(idx)

fit <- stan(file = 'odem.stan', data = dummyinput, chains = 4, iter = 100, control=list(adapt_delta = 0.8))

fit@stanmodel@dso <- new("cxxdso")
saveRDS(fit, file = "fit_Mendota.rds")
