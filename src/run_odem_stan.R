cat('\f')
rm(list= ls())

setwd('/home/robert/Projects/DSI/odem-bayes-ideas')

library(tidyverse)

lake.list <- c('Allequash', 'BigMuskellunge', 'Crystal', 'Fish', 'Mendota',
               'Monona', 'Sparkling', 'Trout', 'Wingra')

for (lake.id in lake.list){
setwd(paste0(lake.id))

oneyear <- 2008
twoyear <- 2007:2008
fiveyear <- 1979:2019
 fiveyear <- 2010:2013

input <- readr::read_csv(
  paste0('input_',lake.id,'.txt'),
  #'input_Allequash.txt',#'inst/extdata/input.txt', inst/extdata/
  col_names=c('datetime', 'thermocline_depth', 'temperature_epi', 'temperature_hypo', 'temperature_total', 'volume_total', 'volume_epi', 'volume_hypo', 'area_thermocline', 'area_surface', 'upper_meta', 'lower_meta', 'year', 'day_of_year', 'wind','airtemp'),
  col_types=cols(datetime=col_datetime(), year=col_integer(), day_of_year=col_integer(), .default=col_double()))
in1yr <- filter(input, year %in% oneyear)
in2yr <- filter(input, year %in% twoyear)
in5yr <- filter(input, year %in% fiveyear)

ggplot(in2yr, aes(x=datetime)) +
  geom_line(aes(y=temperature_epi), color='seagreen') +
  geom_line(aes(y=temperature_hypo), color='navy') +
  theme_bw()
ggplot(in2yr, aes(x=datetime, y=thermocline_depth)) +
  geom_line() +
  theme_bw()

library(zoo)
#  time stamp in the first row, then area weighted average total oxygen in the second row, area weighted average epilimnion oxygen in the third row, and area weighted average hypolimnion oxygen in the fourth row
obs <- read.table(
  paste0('oxy_observed_',lake.id,'.txt'),
  #'oxy_observed_Allequash.txt',# 'inst/extdata/observed.txt',inst/extdata/
  header=FALSE,
  sep=' ',
  as.is=TRUE) %>%
  t() %>%
  as_tibble(.name_repair='minimal') %>%
  setNames(., nm=c('dateint', 'DO_tot', 'DO_epi', 'DO_hypo')) %>%
  mutate(date = zoo::as.Date(dateint, origin='1979-04-01')) %>% # just guessing at origin and therefore at dates
  select(date, everything())
obs1yr <- filter(obs, lubridate::year(date) %in% oneyear)
obs2yr <- filter(obs, lubridate::year(date) %in% twoyear)
obs5yr <- filter(obs, lubridate::year(date) %in% fiveyear)
ggplot(obs, aes(x=date)) +
  geom_line(aes(y=DO_tot), color='black') +
  geom_line(aes(y=DO_epi), color='seagreen') +
  geom_line(aes(y=DO_hypo), color='navy') +
  theme_bw()

#
# in1yr= in2yr
# obs1yr = obs2yr
in1yr= in5yr
obs1yr = obs5yr

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
##################
# New code by Paul
# Create a vector of length d that says which parameter (SED1, SED2, MIN, NEP) value to use
# Assume there's one parameter value per observational data point
# Break into chunks that change every time there's a parameter
# Values in the vector indicate the index of the parameter value to use
# Only change WhenToEstimateParams
# WhenToEstimateParams = idxx # In this case, every time there's observed data
# WhenToEstimateParams = c(150,250,550,650) # About twice per stratified season
# WhenToEstimateParams = c(225,600) # About twice per stratified season
WhenToEstimateParams = sort(c(idxx, idx,strat.pos))
WhenToEstimateParams = sort(c(idxx, idx))

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
  o2satt = o2.at.sat.base(temp = in1yr$temperature_total, altitude = 450) * 1000,
  k600 = k600.2.kGAS.base(k.cole.base(in1yr$wind),temperature = in1yr$temperature_epi, gas = "O2"),
  o2sat = o2.at.sat.base(temp = in1yr$temperature_epi, altitude = 450) * 1000,
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
  n_ParamEst = nParamEstimates, # number of times params are estimated
  airtemp = in1yr$airtemp
)

# simdata$DO_obs_epi = simdata$DO_obs_epi * 1.5
dummyinput$DO_obs_epi = simdata$DO_obs_epi[idxx]
dummyinput$DO_obs_hyp = simdata$DO_obs_hyp[idxx]
dummyinput$DO_obs_tot = simdata$DO_obs_tot[idx]
dummyinput$N_obs = length(dummyinput$ii_obs)
dummyinput$N_obs_mix = length(idx)
#dummyinput$theta0[which(dummyinput$theta0 <= 1.08 ** (4-20))] = 1e-5
dummyinput$k600t[which(in1yr$airtemp <= 0 & in1yr$temperature_total <= 4)] = 1e-5
dummyinput$theta0[which(in1yr$airtemp <= 0 & in1yr$temperature_total <= 4)] = 1e-5

fit <- stan(file = '../src/odem.stan', data = dummyinput, chains = 4, iter = 2000,control=list(adapt_delta = 0.8))


list_of_draws <- extract(fit)
print(names(list_of_draws))

fit_summary <- summary(fit)
print(names(fit_summary))
print(fit_summary$summary)
# fit@stanmodel@dso <- new("cxxdso")
# saveRDS(fit, file = "fit.rds")

# fit <- readRDS(file = 'fit.rds')

# fit <- readRDS(file = paste0('fit_',lake.id,'.rds'))
# an example of extracting parameters for this particular dummy model, i'm
# geting an overestimate of lambda and consequently a much faster modeled drop
# than actual drop in "DO". As an exercise (entirely optional), you might find
# it useful to explore different values of lambda hyperparameters (mu_min,
# mu_max, and sigma), data-generating parameters, MCMC iterations, etc. to see
# if you can achieve a better fit. Then again, it might be better just to jump
# straight into more realistic ODEM equations.
NEP <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$NEP %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
SED1 <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$SED1 %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
MIN <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$MIN %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
SED2 <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$SED2 %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))

ggplot(NEP, aes(x=day, y=mean, col = 'NEP')) + geom_line() +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd, col = 'NEP'), alpha=0.2) +
  geom_line(data = SED1, aes(x=day, y=mean, col = 'SED1')) +
  geom_ribbon(data = SED1, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd, col = 'SED1'), alpha=0.2) +
  geom_line(data = SED2, aes(x=day, y=mean, col = 'SED2')) +
  geom_ribbon(data = SED2, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd, col = 'SED2'), alpha=0.2) +
  geom_line(data = MIN, aes(x=day, y=mean, col = 'MIN')) +
  geom_ribbon(data = MIN, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd, col = 'MIN'), alpha=0.2) +
  theme_bw()


# NEP$doy =  in1yr$day_of_year; NEP$year = in1yr$year
# SED1$doy =  in1yr$day_of_year; SED1$year = in1yr$year
# SED2$doy =  in1yr$day_of_year; SED2$year = in1yr$year
# MIN$doy =  in1yr$day_of_year; MIN$year = in1yr$year
#
# ggplot(NEP, aes(x=doy, y=mean, col = 'NEP')) + geom_line() +
#   geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd, col = 'NEP'), alpha=0.2) +
#   geom_line(data = SED1, aes(x=doy, y=mean, col = 'SED1')) +
#   geom_ribbon(data = SED1, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd, col = 'SED1'), alpha=0.2) +
#   geom_line(data = SED2, aes(x=doy, y=mean, col = 'SED2')) +
#   geom_ribbon(data = SED2, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd, col = 'SED2'), alpha=0.2) +
#   geom_line(data = MIN, aes(x=doy, y=mean, col = 'MIN')) +
#   geom_ribbon(data = MIN, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd, col = 'MIN'), alpha=0.2) +
#   facet_wrap(~ year) + theme_bw()

DO_epi <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$DO_epi %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
DO_hyp <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$DO_hyp %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
ggplot(DO_epi, aes(x=day, y=mean)) + geom_line(col = 'red') +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_point(data=simdata, aes(x=day, y=DO_obs_epi), col ='red') +
  geom_line(data = DO_hyp, aes(x=day, y=mean), col = 'green') +
  geom_ribbon(data = DO_hyp, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_point(data=simdata, aes(x=day, y=DO_obs_hyp),col ='green') +
  geom_point(data=simdata, aes(x=day, y=DO_obs_tot),col ='blue') +
  ylim(c(0,20000)) + theme_minimal()

DO_epi$doy =  in1yr$day_of_year; DO_epi$year = in1yr$year
DO_hyp$doy =  in1yr$day_of_year; DO_hyp$year = in1yr$year
simdata$doy =  in1yr$day_of_year; simdata$year = in1yr$year

calc_fit <- function(mod_data, obs_data){
  obs <- obs_data #cbind(obs_data[3,], obs_data[4,])
  mod <- mod_data #cbind(input.values$o2_epil[proc.obs[1,]]/1000,
  #     input.values$o2_hypo[proc.obs[1,]]/1000)
  return (sqrt(mean((obs-mod)**2,na.rm = TRUE))) # RMSE
}

calc_fit(mod_data = DO_epi$mean, obs_data = simdata$DO_obs_epi)/1000
calc_fit(mod_data = DO_hyp$mean, obs_data = simdata$DO_obs_hyp)/1000
calc_fit(mod_data = cbind(DO_epi$mean,DO_hyp$mean), obs_data = cbind(simdata$DO_obs_epi, simdata$DO_obs_hyp))/1000
rmse <- round(calc_fit(mod_data = cbind(DO_epi$mean,DO_hyp$mean), obs_data = cbind(simdata$DO_obs_epi, simdata$DO_obs_hyp))/1000
              ,2)

g1 <- ggplot(DO_epi, aes(x=doy, y=mean)) + geom_line(col = 'red') +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_point(data=simdata, aes(x=doy, y=DO_obs_epi), col ='red') +
  geom_line(data = DO_hyp, aes(x=doy, y=mean), col = 'green') +
  geom_ribbon(data = DO_hyp, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_point(data=simdata, aes(x=doy, y=DO_obs_hyp),col ='green') +
  geom_point(data=simdata, aes(x=doy, y=DO_obs_tot),col ='blue') +
  facet_wrap(~ year) +
  ggtitle(paste0(lake.id,'_',rmse))+
  ylim(c(0,20000)) + theme_minimal();g1
ggsave(file = paste0(lake.id,'_conc.png'), g1, dpi = 300, width =400, height = 200,
       units='mm')



Ftotepi <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$Ftotepi %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
Ftotepi2 <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$Ftotepi2 %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
Fnep <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$Fnep %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
Fatm <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$Fatm %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
Fsed1 <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$Fsed1 %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
Fmin <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$Fmin %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
Fsed2 <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$Fsed2 %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
ENTR1 <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$Fentrain_epi %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
ENTR2 <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$Fentrain_hyp%>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))

ENTR1$doy =  in1yr$day_of_year; ENTR1$year = in1yr$year
ENTR2$doy =  in1yr$day_of_year; ENTR2$year = in1yr$year
Fsed2$doy =  in1yr$day_of_year; Fsed2$year = in1yr$year
Fsed1$doy =  in1yr$day_of_year; Fsed1$year = in1yr$year
Fmin$doy =  in1yr$day_of_year; Fmin$year = in1yr$year
Fatm$doy =  in1yr$day_of_year; Fatm$year = in1yr$year
Fnep$doy =  in1yr$day_of_year; Fnep$year = in1yr$year

g1 <- ggplot(Fnep, aes(x=doy, y=mean, col = 'NEP')) + geom_line() +
  # geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fmin, aes(x=doy, y=mean, col = 'MIN')) +
  # geom_ribbon(data = Fmin, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fsed1, aes(x=doy, y=mean, col = 'SED1')) +
  # geom_ribbon(data = Fsed1, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fsed2, aes(x=doy, y=mean, col = 'SED2')) +
  # geom_ribbon(data = Fsed2, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fatm, aes(x=doy, y=mean, col = 'ATM')) +
  # geom_ribbon(data = Fatm, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  # geom_line(data = ENTR1, aes(x=doy, y=mean, col = 'ENTR1')) +
  # geom_ribbon(data = ENTR1, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2)+
  # geom_line(data = ENTR2, aes(x=doy, y=mean, col = 'ENTR2')) +
  # geom_ribbon(data = ENTR2, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2)+
  ylim(-1500,1500)+
  ggtitle(paste0(lake.id,'_',rmse))+
  facet_wrap(~year) + theme_bw();g1
ggsave(file = paste0(lake.id,'_fluxes.png'), g1, dpi = 300, width =400, height = 200,
       units='mm')

ggplot(Fnep, aes(x=day, y=mean, col = 'NEP')) + geom_line() +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fsed1, aes(x=day, y=mean, col = 'SED1')) +
  geom_ribbon(data = Fsed1, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fatm, aes(x=day, y=mean, col = 'ATM')) +
  geom_ribbon(data = Fatm, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = DO_epi, aes(x=day, y=mean, col = 'epi')) + geom_line() +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = ENTR1, aes(x=day, y=mean, col = 'ENTR1')) +
  geom_ribbon(data = ENTR1, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2)+ theme_bw()

ggplot(Fsed2, aes(x=day, y=mean, col = 'SED2')) + geom_line() +
  geom_ribbon(data = Fsed2, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fmin, aes(x=day, y=mean, col = 'MIN')) +
  geom_ribbon(data = Fmin, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = DO_hyp, aes(x=day, y=mean, col = 'hyp')) + geom_line() +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2)  +
  geom_line(data = ENTR2, aes(x=day, y=mean, col = 'ENTR2')) +
  geom_ribbon(data = ENTR2, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2)+ theme_bw()

Fepi <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$flux_epi %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
Fhyp <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$flux_hyp %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))

ggplot(Fepi, aes(x=day, y=mean, col = 'Fepi')) + geom_line() +
  geom_ribbon(data = Fepi, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fhyp, aes(x=day, y=mean, col = 'Fhyp')) +
  geom_ribbon(data = Fhyp, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) + theme_bw()





odem_stan <- data.frame('datetime' = input$datetime,
                        in1yr, simdata[,-c(4,5,6)],
                        'DO_epi' = DO_epi$mean, 'DO_hyp' = DO_hyp$mean,
                        'Fatm' = Fatm$mean, 'Fnep' = Fnep$mean, 'Fsed1' = Fsed1$mean, 'Fentr1' = ENTR1$mean,
                        'Fsed2' = Fsed2$mean, 'Fmin' = Fmin$mean, 'Fentr2' = ENTR2$mean)
save(odem_stan, file = paste0(lake.id,'.Rda'))#'/home/robert/Desktop/odem_stan_mendota.Rda')

odem_stan <- data.frame(                        in1yr, simdata[,-c(4,5,6)],
                        'DO_epi' = DO_epi$mean, 'DO_hyp' = DO_hyp$mean,
                        'Fatm' = Fatm$mean, 'Fnep' = Fnep$mean, 'Fsed1' = Fsed1$mean, 'Fentr1' = ENTR1$mean,
                        'Fsed2' = Fsed2$mean, 'Fmin' = Fmin$mean, 'Fentr2' = ENTR2$mean)
save(odem_stan, file = paste0('sparkling_ice','.Rda'))#'/home/robert/Desktop/odem_stan_mendota.Rda')

setwd('..')
}


library(WaveletComp) #http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf
my.data <- data.frame('date' = input$datetime, 'Fnep' = scale(Fnep$mean))
head(my.data,3)
my.w <- analyze.wavelet(my.data, 'Fnep',
                        loess.span = 0,
                        dt = 1, dj = 1/50,
                        lowerPeriod = 2, upperPeriod = 8192,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, color.key = 'interval', n.levels = 250,
         legend.params = list(lab = 'wavelet power levels'),
         periodlab = 'period (days)',
         show.date = TRUE, date.format = '%F', timelab = '')
wt.image(my.w, color.key = 'interval', n.levels = 250,
         legend.params = list(lab = 'wavelet power levels', label.digits = 2),
         periodlab = 'periods (days)', maximum.level = 8,
         # Concerning item 1 above --- plot the square root of power:
         exponent = 0.5,
         # Concerning item 2 above --- time axis:
         show.date = TRUE, date.format = '%F', timelab = '',
         spec.time.axis = list(at = c(paste(1979:2019, '-01-01', sep = '')),
                               labels = c(1979:2015)),
         timetcl = -0.5,  # draws outward ticks
         # Concerning item 3 above --- period axis:
         spec.period.axis =
           list(at = c(32, 64, 128, 365, 1024)),
         periodtck = 1, periodtcl = NULL  # draws horizontal lines
)


library(biwavelet) # https://rstudio-pubs-static.s3.amazonaws.com/152496_026c3ac97f7d40e5ba0cadf757730fce.html

t1 <- cbind(seq_len(length(Fnep$mean)),scale(Fnep$mean))
t2 <- cbind(seq_len(length(Fnep$mean)),scale(Fatm$mean))
nrands = 10
wtc.AB = wtc(t1, t2, nrands = nrands)
jpeg(file="/home/robert/Desktop/odem_wavelet_atm.jpeg",width=1200, height=873)
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = 'grey', lwd.coi = 2,
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = 'Frequency', xlab = 'Time',
     plot.cb = TRUE, main = 'Wavelet Coherence: Fnep vs Fatm')
# Adding grid lines
n = length(t1[, 1])
abline(v = seq(260, n, 260), h = 1:16, col = 'brown', lty = 1, lwd = 1)
# Defining x labels
axis(side = 3, at = c(seq(0, n, 365.35)), labels = c(seq(1979, 2019, 1)))
dev.off()


# Paul's new code to look at parameters
plot(idxx, SED1$mean, type = 'l')
points(idxx,SED1$mean[idxx])

##################
# New code by Paul
# Plot atmospheric flux
# The O2 deficit is saturation - modeled, so if modeled>sat, the value is negative
# Scaler adjusts the deficit so that it can be viewed on the same y axis
Scaler = 1/8
o2sat = o2.at.sat.base(temp = in1yr$temperature_epi, altitude = 300) * 1000
o2def = (o2sat - DO_epi$mean) * Scaler
o2defObs = (o2sat - simdata$DO_obs_epi) * Scaler
plot(Fatm$mean,type='l',ylim=c(-500,500),
     xlab='Day of sim',main="")
legend('topright',col=c('black','red','blue'),lty=c(1,1,1),c('Modeled Fatm','Modeled O2 def scaled','Observed O2 def'))
lines(o2def,col='red')
points(o2defObs,col='blue')
abline(h=0)



#### clustering test
# https://rpubs.com/AlgoritmaAcademy/som
# https://en.proft.me/2016/11/29/modeling-self-organising-maps-r/



load('Allequash/Allequash.Rda')
allequash <- odem_stan
load('BigMuskellunge/BigMuskellunge.Rda')
bigmuskellunge <- odem_stan
load('Crystal/Crystal.Rda')
crystal <- odem_stan
load('Fish/Fish.Rda')
fish <- odem_stan
load('Mendota/Mendota.Rda')
mendota <- odem_stan
load('Monona/Monona.Rda')
monona <- odem_stan
load('Sparkling/Sparkling.Rda')
sparkling <- odem_stan
load('Trout/Trout.Rda')
trout <- odem_stan
load('Wingra/Wingra.Rda')
wingra <- odem_stan

library(kohonen)
library(RColorBrewer)

df.measures <- c('thermocline_depth',
                 'temperature_total',
                 'year',
                 'day_of_year',
                 'wind',
                 'DO_epi',
                 'DO_hyp',
                 'Fatm',
                 'Fnep',
                 'Fsed1',
                 'Fentr1',
                 'Fsed2',
                 'Fmin',
                 'Fentr2')
test.df <- data.frame('id' = c(rep(lake.list, each = nrow(mendota))),
                      'Fatm' = c(allequash$Fatm, bigmuskellunge$Fatm, crystal$Fatm, fish$Fatm,
                                 mendota$Fatm, monona$Fatm, sparkling$Fatm, trout$Fatm, wingra$Fatm),
                      'Fnep' = c(allequash$Fnep, bigmuskellunge$Fnep, crystal$Fnep, fish$Fnep,
                                 mendota$Fnep, monona$Fnep, sparkling$Fnep, trout$Fnep, wingra$Fnep),
                      'Fsed1' = c(allequash$Fsed1, bigmuskellunge$Fsed1, crystal$Fsed1, fish$Fsed1,
                                  mendota$Fsed1, monona$Fsed1, sparkling$Fsed1, trout$Fsed1, wingra$Fsed1),
                      'Fentr1' = c(allequash$Fentr1, bigmuskellunge$Fentr1, crystal$Fentr1, fish$Fentr1,
                                   mendota$Fentr1, monona$Fentr1, sparkling$Fentr1, trout$Fentr1, wingra$Fentr1),
                      'Fsed2' = c(allequash$Fsed2, bigmuskellunge$Fsed2, crystal$Fsed2, fish$Fsed2,
                                  mendota$Fsed2, monona$Fsed2, sparkling$Fsed2, trout$Fsed2, wingra$Fsed2),
                      'Fmin' = c(allequash$Fmin, bigmuskellunge$Fmin, crystal$Fmin, fish$Fmin,
                                 mendota$Fmin, monona$Fmin, sparkling$Fmin, trout$Fmin, wingra$Fmin),
                      'Fentr2' = c(allequash$Fentr2, bigmuskellunge$Fentr2, crystal$Fentr2, fish$Fentr2,
                                  mendota$Fentr2, monona$Fentr2, sparkling$Fentr2, trout$Fentr2, wingra$Fentr2))

ads.train <- scale(test.df[,-c(1)])
ads.grid <- somgrid(xdim = 10, ydim = 10, topo = 'hexagonal')
test.som <- som(ads.train,
                   grid = ads.grid, dist.fcts = 'euclidean')
plot(test.som, palette.name = rainbow)

heatmap.som <- function(model){
  for (i in 1:7) {
    plot(model, type = "property", property = getCodes(model)[,i],
         main = colnames(getCodes(model))[i])
  }
}
heatmap.som(test.som)

colors <- function(n, alpha =1){
  rev(heat.colors(n, alpha))
}
plot(test.som, type = 'counts', palette.name = colors,
     heatkey = TRUE)

library(factoextra)
set.seed(100)
fviz_nbclust(test.som$codes[[1]], kmeans, method = 'wss')
clust <- kmeans(test.som$codes[[1]], 6)

int <- sample(nrow(test.df), nrow(test.df) * 0.8)
train <- test.df[int, ]
test <- test.df[-int, ]

trainX <- scale(train[,-1])
testX <- scale(test[,-1], center = attr(trainX, 'scaled:center'))

train.label <- factor(train[,1])
test.label <- factor(test[,1])
test[,1]
testXY <- list(independent = testX, dependent = test.label)

class <- xyf(trainX, classvec2classmat(train.label), ads.grid, rlen = 500)
plot(class, type = 'changes')

pred <- predict(class, newdata = testXY)
table(Predict = pred$predictions[[2]], Actual = test.label)

plot(test.som, type = 'codes', bgcol = rainbow(9)[clust$cluster], main = 'Cluster SOM')
add.cluster.boundaries(test.som, clust$cluster)

pretty_palette <- c('#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#E377C2', '#99FFFF', '#FFCCCC', '#FFFF00')
coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}

fviz_nbclust(class$codes[[1]], kmeans, method = 'wss')
c.class <- kmeans(class$codes[[2]], 9)
par(mfrow = c(1,2))
plot(class, type = 'codes', main = c('Unsupervised SOM', 'Supervised SOM'), bgcol = pretty_palette[c.class$cluster],
     palette.name = coolBlueHotRed)#palette.name = pretty_palette)
add.cluster.boundaries(class, c.class$cluster)

#### PAULS C calcs

# Must first load project R data file
BeginYear = 2010
EndYear = 2013

epiV = odem_stan$volume_epi # m3
hypoV = odem_stan$volume_hypo # m3
totalV = odem_stan$volume_total # m3
Stratified = odem_stan$strat
LakeArea = odem_stan$area_surface # m2
Year = odem_stan$year
DayOfYear = odem_stan$day_of_year

iStrat = which(Stratified==1)
iNotStrat = which(Stratified==0)
Fatm = rep(NA,dim(odem_stan)[1])
Fnep = rep(NA,dim(odem_stan)[1])
FentrEpi = rep(NA,dim(odem_stan)[1])
FentrHypo = rep(NA,dim(odem_stan)[1])
Fmin = rep(NA,dim(odem_stan)[1])
Fsed = rep(NA,dim(odem_stan)[1])

#####################
# Convert fluxes to areal (lake area)

Fatm[iStrat] = odem_stan$Fatm[iStrat]/1000*epiV[iStrat]/LakeArea[iStrat]
Fatm[iNotStrat] = odem_stan$Fatm[iNotStrat]/1000*totalV[iNotStrat]/LakeArea[iNotStrat]
Fnep[iStrat] = odem_stan$Fnep[iStrat]/1000*epiV[iStrat]/LakeArea[iStrat]
Fnep[iNotStrat] = odem_stan$Fnep[iNotStrat]/1000*totalV[iNotStrat]/LakeArea[iNotStrat]
FentrEpi[iStrat] = odem_stan$Fentr1[iStrat]/1000*epiV[iStrat]/LakeArea[iStrat]
FentrEpi[iNotStrat] = odem_stan$Fentr1[iNotStrat]/1000*totalV[iNotStrat]/LakeArea[iNotStrat]

Fmin[iStrat] = odem_stan$Fmin[iStrat]/1000*hypoV[iStrat]/LakeArea[iStrat]
Fmin[iNotStrat] = odem_stan$Fmin[iNotStrat]/1000*totalV[iNotStrat]/LakeArea[iNotStrat]
Fsed[iStrat] = odem_stan$Fsed2[iStrat]/1000*hypoV[iStrat]/LakeArea[iStrat]
Fsed[iNotStrat] = odem_stan$Fsed2[iNotStrat]/1000*totalV[iNotStrat]/LakeArea[iNotStrat]
FentrHypo[iStrat] = odem_stan$Fentr2[iStrat]/1000*hypoV[iStrat]/LakeArea[iStrat]
FentrHypo[iNotStrat] = odem_stan$Fentr2[iNotStrat]/1000*totalV[iNotStrat]/LakeArea[iNotStrat]

#####################
# Calculate cumulatives

# Calculate the yearly cumulative fluxes
uYears = unique(Year)
yFatm = rep(NA,length(uYears))
yFnep = rep(NA,length(uYears))
yFentrEpi = rep(NA,length(uYears))
yFmin = rep(NA,length(uYears))
yFsed = rep(NA,length(uYears))
yFentrHypo = rep(NA,length(uYears))
yFnepTot = rep(NA,length(uYears))
for (i in 1:length(uYears)){
  thisYear = uYears[i]
  iYear = which(Year==thisYear)
  yFatm[i] = sum(Fatm[iYear])
  yFnep[i] =  sum(Fnep[iYear])
  yFentrEpi[i] =  sum(FentrEpi[iYear])
  yFmin[i] =  sum(Fmin[iYear])
  yFsed[i] =  sum(Fsed[iYear])
  yFentrHypo[i] =  sum(FentrHypo[iYear])
  yFnepTot[i] = yFnep[i]+yFmin[i]+yFsed[i]
}

# Subset data by year indices
iYears = which(Year>=BeginYear & Year<=EndYear) # For all days in all years
iYearsOnly = which(uYears>=BeginYear & uYears<=EndYear) # Just the years

# Cumulative fluxes
FatmCum = cumsum(Fatm[iYears])
FnepCum = cumsum(Fnep[iYears])
FminCum = cumsum(Fmin[iYears])
FsedCum = cumsum(Fsed[iYears])
FAllCum = FatmCum+FnepCum+FminCum+FsedCum
FnepTotCum = FnepCum+FminCum+FsedCum

#####################
# Plots

# For plotting, setup year fractions
YearFrac = Year+(DayOfYear/366)

# Plot Fatm and NEP time series
plot(YearFrac[iYears],Fatm[iYears],type='l',
     xlab = 'Year',ylab = 'Flux (g/m2/d)')
abline(h=0,lty=2)
lines(YearFrac[iYears],Fnep[iYears],col='green',)
lText = c("Fatm","epiNEP")
legend('topleft',lText,cex=0.7,lty=c(1,1,1,1,1),lwd=c(1,1,1,1,1),
       col=c('black','green','red','blue','grey'))

# Plot Fatm against NEP
plot(Fnep[iYears],-Fatm[iYears],
     xlab = 'Fnep (g/m2/d)',ylab = '(-) Fatm (g/m2/d)')
abline(a=0,b=1,lty=2)
abline(lm(-Fatm[iYears]~Fnep[iYears]), col="red") # regression line (y~x)

# Plot yearly total fluxes
maxY = max(c(yFatm,yFnep,yFmin))
myYLim = c(-maxY,maxY)
plot(uYears[iYearsOnly],yFatm[iYearsOnly],type='l',ylim = myYLim,
     xlab = 'Year',ylab = 'Total yearly flux (g/m2/y)')
abline(h=0,lty=2)
lines(uYears[iYearsOnly],yFnep[iYearsOnly],col='green')
lines(uYears[iYearsOnly],yFmin[iYearsOnly],col='red')
lines(uYears[iYearsOnly],yFsed[iYearsOnly],col='blue')
lines(uYears[iYearsOnly],yFnepTot[iYearsOnly],col='grey',lwd=2)

lText = c("Fatm","epiNEP","Min","Sed","totNEP")
legend('topleft',lText,cex=0.7,lty=c(1,1,1,1,1),lwd=c(1,1,1,1,1),
       col=c('black','green','red','blue','grey'))

# Plot total cumulative fluxes
maxY = max(c(FatmCum,FnepCum,FminCum,FsedCum))
myYLim = c(-maxY,maxY)
plot(YearFrac[iYears],FatmCum,type='l',ylim = myYLim,
     xlab = 'Year',ylab = 'Cumultative flux (g/m2)')
abline(h=0,lty=2)
lines(YearFrac[iYears],FnepCum,col='green')
lines(YearFrac[iYears],FminCum,col='red')
lines(YearFrac[iYears],FsedCum,col='blue')
lines(YearFrac[iYears],FnepTotCum,col='grey',lwd=2)
grid(NA, 5, lwd = 1)

lText = c("Fatm","epiNEP","Min","Sed","totNEP")
legend('topleft',lText,cex=0.7,lty=c(1,1,1,1,1),lwd=c(1,1,1,1,1),
       col=c('black','green','red','blue','grey'))

# Convert to annual carbon units
TYears = EndYear-BeginYear+1
Oxygen2Carbon = 1/32 * 12/1
meanFatmC = -FatmCum[length(FatmCum)]/TYears * Oxygen2Carbon # gC/m2/y
meanFnepC = FnepCum[length(FnepCum)]/TYears * Oxygen2Carbon # gC/m2/y
meanFnepTotC = FnepTotCum[length(FnepTotCum)]/TYears * Oxygen2Carbon # gC/m2/y
meanFminC = FminCum[length(FminCum)]/TYears * Oxygen2Carbon # gC/m2/y
meanFsedC = FsedCum[length(FsedCum)]/TYears * Oxygen2Carbon # gC/m2/y
print('*************')
print('Mean annual C units')
print(paste('Fatm: ',signif(meanFatmC,3),' (gC/m2/y)',sep=""))
print(paste('Fnep: ',signif(meanFnepC,3),' (gC/m2/y), epilimnion',sep=""))
print(paste('Fmin: ',signif(meanFminC,3),' (gC/m2/y)',sep=""))
print(paste('Fsed: ',signif(meanFsedC,3),' (gC/m2/y)',sep=""))
print('-------------')
print(paste('Fnep total: ',signif(meanFnepTotC,3),' (gC/m2/y), fate as burial+export',sep=""))
print('If we assume approximate steady state for water column OC,')
print('then FnepTotal is ~ contribution of allocthony to total burial.')
print('Note that total burial+export is greater than FnepTotal,')
print('because burial+export also includes allocthony not mineralized.')

# Plot atmospheric flux plus other estimates of Fatm
# The O2 deficit is saturation - modeled, so if modeled>sat, the value is negative
# Scaler adjusts the deficit so that it can be viewed on the same y axis
TDepth = odem_stan$thermocline_depth
iZeros = which(TDepth == 0)
TDepth[iZeros] = max(TDepth) # set zeros to maximum TDepth
myK = 0.8 # m/d # piston velocity for atm exchange
Scaler = 1/8 # 1/m # mixed layer depth
EpiT = odem_stan$temperature_epi
EpiT[which(EpiT==0)] = odem_stan$temperature_total[which(EpiT==0)]
o2sat = o2.at.sat.base(temp = EpiT, altitude = 480) * 1000
o2def = myK * (o2sat - odem_stan$DO_epi) * 1/TDepth / 1000 # convert to mg/L

o2defObs = myK * (o2sat - odem_stan$DO_obs_epi) * 1/TDepth / 1000
plot(YearFrac[iYears], odem_stan$Fatm[iYears]/1000,type='l',ylim=c(-1,1),
     xlab='Day of sim',ylab='Atmospheric flux (mgO2/L/d)',main="")
lines(YearFrac[iYears],o2def[iYears],col='red')
points(YearFrac[iYears],o2defObs[iYears],col='blue')
abline(h=0)
legend('topright',col=c('black','red','blue'),lty=c(1,1,1),c('Modeled Fatm','Fatm, modeled O2 def scaled','Fatm, observed O2 def scaled'))


# Plot observations and fits in model output units
epiMod = odem_stan$DO_epi
epiObs = odem_stan$DO_obs_epi
hypoMod = odem_stan$DO_hyp
hypoObs = odem_stan$DO_obs_hyp

plot(YearFrac[iYears],epiMod,type='l',col='red',
     xlab = 'Year fraction',ylab='DO (mg/m3)')
points(YearFrac[iYears],epiObs,col='red')
lines(YearFrac[iYears],hypoMod,col='green')
points(YearFrac[iYears],hypoObs,col='green')
legend('topright',col=c('red','green'),lty=c(1,1,1),c('Epi O2','Hypo O2'))

# Plot time series of modeled, observed, and saturated epi DO
plot(YearFrac[iYears],epiMod,ylim=c(6000,13000),type='l',col='red',
     xlab = 'Year fraction',ylab='DO (mg/m3)')
lines(YearFrac[iYears],o2sat,type='l',col='blue')
points(YearFrac[iYears],epiObs,col='red')
legend('topright',col=c('red','red','blue'),lty=c(1,NaN,1),pch=c(NaN,21,NaN),c('Modeled O2','Obs O2','O2sat'))




