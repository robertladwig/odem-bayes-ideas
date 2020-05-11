cat('\f')
rm(list= ls())

setwd('/home/robert/Projects/DSI/odem-bayes-ideas')

library(tidyverse)

oneyear <- 2008
twoyear <- 2007:2008

input <- readr::read_csv(
  'inst/extdata/input.txt',
  col_names=c('datetime', 'thermocline_depth', 'temperature_epi', 'temperature_hypo', 'temperature_total', 'volume_total', 'volume_epi', 'volume_hypo', 'area_thermocline', 'area_surface', 'year', 'day_of_year', 'wind'),
  col_types=cols(datetime=col_datetime(), year=col_integer(), day_of_year=col_integer(), .default=col_double()))
in1yr <- filter(input, year %in% oneyear)
in2yr <- filter(input, year %in% twoyear)

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
  'inst/extdata/observed.txt',
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
ggplot(obs, aes(x=date)) +
  geom_line(aes(y=DO_tot), color='black') +
  geom_line(aes(y=DO_epi), color='seagreen') +
  geom_line(aes(y=DO_hypo), color='navy') +
  theme_bw()


idx1 = which(!is.na(in1yr$thermocline_depth))[1]
idx2 = rev(which(!is.na(in1yr$thermocline_depth)))[1]
in1yr = in1yr[idx1:idx2,]

idy = match(obs1yr$date,zoo::as.Date(in1yr$datetime))
idx = idy[!is.na(idy)]
idz = which(!is.na(idy))
DO_obs_epi = rep(NA, length(in1yr$datetime))
DO_obs_epi[idx] = obs1yr$DO_epi[idz]
DO_obs_hyp = rep(NA, length(in1yr$datetime))
DO_obs_hyp[idx] = obs1yr$DO_hypo[idz]

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(LakeMetabolizer)
# we started with 10 days but I've played with 200 to get more "data"

simdata <- tibble(
  DO_obs_epi = DO_obs_epi * 1000,
  DO_obs_hyp = DO_obs_hyp * 1000,
  day = seq(1, length(in1yr$datetime))
  # DO_obs_tot_true = 10 * exp(-0.04*day) + rnorm(length(day), 0, 0.0002),
  # have_obs = ifelse(day == 1, 0, round(runif(10))),
  # DO_obs_tot = ifelse(have_obs == 1, DO_obs_tot_true, NA)
)
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
  theta1 = 1.08^(in1yr$temperature_epi - 20),
  theta2 = 1.08^(in1yr$temperature_hypo - 20),
  k600 = k600.2.kGAS.base(k.cole.base(in1yr$wind),temperature = in1yr$temperature_epi, gas = "O2"),
  o2sat = o2.at.sat.base(temp = in1yr$temperature_epi, altitude = 300) * 1000,
  volume_epi = in1yr$volume_epi,
  area_epi = in1yr$area_surface,
  volume_hyp = in1yr$volume_hypo,
  area_hyp = in1yr$area_thermocline,
  tddepth = in1yr$thermocline_depth,
  ii_obs = idx,
  wtr_epi = in1yr$temperature_epi,
  wtr_hyp = in1yr$temperature_hypo,
  khalf = 3000,
  err_sigma = 0.0003,
  d = nrow(simdata),
  DO_epi_init = 10 * 1000, #simdata$DO_obs[1],
  DO_hyp_init = 10 * 1000
)

dummyinput$N_obs = length(dummyinput$DO_epi_init)
dummyinput$DO_obs_epi = simdata$DO_obs_epi[idx]
dummyinput$DO_obs_hyp = simdata$DO_obs_hyp[idx]
dummyinput$N_obs = length(dummyinput$ii_obs)

fit <- stan(file = 'src/odem.stan', data = dummyinput, chains = 3, iter = 500)

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
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = SED1, aes(x=day, y=mean, col = 'SED1')) +
  geom_ribbon(data = SED1, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = SED2, aes(x=day, y=mean, col = 'SED2')) +
  geom_ribbon(data = SED2, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = MIN, aes(x=day, y=mean, col = 'MIN')) +
  geom_ribbon(data = MIN, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) + theme_bw()
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
  ylim(c(0,15000))

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
ggplot(Fnep, aes(x=day, y=mean, col = 'NEP')) + geom_line() +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fsed1, aes(x=day, y=mean, col = 'SED1')) +
  geom_ribbon(data = Fsed1, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fatm, aes(x=day, y=mean, col = 'ATM')) +
  geom_ribbon(data = Fatm, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = DO_epi, aes(x=day, y=mean, col = 'epi')) + geom_line() +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) + theme_bw()

ggplot(Fsed2, aes(x=day, y=mean, col = 'SED2')) + geom_line() +
  geom_ribbon(data = Fsed2, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = Fmin, aes(x=day, y=mean, col = 'MIN')) +
  geom_ribbon(data = Fmin, aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) +
  geom_line(data = DO_hyp, aes(x=day, y=mean, col = 'hyp')) + geom_line() +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) + theme_bw()

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
