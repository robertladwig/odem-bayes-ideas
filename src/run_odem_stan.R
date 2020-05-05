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


#  time stamp in the first row, then area weighted average total oxygen in the second row, area weighted average epilimnion oxygen in the third row, and area weighted average hypolimnion oxygen in the fourth row
obs <- read.table(
  'inst/extdata/observed.txt',
  header=FALSE,
  sep=' ',
  as.is=TRUE) %>%
  t() %>%
  as_tibble(.name_repair='minimal') %>%
  setNames(., nm=c('dateint', 'DO_tot', 'DO_epi', 'DO_hypo')) %>%
  mutate(date = as.Date(dateint, origin='1970-01-01')) %>% # just guessing at origin and therefore at dates
  select(date, everything())
obs1yr <- filter(obs, lubridate::year(date) %in% oneyear)
obs2yr <- filter(obs, lubridate::year(date) %in% twoyear)
ggplot(obs, aes(x=date)) +
  geom_line(aes(y=DO_tot), color='black') +
  geom_line(aes(y=DO_epi), color='seagreen') +
  geom_line(aes(y=DO_hypo), color='navy') +
  theme_bw()

library(rstan)
# we started with 10 days but I've played with 200 to get more "data"
simdata <- tibble(
  day = 1:200,
  DO_obs_tot_true = 10 * exp(-0.04*day) + rnorm(length(day), 0, 0.0002),
  have_obs = ifelse(day == 1, 0, round(runif(10))),
  DO_obs_tot = ifelse(have_obs == 1, DO_obs_tot_true, NA)
)
dummyinput <- list(
  lambda_mu_min = 0,
  lambda_mu_max = 5,
  lambda_sigma = 0.000001,
  err_sigma = 0.0003,
  d = nrow(simdata),
  ii_obs = which(simdata$have_obs == 1),
  DO_tot_init = simdata$DO_obs_tot_true[1])

dummyinput$N_obs = length(dummyinput$ii_obs)
dummyinput$DO_obs_tot = simdata$DO_obs_tot[dummyinput$ii_obs]

fit <- stan(file = 'src/odem.stan', data = dummyinput, chains = 4, iter = 10000)

# an example of extracting parameters for this particular dummy model, i'm
# geting an overestimate of lambda and consequently a much faster modeled drop
# than actual drop in "DO". As an exercise (entirely optional), you might find
# it useful to explore different values of lambda hyperparameters (mu_min,
# mu_max, and sigma), data-generating parameters, MCMC iterations, etc. to see
# if you can achieve a better fit. Then again, it might be better just to jump
# straight into more realistic ODEM equations.
lambda <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$lambda %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
ggplot(lambda, aes(x=day, y=mean)) + geom_line() +
  geom_ribbon(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), alpha=0.2) + theme_bw()
DO_tot <- rstan::extract(fit, permuted = TRUE, inc_warmup=FALSE)$DO_tot %>%
  as_tibble() %>%
  mutate(iter = 1:n()) %>%
  pivot_longer(names_to='Vday', cols = -iter) %>%
  tidyr::extract(Vday, into='day', regex='V([[:digit:]]+)', convert=TRUE) %>%
  group_by(day) %>%
  summarize(mean = mean(value), sd = sd(value))
ggplot(DO_tot, aes(x=day, y=mean)) + geom_line() +
  geom_point(data=simdata, aes(x=day, y=DO_obs_tot))
