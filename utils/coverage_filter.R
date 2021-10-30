N <- 593

x <- seq(0, 10, length.out = 5000)
 
plot(x, dlnorm(x, log(20), log(10)), type = 'l')

mean(log(rlnorm(5000, log(20), log(10))))

#remove initial individual filter - let it get caught up in next filter


#Goal - figure out how may individuals will have <3 reads with a given minimum mean coverage at a locus

library(tidyverse)
library(magrittr)
library(fitdistrplus)

mean_depth <- read_lines('../tmp_dir/MiSeq_CHYA_chyaI.2.1.Fltr04.6.site.depth.mean') %>%
  extract(-1) %>%
  as.numeric()

plotdist(mean_depth, histo = TRUE, demp = TRUE)
descdist(mean_depth, discrete=FALSE, boot=500)

lnorm_fit <- fitdist(mean_depth, "lnorm")
gamma_fit <- fitdist(mean_depth, "gamma")


par(mfrow=c(2,2))
plot.legend <- c("lognormal", "gamma")
denscomp(list(lnorm_fit, gamma_fit), legendtext = plot.legend)
cdfcomp(list(lnorm_fit, gamma_fit), legendtext = plot.legend)
qqcomp(list(lnorm_fit, gamma_fit), legendtext = plot.legend)
ppcomp(list(lnorm_fit, gamma_fit), legendtext = plot.legend)
par(mfrow=c(1,1))

gofstat(list(lnorm_fit, gamma_fit), fitnames = c("lnorm", "gamma"))

#Log Normal
lnorm_param <- bootdist(lnorm_fit, niter = 1e3, parallel = 'snow', ncpus = 3)
summary(lnorm_param)
plot(lnorm_param)



#https://stats.stackexchange.com/questions/95498/how-to-calculate-log-normal-parameters-using-the-mean-and-std-of-the-given-distr
library(multidplyr)
library(broom)


tibble(x = rlnorm(5000, 2.20, 1.52)) %>%
  filter(x < 50) %>%
  ggplot(aes(x = x)) +
  geom_histogram()

cluster <- new_cluster(4)
cluster_library(cluster, 'broom')

out <- read_delim('../tmp_dir/meandepthVSvariance', delim = '\t',
           col_names = c('mean', 'variance'),
           col_types = cols(.default = col_number())) %>%
  mutate(sd = sqrt(variance)) %>%
  dplyr::select(-variance) %>%
  dplyr::rename(m = mean, s = sd) %>%
  mutate(id = str_c('R', row_number()),
         l_m = log(m) - 0.5 * log((s/m)^2 + 1),
         l_s = sqrt(log((s/m)^2 + 1)),
         tries = 100) %>%
  rowwise %>%
  partition(cluster) %>%
  mutate(n_less = sum(rlnorm(tries, l_m, l_s) <= 3),
         tidy(binom.test(x = n_less, n = tries))) %>%
  collect

out %>%
  filter(m > 100) %>%
  sample_n(500) %>%
  ggplot(aes(x = m, y = n_less / tries)) +
  geom_point()

library(emmeans)
library(brms)

filter(out, !(n_less == 0 | n_less == tries))



tst <- glm(cbind(n_less, tries - n_less) ~ m * s, 
           data = mutate(out, across(c(n_less, tries), ~. + 1)),
           family = 'binomial') 
summary(tst)


get_prior(n_less | trials(tries) ~ m * s,
          family = 'binomial',
          data = out)

inits <- list(
  Intercept = -4.5e-1,
  beta     = c(-6e-2, 1e-2, 2e-5)
)

tst <- brm(n_less | trials(tries) ~ m * s,
           family = 'binomial',
           prior = prior(normal(-4.5e-1, 1e-3), class = 'Intercept') +
             prior(normal(-6e-2, 3e-5), class = 'b', coef = 'm') +
             prior(normal(1e-2, 1.5e-5), class = 'b', coef = 's') +
             prior(normal(2e-5, 9e-8), class = 'b', coef = 'm:s'),
           data = filter(out, !(n_less == 0 | n_less == tries)),
           
           inits = list(inits, inits),
           
           backend = 'cmdstanr',
           chains = 2,
           cores = 2,
           iter = 500, 
           warmup = 500)

library(lme4)
tst <- glmer(cbind(n_less, tries - n_less) ~ m * s + (1 | id), 
             data = mutate(out, across(c(n_less, tries), ~. + 1)),
             family = 'binomial')

tst %>%
  emmeans(., ~ m, by = 's', type = 'response',
          at=list(m = seq(1, max(out$m), length.out = 1000))) %>%
  tidy(conf.int = TRUE) %>%
  mutate(across(c(prob, asymp.LCL, asymp.UCL), ~.*593)) %>%
  ggplot(aes(x = m, y = prob, ymin = asymp.LCL, ymax = asymp.UCL)) +
  geom_point(data = sample_n(out, 1000), aes(x = m, y = n_less/tries * 593),
             inherit.aes = FALSE) +
  geom_ribbon() +
  geom_line(colour = 'red') +
  scale_x_continuous(limits = c(NA, NA)) +
  theme_classic()
  


