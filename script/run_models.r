# Reproduce disease map model, as developed by Aki and colleagues

rm(list = ls())
gc()
set.seed(1954)

# Adjust to your setting
setwd("~/Desktop/Code/laplace_approximation/writeup/StanCon2020/script")
.libPaths("~/Rlib")

# make sure cmdstanr is linked to cmdstan with the laplace functions
# installed. See install.sh.
library(cmdstanr)
set_cmdstan_path(file.path(getwd(), "cmdstan"))
library(parallel)

# We will use some of rstan's features, but we won't use it
# to fit the model
library(rstan)

#####################################################################
## Disease map

data <- read_rdump("data/disease_100.data.r")

modelName <- "disease_map_ela"
mod <- cmdstan_model(paste0("model/", modelName, ".stan"))

num_chains <- 4
num_cores <- min(num_chains, detectCores())

fit <- mod$sample(
  data = data, chains = num_chains,
  parallel_chains = num_cores,
  iter_warmup = 500, iter_sampling = 500, seed = 123)


theta_mean <- fit$summary()[6:105, 2]$mean

plot_data <- data.frame(x1 = data$x[, 1],
                        x2 = data$x[, 2],
                        theta_mean = theta_mean)

plot <- ggplot(data = plot_data,
               aes(x = x1, y = x2, color = theta_mean)) +
  geom_point() + theme_bw() +
  scale_color_gradient2(low = "black", mid = "blue", high = "red")
plot


# for convenience, let's convert the fit object to an
# rstan object.
stan_fit <- read_stan_csv(fit$output_files())
theta <- rstan::extract(stan_fit, pars = c("theta"))$theta


#####################################################################
## SKIM

source("tools.r")

data <- read_rdump("data/prostate_200.data.r")

modelName <- "skim_logit_ela"
mod <- cmdstan_model(paste0("model/", modelName, ".stan"))

num_chains <- 4
num_cores <- 4

fit <- mod$sample(
  data = data, chains = num_chains,
  parallel_chains = num_cores,
  iter_warmup = 1000, iter_sampling = 2000, seed = 123
)

pars = c("lp__", "eta_two", "tau", "lambda[1]", "lambda[2]")
# pars = c("p[1]", "p[2]")

fit$summary(pars)
stan_fit <- read_stan_csv(fit$output_files())

saveRDS(stan_fit, paste0("saved_fit/", modelName, ".RSave"))

# in case one iteration returns an NA.
p <- na.omit(extract(stan_fit, pars = c("p"))$p)

# plot the 90th quantiles.
quant = 0.9
lambda <- rstan::extract(stan_fit, pars = c("lambda"))$lambda
log_lambda <- log(lambda)
quant_select_plot(log_lambda, quant, threshold = 2.4) +
  ylab("90th quantile for \n log lambda")

select_lambda(log_lambda, quant, n = 6)


## Fit the model with full HMC
modelName <- "skim_logit"
mod <- cmdstan_model(paste0("model/", modelName, ".stan"))

fit <- mod$sample(
  data = data, chains = num_chains,
  parallel_chains = num_cores,
  iter_warmup = 1000, iter_sampling = 2000, seed = 123,
  adapt_delta = 0.99
)

pars = c("lp__", "eta_two", "tau", "lambda[1]", "lambda[2]")
fit$summary(pars)

stan_fit2 <- read_stan_csv(fit$output_files())
saveRDS(stan_fit2, paste0("saved_fit/", modelName, ".RSave"))

# stan_read <- readRDS(paste0("saved_fit/", modelName, ".RSave"))

p2 <- na.omit(extract(stan_fit2, pars = c("p"))$p)

# Compare the log lambdas
log_lambda2 <- log(extract(stan_fit2, pars = c("lambda"))$lambda)

quant_select_plot2(log_lambda2, log_lambda, quant, threshold = 2.4,
                   alpha = 0.5) +
  ylab("90th quantile for \n log lambda")

# plot the mean probability
p_mean <- colMeans(p)
p2_mean <- colMeans(p2)
plot_data <- data.frame(p_mean, p2_mean)

plot <- ggplot(data = plot_data, aes(x = p2_mean, y = p_mean)) +
  geom_point(size = 1) + theme_bw() +
  geom_abline(intercept = 0, slope = 1, 
              color = "red", 
              linetype = "dashed", size = 0.5) +
  xlim(0, 1) + ylim(0, 1) + xlab("Probability (full HMC)") +
  ylab("Probability (HMC + Laplace)") +
  theme(text = element_text(size = 15))
plot

