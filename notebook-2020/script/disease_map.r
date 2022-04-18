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

data <- read_rdump("data/disease_data_100.r")

modelName <- "disease_map_ela.stan"
mod <- cmdstan_model(paste0("model/", modelName))

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














