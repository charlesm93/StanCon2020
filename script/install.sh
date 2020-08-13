# install
# v 2.0
# Run this bash file the directory where you want to install cmdStan
# with prototype functions for the embedded Laplace approximation
# to run the code in Stan Con 2020 notebook:
# "Approximate Bayesian inference for latent Gaussian models in Stan".
#

# set your working directory

git clone https://github.com/stan-dev/stanc3.git
cd stanc3
git checkout try-laplace_approximation2
cd ..
git clone https://github.com/stan-dev/cmdstan.git
cd cmdstan
make stan-update
cd stan
cd lib/stan_math
git checkout try-laplace_approximation2
cd ..
cd ..
cd ..
# make build
STANC3=../stanc3 make examples/bernoulli/bernoulli
