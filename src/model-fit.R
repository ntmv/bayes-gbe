# Libraries
library(rstan)

# STAN code from brms 
scode <- 
 "functions {
  }
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_shape;  // number of population-level effects
  matrix[N, K_shape] X_shape;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  int prior_only;  
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  int Kc_shape = K_shape - 1;
  matrix[N, Kc_shape] Xc_shape;  // centered version of X_shape without an intercept
  vector[Kc_shape] means_X_shape;  // column means of X_shape before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_shape) {
    means_X_shape[i - 1] = mean(X_shape[, i]);
    Xc_shape[, i - 1] = X_shape[, i] - means_X_shape[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=0>[Kc_shape] b_shape;  // population-level effects
  real Intercept_shape;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  r_1_1 = (sd_1[1] * (z_1[1]));
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    // initialize linear predictor term
    vector[N] shape = Intercept_shape + Xc_shape * b_shape;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
    for (n in 1:N) {
      // apply the inverse link function
      shape[n] = exp(shape[n]);
    }
    for (n in 1:N) {
      // apply the inverse link function
      mu[n] = shape[n] * exp(-(mu[n]));
    }
    target += gamma_lpdf(Y | shape, mu);
  }
  // priors including constants
  target += normal_lpdf(b | 0, 2);
  target += normal_lpdf(Intercept | 0, 5);
  target += normal_lpdf(b_shape | 1, 0.7)
  - 1 * normal_lccdf(0 | 1, 0.7);
  target += student_t_lpdf(Intercept_shape | 3, 0, 2.5);
  target += cauchy_lpdf(sd_1 | 1, 0.7)
  - 1 * cauchy_lccdf(0 | 1, 0.7);
  target += std_normal_lpdf(z_1[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_shape_Intercept = Intercept_shape - dot_product(means_X_shape, b_shape);
}"

# Specifying priors
priors <- c(prior(normal(0,5), class = Intercept),
            prior(cauchy(1, 0.7), class = sd), prior(normal(0,2), class = "b"), 
            prior(normal(1, 0.7), class = "b", dpar = "shape", lb = 0))


# Fit using brms 
fit1 <- brm(data = mouse_subsetCVR, family = Gamma(link = "log"),
                  bf(logCVR  ~ 1 + alleletype + zygosity + (1|study), 
                     shape ~ 1 + logCVR_var), iter = 30000, prior = priors,
                  warmup = 1000, cores = 4, chains = 4, seed = 111, inits  = 'random', 
                   file = "model_fit.rds") 

png(filename = "posterior_trace.png", height = 8, width = 12, res = 300, unit = "cm")
# Traceplots and posterior distributions
plot(fit1)
dev.off()

# Posterior predictive check
png(filename = "pp_check.png", height = 8, width = 12, res = 300, unit = "cm")
pp_check(fit1, ndraws = 100) + xlim(c(0, 5))
dev.off()

# Posterior predictive checks 
# y and yrep for more plots
y <- mouse_subsetCVR$logCVR
yrep <- posterior_predict(fit1, draws = 1000)

# LOO-CV
# Checking calibration by comparing transformed variables to uniform distribution 
model_loo <- loo(fit1, save_psis = TRUE)
w <- weights(model_loo$psis_object)
color_scheme_set("gray")
png(filename = "loo-pit.png", height = 8, width = 12, res = 300, unit = "cm")
ppc_loo_pit_overlay(y, yrep, lw = w)
dev.off()

