library(cmdstanr)
binomial_sampler <- function(model, data, iter=2000, chains=4, cores=4, warmup=1000, adapt_delta=0.95, max_treedepth=15){

  # generate prediction matrix
  Lp <- predict(model, data, type="lpmatrix")

  dat.Stan.fit <- list(betas=coef(model),
                       V=vcov(model),
                       X = model.matrix(model),
                       ncoef=length(coef(model)),
                       N=nrow(model.matrix(model)),
                       Xp = Lp,
                       Np = nrow(Lp),
                       y=model$y)

  mod <- cmdstan_model("binomial_sampler.stan")

  lapapp <- mod$laplace(data = dat.Stan.fit,draws = 500, refresh = 100,
                        threads=4)
  lapapp
}
