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

  fit <- stan("binomial_sampler.stan", model_name="sampler",
              data = dat.Stan.fit,
              algorithm="HMC",
              iter = iter, chains = chains, cores = cores, warmup = warmup,
              control = list(adapt_delta = adapt_delta,
                             max_treedepth = max_treedepth),
              pars = c("y_rep"))

  # response matrix on logit scale
  as.matrix(fit, pars = "y_rep")
}
