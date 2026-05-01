// what if smoothing, but in Stan?
data {
  // dimensions
  // number of data
  int<lower=0> N;
  // number of preds
  int<lower=0> Np;
  // number of non-hyperparameters
  int<lower=0> ncoef;

  // response
  array[N] int<lower=0, upper=1> y;
  // design matrix
  matrix[N, ncoef] X;
  // design matrix for predictions
  matrix[Np, ncoef] Xp;
  // coefficients
  vector[ncoef] betas;
  // variance matrix
  matrix[ncoef, ncoef] V;
}

parameters {
  // coefs (all together)
  vector[ncoef] b;
}

transformed parameters{

  // calculate the linear predictor
  vector[N] mu;
  mu = X*b;

  // for the predictions
  vector[Np] mup;
  mup = Xp*b;

}

model {

  // generate the coefficients for our model
  // sample from prior of the spline coefficients for this term
  b ~ multi_normal(betas, V);

  // likelihood!
  y ~ bernoulli_logit(mu);
}

generated quantities {
  vector[Np] y_rep;

  y_rep = inv_logit(mup);
}
