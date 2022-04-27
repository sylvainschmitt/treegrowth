data {
  int<lower=1> N ; // # of observations
  int<lower=1> I ; // # of individuals
  int<lower=1> Y ; // # of years
  vector<lower=1, upper=Y>[N] year ; // years
  int<lower=1, upper=I> ind[N] ; // individuals
  vector[N] dbh ; // diameter
  int<lower=1> Np ; // # of observations for prediction
  vector[Np] dbhp ; // diameter diameter for prediction
  int<lower=1, upper=I> indp[Np] ; // individuals for prediction
  vector<lower=1, upper=Y>[Np] yearp ; // years for prediction
}
parameters {
  vector<lower=0, upper=200>[I] alpha ;
  vector<lower=0, upper=10>[I] gamma ;
  real<lower=0> sigma ;
}
transformed parameters {
  vector[N] mu = 10 + 1 ./ (1 ./ alpha[ind] + 1 ./ (gamma[ind] .* year)) ;
}
model {
  dbh ~ normal(mu, sigma) ;
  alpha ~ normal(0, 100) ;
  gamma ~ normal(0, 10) ;
  sigma ~ normal(0, 1) ;
}
generated quantities {
  real RMSEP = sqrt(mean(square(dbhp - (10 + 1 ./ (1 ./ alpha[indp] + 1 ./ (gamma[indp] .* yearp)))))) ;
  vector[N] log_lik ;
  for(n in 1:N)
    log_lik[n] = normal_cdf(dbh[n], mu[n], sigma) ;
}
