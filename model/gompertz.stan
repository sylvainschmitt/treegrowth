data {
  int<lower=1> N ; // # of observations
  int<lower=1> I ; // # of individuals
  vector[N] dbhprev ; // diameter
  vector[N] dt ; // diameter
  vector[N] dbh ; // diameter
  int<lower=1, upper=I> ind[N] ; // individuals
  int<lower=1> Np ; // # of observations for prediction
  vector[Np] dbhp ; // diameter diameter for prediction
  int<lower=1, upper=I> indp[Np] ; // individuals for prediction
  vector[Np] dtp ;
  vector[Np] dbhprevp ;
}
parameters {
  vector<lower=0.001, upper=3>[I] gmax ;
  vector<lower=0.001, upper=3>[I] dopt ;
  vector<lower=0.001, upper=3>[I] ks ;
  real<lower=0> sigma ;
}
transformed parameters {
  vector[N] mu = dbhprev + gmax[ind] .* exp(-0.5* square(log(dbhprev ./ (100*dopt[ind])) ./ ks[ind])) .* dt ;
}
model {
  dbh ~ normal(mu, sigma) ;
  gmax ~ normal(0, 1) ;
  dopt ~ normal(0, 1) ;
  ks ~ normal(0, 1) ;
  sigma ~ normal(0, 1) ;
}
generated quantities {
  real RMSEP = sqrt(mean(square(dbhp - dbhprevp - gmax[indp] .* exp(-0.5* square(log(dbhprevp ./ (100*dopt[indp])) ./ ks[indp]) .* dtp)))) ;
  vector[N] log_lik ;
  for(n in 1:N)
    log_lik[n] = normal_cdf(dbh[n], mu[n], sigma) ;
}
