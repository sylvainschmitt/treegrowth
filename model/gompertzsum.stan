data {
  int<lower=1> N ; // # of observations
  int<lower=1> I ; // # of individuals
  int<lower=1> Y ; // max # of years
  array[N] int<lower=1> year ; // years
  vector[N] dbh ; // diameter
  array[N] int<lower=1, upper=I> ind ; // individuals
  int<lower=1> Np ; // # of observations for prediction
  vector[Np] dbhp ; // diameter diameter for prediction
  array[Np] int<lower=1, upper=I> indp ; // individuals for prediction
  array[Np] int<lower=1> yearp ; // years for prediction
}
parameters {
  vector<lower=0.001, upper=3>[I] gmax ;
  vector<lower=0.001, upper=3>[I] dopt ;
  vector<lower=0.001, upper=3>[I] ks ;
  real<lower=0> sigma ;
}
transformed parameters {
  matrix<lower=10>[I,Y] DBH ;
  vector[N] mu ;
  DBH[,1] = rep_vector(10, I) ;
  for(t in 2:Y)
    DBH[,t] = DBH[,t-1] + gmax .* exp(-0.5 * square(log(DBH[,t-1] ./ (100*dopt)) ./ ks)) ;
  for(n in 1:N)
    mu[n] = DBH[ind[n],year[n]] ;
}
model {
  dbh ~ normal(mu, sigma) ;
  gmax ~ normal(0, 1) ;
  dopt ~ normal(0, 1) ;
  ks ~ normal(0, 1) ;
  sigma ~ normal(0, 1) ;
}
generated quantities {
  vector[Np] mup ;
  real RMSEP ;
  vector[N] log_lik ;
  for(n in 1:N)
    log_lik[n] = normal_cdf(dbh[n] | mu[n], sigma) ;
  for(n in 1:Np)
    mup[n] = DBH[indp[n],yearp[n]] ;
  RMSEP = sqrt(mean(square(dbhp - mup)));
}
