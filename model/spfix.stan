data {
  int<lower=1> N ; // # of observations
  int<lower=1> I ; // # of individuals
  int<lower=1> S ; // # of individuals
  int<lower=1> Y ; // max # of years
  array[N] int<lower=1> year ; // years
  vector[N] dbh ; // diameter
  array[N] int<lower=1, upper=I> ind ; // individuals
  array[N] int<lower=1, upper=S> sp ; // species
  array[I] int<lower=1, upper=S> indsp ; // species corresponding to individual
  int<lower=1> Np ; // # of observations for prediction
  vector[Np] dbhp ; // diameter diameter for prediction
  array[Np] int<lower=1, upper=I> indp ; // individuals for prediction
  array[Np] int<lower=1> yearp ; // years for prediction
}
parameters {
  vector<lower=0.001, upper=3>[I] gmax ;
  vector<lower=0.001, upper=3>[S] dopt_s ;
  vector<lower=0.001, upper=3>[S] ks_s ;
  real<lower=0> sigma ;
  real<lower=0, upper=1> sigmaD ;
  real<lower=0, upper=1> sigmaK ;
}
transformed parameters {
  matrix<lower=10>[I,Y] DBH ;
  vector[N] mu ;
  DBH[,1] = rep_vector(10, I) ;
  for(t in 2:Y)
    DBH[,t] = DBH[,t-1] + gmax .* exp(-0.5 * square(log(DBH[,t-1] ./ (100*dopt_s[indsp])) ./ ks_s[indsp])) ;
  for(n in 1:N)
    mu[n] = DBH[ind[n],year[n]] ;
}
model {
  dbh ~ normal(mu, sigma) ;
  gmax ~ normal(0, 1) ;
  dopt_s ~ normal(0, 1) ;
  ks_s ~ normal(0, 1) ;
  sigma ~ normal(0, 1) ;
  sigmaD ~ normal(0, 1) ;
  sigmaK ~ normal(0, 1) ;
}
generated quantities {
  vector[Np] mup ;
  real RMSEP ;
  vector[N] log_lik ;
  for(n in 1:N)
    log_lik[n] = normal_lpdf(dbh[n] | mu[n], sigma) ;
  for(n in 1:Np)
    mup[n] = DBH[indp[n],yearp[n]] ;
  RMSEP = sqrt(mean(square(dbhp - mup)));
}
