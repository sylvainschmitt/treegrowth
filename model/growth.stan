data {
  int<lower=1> N ; // # of observations
  int<lower=1> I ; // # of individuals
  int<lower=1> S ; // # of individuals
  int<lower=1> Y ; // max # of years
  array[N] int<lower=1> year ; // years
  vector[N] dbh ; // diameter
  array[N] int<lower=1, upper=I> ind ; // individuals
  array[I] int<lower=1, upper=S> indsp ; // species corresponding to individual
}
parameters {
  vector<lower=0.001, upper=3>[I] gmax ;
  vector<lower=0.001, upper=3>[I] dopt ;
  vector<lower=0.001, upper=3>[I] ks ;
  vector<lower=0.001, upper=3>[S] dopt_s ;
  vector<lower=0.001, upper=3>[S] ks_s ;
  real<lower=0> sigma ;
  real<lower=0> sigmaD ;
  real<lower=0> sigmaK ;
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
  dopt ~ normal(dopt_s[indsp], sigmaD) ;
  ks ~ normal(ks_s[indsp], sigmaK) ;
  gmax ~ normal(0, 1) ;
  dopt ~ normal(0, 1) ;
  ks ~ normal(0, 1) ;
  sigma ~ normal(0, 1) ;
  sigmaD ~ normal(0, 1) ;
  sigmaK ~ normal(0, 1) ;
}
