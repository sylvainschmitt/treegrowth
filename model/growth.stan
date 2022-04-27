data {
  int<lower=1> N ; // # of observations
  int<lower=1> I ; // # of individuals
  int<lower=1> S ; // # of species
  int<lower=1> Y ; // max # of years
  array[N] int<lower=1> year ; // years
  vector[N] dbh ; // diameter
  vector[I] dbh0 ; // diameter at recruitment
  array[N] int<lower=1, upper=I> ind ; // individuals
  array[N] int<lower=1, upper=S> sp ; // species
}
parameters {
  vector<lower=0.001, upper=3>[I] gmax ;
  vector<lower=0.001, upper=3>[I] dopt ;
  vector<lower=0.001, upper=3>[I] ks ;
  real<lower=0> sigma ;
}
transformed parameters {
  matrix[I,Y] DBH ;
  vector[N] mu ;
  DBH[,1] = dbh0 ;
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
