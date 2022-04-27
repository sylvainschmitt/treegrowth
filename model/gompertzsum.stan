functions {
  real dbhtoday(int years, real gmax, real dopt, real ks) {
    real agr = 0 ;
    real dbh = 10 ;
    for(t in 1:years) {
      agr = gmax*exp(-0.5* square(log(dbh / (100*dopt)) / ks)) ;
      dbh += agr ;
    }
    return dbh ;
  }
}
data {
  int<lower=1> N ; // # of observations
  int<lower=1> I ; // # of individuals
  int<lower=1> year[N] ; // years
  vector[N] dbh ; // diameter
  int<lower=1, upper=I> ind[N] ; // individuals
  int<lower=1> Np ; // # of observations for prediction
  vector[Np] dbhp ; // diameter diameter for prediction
  int<lower=1, upper=I> indp[Np] ; // individuals for prediction
  int<lower=1> yearp[Np] ; // years for prediction
}
parameters {
  vector<lower=0.001, upper=3>[I] gmax ;
  vector<lower=0.001, upper=3>[I] dopt ;
  vector<lower=0.001, upper=3>[I] ks ;
  real<lower=0> sigma ;
}
transformed parameters {
  vector[N] mu ;
  for(n in 1:N)
    mu[n] = dbhtoday(year[n], gmax[ind[n]], dopt[ind[n]], ks[ind[n]]) ;
}
model {
  dbh ~ normal(mu, sigma) ;
  gmax ~ normal(0, 1) ;
  dopt ~ normal(0, 1) ;
  ks ~ normal(0, 1) ;
  sigma ~ normal(0, 1) ;
}
generated quantities {
  vector[N] SEP ;
  real RMSEP ;
  vector[N] log_lik ;
  for(n in 1:Np)
    SEP[n] = square(dbhp[n] - dbhtoday(yearp[n], gmax[indp[n]], dopt[indp[n]], ks[indp[n]])) ;
  RMSEP = sqrt(mean(SEP)) ;
  for(n in 1:N)
    log_lik[n] = normal_cdf(dbh[n], mu[n], sigma) ;
}
