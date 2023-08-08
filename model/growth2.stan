data {
  int<lower=1> N ; // # of observations
  int<lower=1> I ; // # of individuals
  int<lower=1> S ; // # of species
  int<lower=1> G ; // # of genera
  int<lower=1> F ; // # of families
  int<lower=1> Y ; // max # of years
  array[N] int<lower=1> year ; // years
  vector[N] dbh ; // diameter
  vector[I] nci ; // neighbourhood
  vector[I] twi ; // topography
  vector[S] dmax ; // diameter
  array[N] int<lower=1, upper=I> ind ; // individuals
  array[I] int<lower=1, upper=S> indsp ; // species corresponding to individual
  array[I] int<lower=1, upper=G> indgen ; // genus corresponding to individual
  array[I] int<lower=1, upper=F> indfam ; // family corresponding to individual
}
parameters {
  real<lower=0.001, upper=10> gmax ;
  vector<lower=0.1, upper=1>[I] d ;
  vector<lower=0.001, upper=3>[I] ks ;
  vector<lower=0.1, upper=1>[S] d_s ;
  vector<lower=0.001, upper=3>[S] ks_s ;
  vector[2] beta ;
  real<lower=0.001> sigma ;
  real<lower=0.001> sigmaD ;
  real<lower=0.001> sigmaK ;
  real<lower=0.001> sigma_species ;
  real<lower=0.001> sigma_genera ;
  real<lower=0.001> sigma_families ;
  vector[S] epsilon_species ;
  vector[G] epsilon_genera ;
  vector[F] epsilon_families ;
}
transformed parameters {
  vector[I] dopt = d .* dmax[indsp] ;
  vector[S] dopt_s = d_s .* dmax ;
  vector[I] log_gmax_i = log(gmax) + 
                          epsilon_families[indfam]*sigma_families  + 
                          epsilon_genera[indgen]*sigma_genera  + 
                          epsilon_species[indsp]*sigma_species + 
                          beta[1]*log(nci) + beta[2]*log(twi+1) ;
  matrix<lower=10>[I,Y] DBH ;
  vector[N] mu ;
  DBH[,1] = rep_vector(10, I) ;
  for(t in 2:Y)
    DBH[,t] = DBH[,t-1] + exp(log_gmax_i) .* exp(-0.5 * square(log(DBH[,t-1] ./ dopt) ./ ks)) ;
  for(n in 1:N)
    mu[n] = DBH[ind[n],year[n]] ;
}
model {
  dbh ~ normal(mu, sigma) ;
  epsilon_species ~ std_normal() ;
  epsilon_genera ~ std_normal() ;
  epsilon_families ~ std_normal() ;
  d ~ normal(d_s[indsp], sigmaD) ;
  ks ~ normal(ks_s[indsp], sigmaK) ;
  gmax ~ normal(0, 1) ;
  d ~ normal(0.4, 1) ;
  ks ~ normal(0, 1) ;
  beta ~ normal(0, 1) ;
  sigma ~ normal(0, 1) ;
  sigmaD ~ normal(0, 1) ;
  sigmaK ~ normal(0, 1) ;
  sigma_species ~ normal(0, 1) ;
  sigma_genera ~ normal(0, 1) ;
  sigma_families ~ normal(0, 1) ;
}
generated quantities {
  real Vt = variance(log_gmax_i) ;
  real Vf = square(sigma_families) ;
  real Vg = square(sigma_genera) ;
  real Vs = square(sigma_species) ;
  real Ve = variance(beta[1]*log(nci) + beta[2]*log(twi+1)) ;
  real Vr = Vt - Vf - Vg - Vs - Ve ;
}
