data {
  int<lower=1> I ; // # of individuals
  int<lower=1> Y ; // # of census years
  int<lower=1> S ; // # of species
  vector[Y] years ; // years
  vector[I] DBH0 ; // recruitment DBH
  vector[I] Y0 ; // recruitment year
  int<lower=1, upper=S> species[I] ; // gene pools
  vector[I] DBHtoday ; // today DBH
}
parameters {
  matrix<lower=0, upper=1>[S,3] thetas ;
  real<lower=0> sigma ;
}
transformed parameters {
  vector<lower=0>[I] Sigma = rep_vector(1, I) ;
  for(t in 1:Y-1) {
    for(i in 1:I) {
      if(years[t] == Y0[i])
        Sigma[i] = DBH0[i] ;
    }
    Sigma += exp(-0.5* square(log(Sigma ./ (100*thetas[species,2])) ./ thetas[species,3])) ;
  }
  Sigma = Sigma - DBH0 ;
}
model {
  DBHtoday - DBH0 ~ lognormal(log(thetas[species,1] .* Sigma), sigma) ;
  for(j in 1:3)
    thetas[,j] ~ lognormal(0, 1) ;
  sigma ~ normal(0, 1) ;
}
generated quantities {
  vector[I] thetai = DBHtoday ./ Sigma ; 
}
