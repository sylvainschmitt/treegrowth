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
  vector<lower=0, upper=1>[S] thetas1 ;
  vector<lower=0, upper=1>[S] thetas2 ;
  vector<lower=0, upper=1>[S] thetas3 ;
  real<lower=0> sigma ;
}
transformed parameters {
  vector<lower=0>[I] DBH = rep_vector(1, I) ;
  for(t in 1:Y-1) {
    for(i in 1:I) {
      if(years[t] == Y0[i])
        DBH[i] = DBH0[i] ;
    }
    DBH += thetas1[species] .* exp(-0.5* square(log(DBH ./ (100*thetas2[species])) ./ thetas3[species])) ;
  }
}
model {
  log(DBHtoday) ~ normal(log(DBH), sigma) ;
  thetas1 ~ lognormal(0, 1) ;
  thetas2 ~ lognormal(0, 1) ;
  thetas3 ~ lognormal(0, 1) ;
  sigma ~ normal(0, 1) ;
}
generated quantities {
  vector[I] thetai1 = thetas1[species] + DBHtoday - DBH ; 
}
