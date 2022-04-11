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
  real<lower=0, upper=1> theta2 ;
  real<lower=0, upper=1> theta3 ;
  vector[S] epsilon_s2 ;
  vector[S] epsilon_s3 ;
  vector<lower=0>[3] sigma ;
}
transformed parameters {
  vector<lower=0>[I] DBH = rep_vector(1, I) ;
  vector[S] thetas2 = exp(log(theta2) + sigma[2]*epsilon_s2) ; 
  vector[S] thetas3 = exp(log(theta3) + sigma[3]*epsilon_s3) ; 
  for(t in 1:Y-1) {
    for(i in 1:I) {
      if(years[t] == Y0[i])
        DBH[i] = DBH0[i] ;
    }
    DBH += thetas1[species] .* exp(-0.5* square(log(DBH ./ (100*thetas2[species])) ./ thetas3[species])) ;
  }
}
model {
  log(DBHtoday) ~ normal(log(DBH), sigma[1]) ;
  epsilon_s2 ~ std_normal() ;
  epsilon_s3 ~ std_normal() ;
  thetas1 ~ lognormal(0, 1) ;
  theta2 ~ lognormal(0, 1) ;
  theta3 ~ lognormal(0, 1) ;
  sigma ~ normal(0, 1) ;
}
generated quantities {
  vector[I] thetai1 = thetas1[species] + DBHtoday - DBH ; 
}
