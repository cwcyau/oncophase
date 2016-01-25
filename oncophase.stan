data {
  int<lower=0> Q;
  real<lower=0> mu_h;
  real<lower=0> sigma_h;
  real<lower=0> e_a;
  real<lower=0> e_b;
  vector<lower=0>[Q] alpha;
  int<lower=0> k_snv;
  int<lower=0> d_snv;
  int<lower=0> k_snp;
  int<lower=0> d_snp;
  vector[Q] M_snv;
  vector[Q] C_snv;
  vector[Q] M_snp;
  vector[Q] C_snp;
}
parameters {
  real<lower=0,upper=1e6> h;
  real<lower=0,upper=1e6> sigma_d;
  simplex[Q] theta;
  real<lower=0,upper=1> e;
}
transformed parameters {
  
  real<lower=0,upper=1> phi_snv;
  real<lower=0,upper=1> phi_snp;
  
  real<lower=0,upper=1> p_snv;
  real<lower=0,upper=1> p_snp;
  
  real<lower=0> l_snv;
  real<lower=0> l_snp;
  
  phi_snv <- (M_snv'*theta)/(C_snv'*theta);
  phi_snp <- (M_snp'*theta)/(C_snp'*theta);
  l_snv <- h*(C_snv'*theta);
  l_snp <- h*(C_snp'*theta);
  p_snv <- (1-e)*phi_snv + e*(1-phi_snv);
  p_snp <- (1-e)*phi_snp + e*(1-phi_snp);
}
model {
  k_snv ~ binomial(d_snv, p_snv )  ;
  k_snp ~ binomial(d_snp, p_snp ) ;
  d_snv ~ normal(l_snv, sigma_d) ;
  d_snp ~ normal(l_snp, sigma_d) ;
  theta ~ dirichlet(alpha);
  h ~ normal(mu_h, sigma_h);
  e ~ beta(e_a, e_b);
}