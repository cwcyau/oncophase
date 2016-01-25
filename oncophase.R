library(rstan)
library(gtools)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# define parameters 
Q <- 3 # number of cell type per configuration
idx <- 2 # configuration to simulate data from
n <- 10 # number of simulations to perform
mu_h <- 30 # mean haploid read coverage
sigma_h <- sqrt(mu_h) # std dev of haploid read coverage
e_a <- 0.5 # hyperparameters of beta prior on sequence read error parameter
e_b <- 0.5
alpha <- c(1, 1, 1) # hyperparameters of dirichlet prior on cellular prevalences

# read in M and C matrices
Mdat = read.table("example1-M.txt") 
Cdat = read.table("example1-C.txt") 

# extract number of configurations
S = dim(Cdat)[1]

# extract copy numbers for each configuration (SNP/SNV)
M_snp_table = as.matrix(Mdat[,1:3])
M_snv_table = as.matrix(Mdat[,4:6])
C_snp_table = as.matrix(Cdat[,1:3])
C_snv_table = as.matrix(Cdat[,4:6])

M_snp <- M_snp_table[idx, ]
M_snv <- M_snv_table[idx, ]

C_snp <- C_snp_table[idx, ]
C_snv <- C_snv_table[idx, ]


# simulate some random read counts

# preallocate vectors to store read counts
d_snp_all = rep(0, n) # total read count for SNP
d_snv_all = rep(0, n) # total read count for SNV
k_snp_all = rep(0, n) # total variant count for SNP
k_snv_all = rep(0, n) # total variant count for SNV

# generate n observations
theta = matrix(0, Q, n)
for ( i in 1:n ) {

  # generate prevalences
  theta[ , i] = c( 0.1, 0.45, 0.45 )
  theta[ , i] = as.vector(rdirichlet(1, c(1, 1, 1)))
  
  # simulate total read counts
  d_snp_all[i] <- round( rnorm(1, mean = mu_h*C_snp%*%theta[, i], sd = sigma_h ))
  d_snv_all[i] <- round( rnorm(1, mean = mu_h*C_snv%*%theta[, i], sd = sigma_h ))
  
  # compute variant allele frequency
  vaf_snp <- (M_snp%*%theta[, i])/(C_snp%*%theta[, i])
  vaf_snv <- (M_snv%*%theta[, i])/(C_snv%*%theta[, i])

  # simulate variant counts  
  k_snp_all[i] <- rbinom(1, d_snp_all[i], vaf_snp)
  k_snv_all[i] <- rbinom(1, d_snv_all[i], vaf_snv)
}


# preallocate matrices to store results
h_store <- matrix( 0, n, S ) # posterior means of h
lp_store <- matrix( 0, n, S ) # mean log posteriors
res_store <- matrix( 0, n, S ) # expected residual model fits to observations
theta_store <- list() # cellular prevalences

# for each configuration
for ( s in 1:S ) {
  
  # preallocate matrix to store cellular prevalences for this configuration
  theta_store[[s]] <- matrix( 0, Q, n)
  
  # for each observation
  for ( i in 1:n ) {
  
    # set data input for oncophase
    snpdat <- list(   Q = Q, 
                      mu_h = mu_h,
                      sigma_h = sigma_h,
                      e_a = e_a,
                      e_b = e_b,
                      alpha = alpha,
                      k_snv = k_snv_all[i],
                      d_snv = d_snv_all[i],
                      k_snp = k_snp_all[i], 
                      d_snp = d_snp_all[i],
                      M_snv = M_snv_table[s,],
                      C_snv = C_snv_table[s,],
                      M_snp = M_snp_table[s,],
                      C_snp = C_snp_table[s,] )
  
    # call oncophase
    f <- stan(file = 'oncophase.stan', data = snpdat, iter = 10000, chains = 1)
    
    # extract posterior summaries
    h_post_mean = mean(extract(f)$h)
    theta_post_mean = colMeans(extract(f)$theta)
    log_post_max = max(extract(f)$lp__)
    
    # store results
    h_store[i, s] = h_post_mean
    theta_store[[s]][1:Q, i] = theta_post_mean
    lp_store[i, s] = log_post_max
    
    # compute residual errors of model fit to observations
    d_snp_exp = h_post_mean*C_snp_table[s, ]%*%theta_post_mean
    d_snv_exp = h_post_mean*C_snv_table[s, ]%*%theta_post_mean
    vaf_snp_exp = (M_snp_table[s, ]%*%theta_post_mean)/(C_snp_table[s, ]%*%theta_post_mean)
    vaf_snv_exp = (M_snv_table[s, ]%*%theta_post_mean)/(C_snv_table[s, ]%*%theta_post_mean)
    
    d_res = abs( d_snp_all[i] - d_snp_exp )/d_snp_all[i] + abs( d_snv_all[i] - d_snv_exp )/d_snv_all[i]
    
    vaf_snp = k_snp_all[i]/d_snp_all[i]
    vaf_snv = k_snv_all[i]/d_snv_all[i]
    
    vaf_res = abs( vaf_snp - vaf_snp_exp )/vaf_snp + abs( vaf_snv - vaf_snv_exp )/vaf_snv 
    
    # store residual error
    res_store[i, s] = d_res + vaf_res
    
  }
}

boxplot(res_store)

