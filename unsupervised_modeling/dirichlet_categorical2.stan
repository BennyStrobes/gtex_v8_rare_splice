data {
  int<lower=0> N; // sample size
  int<lower=0> K; // number of classes
  vector[K] y[N]; // counts 
  real<lower=0> concShape; // concentration shape
  real<lower=0> concRate; // concentration rate
}
parameters {
  simplex[K] beta_raw; // symmetric K-1 dof encoding of beta
  real beta_scale; 
  real<lower=0> conc[K]; // concentration parameter
}

model {
  // beta reparameterization (Section 5.6 in the Stan manual)
  vector[K] beta;
  for (k in 1:K)
    beta[k] <- beta_scale * (beta_raw[k] - 1.0 / K);

  conc ~ gamma(concShape, concRate);
  for (n in 1:N) {
    vector[K] a; 
    real suma;
    vector[K] aPlusY;
    vector[K] lGaPlusY; 
    vector[K] lGaA ;
    vector[K] s; 
    s <- softmax(beta); 
    for (k in 1:K)
      a[k] <- conc[k] * s[k]; 
    // explicit construction of multinomial dirichlet
    // y ~ multinomial_dirichlet( conc * softmax(beta * x[n]) )
    suma <- sum(a);
    aPlusY <- a + y[n];
    for (k in 1:K) {
      lGaPlusY[k] <- lgamma(aPlusY[k]);
      lGaA[k] <- lgamma(a[k]);
    }
    increment_log_prob(lgamma(suma)+sum(lGaPlusY)-lgamma(suma+sum(y[n]))-sum(lGaA));
  }
}