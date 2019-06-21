data {
  int<lower=0> K; // number of classes
  vector[K] y; // counts 
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
  vector[K] alpha; 
  real suma;
  vector[K] aPlusY;
  vector[K] lGaPlusY; 
  vector[K] lGaA ;
  vector[K] s; 
  for (k in 1:K)
    beta[k] <- beta_scale * (beta_raw[k] - 1.0 / K);
  conc ~ gamma(concShape, concRate);
  s <- softmax(beta); 
  for (k in 1:K)
    alpha[k] <- conc[k] * s[k]; 
  // explicit construction of categorical dirichlet
  suma <- sum(alpha);
  aPlusY <- alpha + y;
  for (k in 1:K) {
    lGaPlusY[k] <- lgamma(aPlusY[k]);
    lGaA[k] <- lgamma(alpha[k]);
  }
  increment_log_prob(lgamma(suma)+sum(lGaPlusY)-lgamma(suma+sum(y))-sum(lGaA));
}