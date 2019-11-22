data {
  int<lower=0> n;
  int<lower=0> m;
  int<lower=0> p;
  int g[n];
  int<lower=1> ng;
  int X[m,n];
}

parameters {
  matrix<lower=0>[p,ng] Theta;
  vector<lower=0>[p] rk;
  vector<lower=0, upper=1>[p] pk;
  simplex[m] simp[p];
}

transformed parameters {
  matrix<lower=0>[m,p] Psi;
  for (i in 1:p) {
    Psi[ ,i] = simp[i];
  }
}

model {
  for (i in 1:p) {
    simp[i] ~ dirichlet(rep_vector(1.0, m));
  }
  for (i in 1:ng) {
    for (j in 1:p) {
      Theta[j,i] ~ gamma(rk[j], pk[j] / (1 - pk[j]));
    }
  }
  for (i in 1:n)
    X[ ,i] ~ poisson(Psi * Theta[ ,g[i]]);
}
