functions {
  matrix GP_LPK(
    int Nx,
    real[] x,
    int Ny,
    real[] y,
    real alpha2,
    real rho2,
    real rho3,
    real period_length){
      matrix[Nx, Ny] K2;
      matrix[Nx, Ny] Sigma;
      for(i in 1:Nx){
        for(j in 1:Ny){
          K2[i, j] = alpha2*exp(-2*square(sin(pi()*fabs(x[i]-y[j]) / period_length))/square(rho2))*
                     exp(-square(x[i]-y[j])/2/square(rho3));
        }
      }
    Sigma = K2;
    return Sigma;
  }
}
data {
  int<lower=1> N;
  int<lower=1> NG;
  int<lower=1> M;
  int<lower=1> P;
  real t[N];
  int<lower=1> G[NG];
  matrix[M, N] X;

  // priors
  int<lower=0> p_dist[3];
  real p_val[6];

  // period length
  real<lower=0> period_length;
}
transformed data {
  real sq_sigma = 0.00000001;
}
parameters {
  matrix[M,P] simp;
  matrix[P,N] fn_pri;
  vector[M] base[NG];

  vector<lower=0>[P] a2[NG];
  vector<lower=0>[P] r2[NG];
  vector<lower=0>[P] r3[NG];

  vector<lower=0>[M] sigma;
}
transformed parameters {
  matrix[M, P] Psi;
  matrix[M, N] mu_fac;
  matrix[P,N] fn;

  for (p in 1:P) {
    int pos;
    pos = 1;
    for (ng in 1:NG) {
      matrix[G[ng], G[ng]] L_K_tmp;
      matrix[G[ng], G[ng]] K;
      K = GP_LPK(G[ng], t[pos:(pos + G[ng] - 1)], G[ng], t[pos:(pos + G[ng] - 1)], a2[ng,p], r2[ng,p], r3[ng,p], period_length);
      K = K + diag_matrix(rep_vector(sq_sigma, G[ng]));
      L_K_tmp = cholesky_decompose(K);
      fn[p,pos:(pos + G[ng] - 1)]    = fn_pri[p, pos:(pos + G[ng] - 1)] * L_K_tmp';
      pos = pos + G[ng];
    }
    for (m in 1:M) {
      if (p == m) {
        Psi[m,p] = 1;
      } else if (p > m) {
        Psi[m,p] = 0;
      } else {
        Psi[m,p] = simp[m,p];
      }
    }
  }
  mu_fac = Psi * fn;
  {
    int pos;
    pos = 1;
    for (ng in 1:NG) {
      mu_fac[ ,pos:(pos + G[ng] - 1)] = mu_fac[ ,pos:(pos + G[ng] - 1)] + base[ng, ] * rep_vector(1, G[ng])';
      pos = pos + G[ng];
    }
  }
}
model {
  for (ng in 1:NG) {
    a2[ng, ]  ~ student_t(3,0,1);
    if (p_dist[1] == 1) {
      r2[ng, ]  ~ inv_gamma(p_val[1], p_val[2]);
    }
    if (p_dist[2] == 1) {
      r3[ng, ]  ~ inv_gamma(p_val[3], p_val[4]);
    }
    base[ng, ] ~ normal(0, 10);
  }
  for (p in 1:P) {
    fn_pri[p] ~ normal(0, 1);
  }
  for (m in 1:M) {
    simp[m]  ~ normal(0, 1);
    sigma[m] ~ student_t(3,0,1);
    X[m]     ~ normal(mu_fac[m], sigma[m]);
  }
}


