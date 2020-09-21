functions {
  vector ghk_cop(row_vector u, matrix L, int[] X, vector mu, vector phi){
    int M = cols(u);
    matrix[M,2] ab;
    vector[M] a;
    vector[M] b;
    vector[M] a_tilde;
    vector[M] b_tilde;
    real tmp_a;
    real tmp_b;
    vector[M] tmp;
    vector[M] lhood;
    for(m in 1:M){
      int m1 = m - 1;
      if(X[m]<0.001){
        a[m] = -1;
        b[m] = inv_Phi(neg_binomial_2_cdf(X[m], mu[m], phi[m]));
        tmp_b = (b[m] - ((m > 1) ? L[m,1:m1] * head(tmp,m1) : 0)) / L[m,m];
        a_tilde[m] = 0;
        b_tilde[m] = Phi(tmp_b);
      }else{
        a[m] = inv_Phi(neg_binomial_2_cdf(X[m]-1, mu[m], phi[m]));
        b[m] = inv_Phi(neg_binomial_2_cdf(X[m], mu[m], phi[m]));
        tmp_a = (a[m] - ((m > 1) ? L[m,1:m1] * head(tmp,m1) : 0)) / L[m,m];
        tmp_b = (b[m] - ((m > 1) ? L[m,1:m1] * head(tmp,m1) : 0)) / L[m,m];
        a_tilde[m] = Phi(tmp_a);
        b_tilde[m] = Phi(tmp_b);
      }
      tmp[m] = inv_Phi(a_tilde[m] + (b_tilde[m] - a_tilde[m]) * u[m]);
      lhood[m] = b_tilde[m] - a_tilde[m];
    }
    return(lhood);
  }
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
      // locally periodic
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
  int X[M, N];
  
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

  vector<lower=0>[M] phi;
  
  // copula parameters
  cholesky_factor_corr[M] L;
  matrix<lower=0,upper=1>[N,M] u;
}
transformed parameters {
  matrix[M, P] Psi;
  matrix[M, N] mu_fac;
  matrix[P,N] fn;
  
  // copula parameters
  matrix[N,M] lhood;
  
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
  mu_fac = exp(mu_fac);
  
  // copula
  for(i in 1:N){
    lhood[i] = ghk_cop(u[i], L, X[ ,i], mu_fac[ ,i], phi)';
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
  L ~ lkj_corr_cholesky(2.0);
  for (m in 1:M) {
    simp[m] ~ normal(0, 1);
    if (p_dist[3] == 1) {
      phi[m]  ~ gamma(p_val[5], p_val[6]);
    } else if (p_dist[3] == 2) {
      phi[m] ~ inv_gamma(p_val[5], p_val[6]);
    }
    for (n in 1:N) {
      target += log(lhood[n,m]);
    }
  }
}


