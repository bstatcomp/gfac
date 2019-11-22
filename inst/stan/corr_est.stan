data {
  int n;
  int m;
  matrix[n,m] X;
}

parameters {
  corr_matrix[m] Sigma;
}

model {
  Sigma ~ lkj_corr(2);
  for (i in 1:n)
    X[i, ] ~ multi_normal(rep_vector(0, m), Sigma);
}

