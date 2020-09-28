# Poisson z low count-i

my_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(my_path)
setwd("..")


# Functions
se_kernel <- function(x, y, alpha1 = 1, rho1) {
  alpha1 * exp(- (x - y)^2 / (2 * rho1^2))
}

lp_kernel <- function(x, y, alpha2, rho2, rho3) {
  alpha2 * exp((-2 * (sin(pi * abs(x - y) / 2))^2) / rho2^2) *
    exp((-(x - y)^2) / (2 * rho3))
}

p_kernel <- function(x, y, alpha2, rho2) {
  alpha2 * exp((-2 * (sin(pi * abs(x - y) / 2))^2) / rho2^2)
}

draw_samples <- function(x, N, seed = 1, kernel_fn, ...) {
  Y <- matrix(NA, nrow = length(x), ncol = N)
  set.seed(seed)
  for (n in 1:N) {
    K <- cov_matrix(x, kernel_fn, ...)
    Y[, n] <- mvrnorm(1, mu = rep(0, times = length(x)), Sigma = K)
  }
  Y
}

cov_matrix <- function(x, kernel_fn, ...) {
  outer(x, x, function(a, b) kernel_fn(a, b, ...))
}


library(MASS)
library(MCMCpack)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)

# 2 latent dimensions, governed by a Gaussian process --------------------------
x  <- seq(0.1, 20, by = 0.1)
F1 <- draw_samples(x, 1, seed = 2, se_kernel, alpha1 = 0.5, rho1 = 4 * 2)
F2 <- draw_samples(x, 1, seed = 10, lp_kernel, alpha2 = 0.5, rho2 = 2 * 2, rho3 = 7 * 2)
F3 <- draw_samples(x, 1, seed = 10, p_kernel, alpha2 = 20, rho2 = 6 * 2)
F4 <- draw_samples(x, 1, seed = 42, se_kernel, alpha1 = 0.5, rho1 = 3 * 2)
F5 <- draw_samples(x, 1, seed = 21, lp_kernel, alpha2 = 0.5, rho2 = 1 * 2, rho3 = 8 * 2)
F6 <- draw_samples(x, 1, seed = 10, p_kernel, alpha2 = 20, rho2 = 5 * 2)

plot(F1, type = "l", ylim = c(-2,2))
plot(F2, type = "l", ylim = c(-2,2))
plot(F3, type = "l", ylim = c(-2,2))
plot(F4, type = "l", ylim = c(-2,2))
plot(F5, type = "l", ylim = c(-2,2))
plot(F6, type = "l", ylim = c(-2,2))


F_mat <- as.matrix(data.frame(D1 = c(F1, F4), D2 = c(F2, F5), D3 = c(F3, F6)))

# we also need base measure


# Generate X from the loading matrix -------------------------------------------
# loadings 10 x 3
L_mat <- matrix(c(1, 0, 0,
                  0.1, 0.1, -0.3,
                  -1, 1, 0,
                  1.5, 0, 0.5,
                  -2, -1, 1,
                  1.2, 1, -0.4,
                  0, -1.5, 1,
                  1, 1, 0.3,
                  0, 0.6, 0),
                ncol = 3,
                byrow = T)
X               <- t(L_mat %*% t(F_mat))
minX            <- apply(X, 2, min)
minX            <- rnorm(9, 50, 5)
# minX[minX >= 0] <- 0
X               <- sweep(X * 4, 2, minX, "+")


X_tmp <- as.data.frame(X)
colnames(X_tmp) <- paste0("X", 1:ncol(X_tmp))
x_out_df        <- cbind(X_tmp, x = x, gp = rep(c(1,2), each = length(x)))
x_out_long      <- gather(x_out_df, key = "key", value = "value", - c(x, gp))
ggplot(x_out_long, aes(x = x, y = value, color = as.factor(gp))) + geom_line() + facet_wrap(~ key)
x1_ol <- cbind(x_out_long, type = "mean")


# Poisson likelihood -----------------------------------------------------------
set.seed(1)
cor_mat <- diag(1, nrow = 9)
cor_mat[1,2]  <- 0.4
cor_mat[2,1]  <- 0.4
cor_mat[1,4]  <- -0.3
cor_mat[4,1]  <- -0.3
cor_mat[2,9] <- 0.6
cor_mat[9,2] <- 0.6
cor_mat[3,6]  <- -0.8
cor_mat[6,3]  <- -0.8
cor_mat[6,7]  <- 0.7
cor_mat[7,6]  <- 0.7


## random covariance
set.seed(1)
A <- matrix(rnorm(81), nrow = 9)
B <- A %*% t(A)
cor_mat <- cov2cor(B)
cor_mat[abs(cor_mat) < 0.3] <- 0
cor_mat <- Matrix::nearPD(cor_mat)$mat

gaussianC_p <- function(cor_mat, mu){
  M  <- length(mu)
  vX <- mvrnorm(1, rep(0, M), cor_mat)
  vU <- pnorm(vX)
  lY <- vector(length = M)
  lY <- qpois(vU, mu)
  return(lY)
}
X_out <- X
for (i in 1:nrow(X_out)) {
  X_out[i, ] <- gaussianC_p(cor_mat, X[i, ])
}


X_tmp <- as.data.frame(X_out)
colnames(X_tmp) <- paste0("X", 1:ncol(X_tmp))
x_out_df        <- cbind(X_tmp, x = x, gp = rep(c(1,2), each = length(x)))
x_out_long      <- gather(x_out_df, key = "key", value = "value", - c(x, gp))
ggplot(x_out_long, aes(x = x, y = value, color = as.factor(gp))) + geom_line() + facet_wrap(~ key) + geom_point()
x2_ol <- cbind(x_out_long, type = "Poisson")




## Poisson vs mean
dft <- rbind(x1_ol, x2_ol)
ggplot(dft, aes(x = x, y = value, color = as.factor(type))) + geom_line() + facet_wrap(~ key + as.factor(gp))
# ggsave("./data_clean/plots/toy_D.pdf", width = 16, height = 14)


## Split -----------------------------------------------------------------------
n      <- nrow(X_tmp)
split1 <- rep(0, n)
split2 <- rep(0, n)
samps  <- sample(1:n, size = n / 2)
split1[samps] <- 1
split2[split1 == 0] <- 1


splits <- data.frame(split1, split2)


# Save data --------------------------------------------------------------------
data_out <- list(
  df     = X_tmp,
  ts     = c(x, x),
  gp     = rep(c("T1", "T2"), each = length(x)),
  splits = splits
)
toy_D <- data_out

save(toy_D, file = "./data/toy_D.RData")

