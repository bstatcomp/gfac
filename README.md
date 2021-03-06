# gfac R package
The package contains several factor analysis methods, with
emphasis to count data. The main methods rely on copulas to estimate the
residual covariance, not explained by the latent structure. Additionally, the
package contains several data sets for empirical evaluation. The methods use RStan for inference.

# Installation
Download the contents of the package to a local folder. Then use

```{r eval = FALSE}
devtools::install("path", quick = TRUE)
``` 
where path is the path to the folder where you saved the contents of the package. The parameter `quick = TRUE` allows to install already compiled 64-bit models in the package. If you want to re-compile the models on your computer, you can install the package directly from git with

```{r eval = FALSE}
devtools::install_github("bstatcomp/gfac")
``` 

If you are installing this way, you may have trouble installing the 32-bit and 64-bit versions on Windows. In 
that case only install one version with

```{r eval = FALSE}
devtools::install_github("bstatcomp/gfac", INSTALL_opts = c("--no-multiarch"))
``` 

which will automatically install the arch that is used in your current RStudio
session. If it still does not install, you may need to check your system Path,
where the compiler for the desired version should be before the compiler for
the other version.

# Static factor analysis
The main functionality of the package is training and predicting with
different factor analysis methods. For that, you only need to use functions
starting with *train_* and *pred_*. Additionally, there are 20 data sets that
can be used for training, evaluation, and latent structure extraction.

For example, to train, predict, and evaluate
the statics negative binomial factor analysis with copula on EPL data, use

```{r eval = FALSE}
data(EPL_data)
train_ind <- EPL_data$splits[ ,1] == 2
test_ind  <- EPL_data$splits[ ,1] == 1
train_dat <- EPL_data$fa_data[train_ind, ]
train_grp <- EPL_data$group[train_ind]
test_dat  <- EPL_data$fa_data[test_ind, ]
test_grp  <- EPL_data$group[test_ind]

my_model <- train_CgNBFA(train_dat, train_grp)
my_eval  <- pred_CgNBFA(my_model, test_grp, test_dat)
``` 



# Factor analysis with smooth trajectories
These methods are an extension of static factor models. The latent factors are further modeled with smooth locally periodic Gaussian processes. Below is an example of fitting gNBGPFA and gNBGPFAC to a toy data set.
```{r eval = FALSE}
data(toy_A)
train_ind <- toy_A$splits[ ,1] == 0


# gNBGPFA
trained_mod <- train_gNBGPFA_LP(X = toy_A$df[train_ind, ], 
                                ts = toy_A$ts[train_ind], 
                                gp = toy_A$gp[train_ind], 
                                nfac = 3,
                                nchain = 1,
                                nit = 2000,
                                period_length = 2,
                                prior_r2 = c(4.6, 22),
                                prior_r3 = c(4.6, 22),
                                prior_phi_type = 1,
                                prior_phi = c(1, 0.01),
                                transform_t = function (t) return (t))
predicted_mod <- pred_gNBGPFA_LP(mod = trained_mod, 
                                 X_test = toy_A$df[!train_ind, ], 
                                 ts_test = toy_A$ts[!train_ind], 
                                 gp_test = toy_A$gp[!train_ind], 
                                 ts_train = toy_A$ts[train_ind], 
                                 gp_train = toy_A$gp[train_ind],
                                 transform_t = function (t) return (t),
                                 period_length = 2)


# gNBGPFAC
# Set initial values for gNBGPFAC
ext <- rstan::extract(trained_mod$samps)
init_list <- list( # initial values from gNBGPFA
  simp = apply(ext$simp, c(2,3), mean),
  fn_pri = apply(ext$fn_pri, c(2,3), mean),
  base = apply(ext$base, c(2,3), mean),
  phi = apply(ext$phi, 2, mean),
  a2 = apply(ext$a2, c(2,3), mean),
  r2 = apply(ext$r2, c(2,3), mean),
  r3 = apply(ext$r3, c(2,3), mean),
  L  = diag(1, nrow = 8)
)

trained_mod2 <- train_gNBGPFAC_LP(X = toy_A$df[train_ind, ], 
                                  ts = toy_A$ts[train_ind], 
                                  gp = toy_A$gp[train_ind], 
                                  nfac = 3,
                                  nchain = 1,
                                  nit = 2000,
                                  period_length = 2,
                                  prior_r2 = c(4.6, 22),
                                  prior_r3 = c(4.6, 22),
                                  prior_phi_type = 1,
                                  prior_phi = c(1, 0.01),
                                  transform_t = function (t) return (t),
                                  init_list = init_list)
predicted_mod2 <- predict_gNBGPFAC_LP(mod = trained_mod2, 
                                      X_test = toy_A$df[!train_ind, ], 
                                      ts_test = toy_A$ts[!train_ind], 
                                      gp_test = toy_A$gp[!train_ind], 
                                      ts_train = toy_A$ts[train_ind], 
                                      gp_train = toy_A$gp[train_ind],
                                      transform_t = function (t) return (t),
                                      period_length = 2)

``` 
