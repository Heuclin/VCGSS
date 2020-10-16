
# Bayesian varying coefficient model with selection: an application to functional mapping

Benjamin Heuclin, Frédéric Mortier, Catherine Trottier, Marie Denis


## Exemple on simulated dataset

### Data simulation 
```{r}
n <- 200
q <- 50
T <- 20

se2 <- 4
rho <- 0.7
nb_noeuds <- round(T/3)
 
#_____ valeurs simulées
sim <- list()


sim$mu <-  1 + sin(pi*1:T/20) 
sim$beta <- matrix(0, T, q)
sim$beta[, 1] <- 4 + cumsum(rnorm(T, 0, 2))
sim$beta[, 2] <- 2*cos(pi * (1:T-25)/15)+1:T/50
sim$beta[, 3] <- c(rep(2, round(T/3)), rep(0, round(T/3)), rep(1, T-2*round(T/3)))
sim$beta[3, 3] <- 15
sim$beta[, 4] <- 2*60 / (25 + (1:T - T/2)^2)

sim$env <-  cos(pi * (1:T-25)/2)+1:T/50  # sin(1:T/10)
sim$f_env <- sim$env*0.5 + sim$env^2 *0.3

par(mfrow= c(4, 1), mar = c(1, 4, 1, 1))
plot(sim$mu, t='l'); abline(0, 0)
plot(sim$f_env, t='l'); abline(0, 0)
plot(sim$mu + sim$f_env, t='l'); abline(0, 0)

sim$rho <-  rho
sim$se2 <- se2

# par(mfrow = c(3,1), mar=c(2, 4,0.5, 1))
# plot(sim$env, t='l')
# plot(sim$env, sim$f_env, t='l')
# plot(1:T, sim$f_env2, t='l')
# 
# par(mfrow = c(5,1), mar=c(2, 4,0.5, 1))
# plot(sim$mu, col =2, ylab = "mu")
# for(j in 1:4) plot(sim$beta[, j], col = 2)

# X <- scale(mvtnorm::rmvnorm(n, rep(0, q), 1*Gamma));
# View(cor(X))
# X <- scale(matrix(runif(n*q, 0,1), n, q))
X <- scale(matrix(sample(0:2, n*q, replace = T), n, q))
Gamma <- matrix(sim$rho, T, T); for(i in 1:T) for(j in 1:T) Gamma[i, j] <- sim$rho^abs(i-j) # matrice intervenant dans la variance résiduelle

Y <- matrix(NA, n, T)
for(i in 1:n){
  Y[i, ] <- sim$mu + sim$f_env + sim$beta %*% as.numeric(X[i, ]) + t(mvtnorm::rmvnorm(1, rep(0, T), sim$se2*Gamma))   
}

sim$f_env1 <- sim$f_env
envv <- matrix(sim$env, ncol = 1)
```

### Source the R and Cpp functions 
```{r}
source('SS_VCM.R')
```
### Fit the varying coefficient model with spike-and-slab prior for markers selection 
```{r}
fit <- VCM_fct(Y, X, ENV = envv, selection = TRUE, core = -1, rep = 6, niter = 10000, burnin = 5000)
```

By default the `VCM_fct` function applies Gelman Rubins diagnostic.
The Gelman Rubins diagnostic may also be calculate using the following functions:

```{r}
# Diagnostic fot alpha, m, e, pi, rho, se2
fit$gelman.diag

# gelman.plot
coda::gelman.plot(fit$mcmc_list)

# traceplot
plot(fit$mcmc_list)

# Diagnostic for beta
plot(fit$gelman.diag.b.psrf.median); abline(1.1, 0, col=2)  # NA is due to variable which have not been selected

# Diagnostic for gamma (probabilities of inclusion of beta)
plot_diagnostic_gamma(fit)

```

`fit$estimation` contains the estimations of different parameters 

```{r}
fit$estimation$se2    # residual variance
fit$estimation$pi     # inclusion probability of variables
fit$estimation$rho    # Autoregresive decay parameter
fit$estimation$alpha  # intercept
# fit$estimation$mean.marginal.probabilities  # posterior marginal probabilities of inclusion of variables
# fit$estimation$curve_mu    # curve over time of the intercept
# fit$estimation$curve_env   # curve of the environment   (not available for this application (ENV = NULL))
# fit$estimation$curve_beta  # curve of the effect over time of each variables
```

We can plot the posterior marginal probabilities of inclusion of variables

```{r}
prob <- fit$estimation$mean.marginal.probabilities
plot(prob, ylim = c(0, 1))
```


```{r}
plot_curve(fit, plot=c("Y", "mu", "beta"), mfrow = c(2, 3), add = c("matplot", "quantile"))

```



