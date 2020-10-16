---
title: "Bayesian varying coefficient model with selection: an application to functional mapping"
author: "B. Heuclin, F. Mortier, C. Trottier, M. Denis"
date: "14/10/2020"
output:
   html_document:
    number_sections: yes
    toc: yes
---



We propose a Bayesian varying coefficient model  allowing, in a single step, the identification of variables involved in the variability of numerical response variable and the estimation of their dynamic effects. Selection is achieve using group spike-and-slab priors with either  P-spline interpolation  or  non-functional techniques to model the dynamic effects. This work is related to the paper: Bayesian varying coefficient model with selection: An application to functional mapping, B. Heuclin, F. Mortier, C. Trottier and M. Denis, Journal of the Royal Statistical Society: Series C (Applied Statistics), 2020,  DOI: 10.1111/rssc.12447



# Installation in R

``` r
library(devtools)
devtools::install_github("Heuclin/VCGSS")
```


# Exemple on simulated dataset

## Data simulation 
```{r}
n <- 200                      # number of individuals
q <- 50                       # number of variables
T <- 20                       # number of time points

se2 <- 4                      # residual variance is fixed to 7
rho <- 0.7                    # autoregressive decay parameter on the residual correlation structure is fixed to 0.7
nb_noeuds <- round(T/3)       # knots number for B-spline basis
 
# Simulation
sim <- list()

# intercept over time
sim$mu <-  1 + sin(pi*1:T/20) 

# plot the intercept onver time
plot(sim$mu, t='l'); abline(0, 0)


# Among the 50 variables, we fixed the first four variables with non zero effect over time
sim$beta <- matrix(0, T, q)
sim$beta[, 1] <- 4 + cumsum(rnorm(T, 0, 2))
sim$beta[, 2] <- 2*cos(pi * (1:T-25)/15)+1:T/50
sim$beta[, 3] <- c(rep(2, round(T/3)), rep(0, round(T/3)), rep(1, T-2*round(T/3)))
sim$beta[3, 3] <- 15
sim$beta[, 4] <- 2*60 / (25 + (1:T - T/2)^2)

# plot variable effects over time
par(mfrow = c(4,1), mar=c(2, 4,0.5, 1))
for(j in 1:4) plot(sim$beta[, j], col = 2)


# one environmental variable 
sim$env <-  cos(pi * (1:T-25)/2)+1:T/50  # sin(1:T/10)

# non linear effect of the environmental variable
sim$f_env <- sim$env*0.5 + sim$env^2 *0.3

# plot the additive part (intercept over time + environment effect)
par(mfrow= c(4, 1), mar = c(1, 4, 1, 1))
plot(sim$mu, t='l'); abline(0, 0)
plot(sim$env, t='l'); abline(0, 0)
plot(sim$f_env, t='l'); abline(0, 0)
plot(sim$mu + sim$f_env, t='l'); abline(0, 0)


sim$rho <-  rho
sim$se2 <- se2

# SNPs simulation
X <- scale(matrix(sample(0:2, n*q, replace = T), n, q))

# Construction of the autoregressive residual correlation structure
Gamma <- matrix(sim$rho, T, T); for(i in 1:T) for(j in 1:T) Gamma[i, j] <- sim$rho^abs(i-j) 

# Simulation of the response variable over time
Y <- matrix(NA, n, T)
for(i in 1:n){
  Y[i, ] <- sim$mu + sim$f_env + sim$beta %*% as.numeric(X[i, ]) + t(mvtnorm::rmvnorm(1, rep(0, T), sim$se2*Gamma))   
}

matplot(Y, t='l')

sim$f_env1 <- sim$f_env
envv <- matrix(sim$env, ncol = 1)
```


## Fit the varying coefficient model with spike-and-slab prior for variable selection with the main function `VCM_fct` 
```{r}
library(VCGSS)

fit <- VCM_fct(Y, X, ENV = envv, selection = TRUE, core = -1, rep = 6, niter = 10000, burnin = 5000)
```

## Check the MCMC convergence

### Gelman-Rubin's diagnostics 

By default the `VCM_fct` function provides Gelman Rubins diagnostic.

Gelman-Rubin's diagnostics may also be calculate using the following functions:

```{r}
# Diagnostic for alpha, m, e, pi, rho, se2
fit$gelman.diag

# gelman.plot
coda::gelman.plot(fit$mcmc_list)

# Diagnostic for beta
plot(fit$gelman.diag.b.psrf.median); abline(1.1, 0, col=2)  
# NA is due to variable which have not been selected

```

### Traceplot

```{r}
# traceplot
coda::plot(fit$mcmc_list)
```

### Visual convergence diagnostic of gamma (probabilities of variable inclusion)

```{r}
plot_diagnostic_gamma(fit)
```



## Parameter estimation
`fit$estimation` contains the estimations of different parameters 

```{r}
# residual variance
fit$estimation$se2 

# inclusion probability of variables
fit$estimation$pi   

# Autoregresive decay parameter
fit$estimation$rho  

# intercept
fit$estimation$alpha 

# posterior marginal probabilities of inclusion of variables
fit$estimation$mean.marginal.probabilities  

# curve over time of the intercept
fit$estimation$curve_mu    

# curve of the environment  
fit$estimation$curve_env   

# curve of the effect over time of each variables
fit$estimation$curve_beta  
```

### Visualisation of the marginal posterior probabilities of variable inclusion 

```{r}
prob <- fit$estimation$mean.marginal.probabilities
plot(prob, ylim = c(0, 1))
```

### Visualisation of the estimated dinamic effects over time

```{r}
plot_curve(fit, plot=c("Y", "mu", "beta"), mfrow = c(2, 3), add = c("matplot", "quantile"))
```

![](logo.png)

