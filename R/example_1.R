


# n <- 50                      # number of individuals
# q <- 20                      # number of variables
# T <- 20                       # number of time points
# 
# se2 <- 4                      # residual variance
# rho <- 0.7                    # autoregressive parameter on the residual correlation       structure
# nb_noeuds <- round(T/3)       # knot number for B-spline basis
# 
# #### Simulation
# sim <- list()
# 
# # intercept over time
# sim$mu <-  1 + sin(pi*1:T/20)
# 
# # plot the intercept onver time
# plot(sim$mu, t='l'); abline(0, 0)
# 
# 
# # among the 50 variables, we fixe the first four variables with non zero effect over time
# sim$beta <- matrix(0, T, q)
# sim$beta[, 1] <- 4 + cumsum(rnorm(T, 0, 2))
# sim$beta[, 2] <- 2*cos(pi * (1:T-25)/15)+1:T/50
# sim$beta[, 3] <- c(rep(2, round(T/3)), rep(0, round(T/3)), rep(1, T-2*round(T/3)))
# sim$beta[, 4] <- 2*60 / (25 + (1:T - T/2)^2)
# 
# # plot variable effects over time
# par(mfrow = c(4,1), mar=c(2, 4,0.5, 1))
# for(j in 1:4) {plot(sim$beta[, j], ylab = paste0("beta_", j)); lines(sim$beta[, j]); abline(0, 0, lty = 2)}
# 
# 
# # one environmental variable
# sim$env <-  cos(pi * (1:T-25)/2)+1:T/50  # sin(1:T/10)
# 
# # non linear effect of the environmental variable
# sim$f_env <- sim$env*0.5 + sim$env^2 *0.3
# 
# # plot the additive part (intercept over time + environment effect)
# par(mfrow= c(4, 1), mar = c(1, 4, 1, 1))
# plot(sim$mu); lines(sim$mu); abline(0, 0, lty = 2)
# plot(sim$env); lines(sim$env); abline(0, 0, lty = 2)
# plot(sim$f_env); lines(sim$f_env); abline(0, 0, lty = 2)
# plot(sim$mu + sim$f_env); lines(sim$mu + sim$f_env); abline(0, 0, lty = 2)
# 
# 
# sim$rho <-  rho
# sim$se2 <- se2
# 
# # SNPs simulation
# X <- scale(matrix(sample(0:2, n*q, replace = T), n, q))
# 
# # construction of the autoregressive residual correlation structure
# Gamma <- matrix(sim$rho, T, T); for(i in 1:T) for(j in 1:T) Gamma[i, j] <- sim$rho^abs(i-j)
# 
# # simulation of the response variable over time
# Y <- matrix(NA, n, T)
# for(i in 1:n){
#   Y[i, ] <- sim$mu + sim$f_env + sim$beta %*% as.numeric(X[i, ]) + t(mvtnorm::rmvnorm(1, rep(0, T), sim$se2*Gamma))
# }
# 
# # plot response variable over time
# par(mfrow= c(1, 1), mar = c(4, 4, 4, 1))
# matplot(t(Y), t='l', main = "matplot of Y", xlab = "time", ylab = "")
# 
# sim$f_env1 <- sim$f_env
# envv <- matrix(sim$env, ncol = 1)
# 
# 
# 
# 
# # Fit Bayesian varying coefficient with variable selection using group spike-and-slab prior
# fit <- VCM_fct(Y, X, ENV = envv, selection = TRUE, save = FALSE, cores = -1, rep = 4, niter = 10000, burnin = 5000)
# 
# 
# # Diagnostic for alpha, m, e, pi, rho, se2
# fit$gelman.diag
# 
# # gelman.plot
# coda::gelman.plot(fit$mcmc_list)
# 
# # Diagnostic for beta
# plot(fit$gelman.diag.b.psrf.median); abline(1.1, 0, col=2)
# # NA is due to variable which have not been selected
# 
# 
# # traceplot
# plot(fit$mcmc_list)
# 
# 
# #### Visual diagnostic of convergence ofmarginal posterior probabilities of variable
# # inclusion (gamma parameters)
# plot_diagnostic_gamma(fit)
# 
# 
# 
# ### Parameter estimation
# # intercept
# fit$estimation$alpha
# 
# # average inclusion probability of variables
# fit$estimation$pi
# 
# # posterior marginal probabilities of inclusion of variables
# fit$estimation$mean.marginal.probabilities
# 
# # curve over time of the intercept
# fit$estimation$curve_mu
# 
# # curve of the environment
# fit$estimation$curve_env
# 
# # curve of the effect over time of each variables
# fit$estimation$curve_beta
# 
# # residual variance
# fit$estimation$se2
# 
# # Autoregressive parameter
# fit$estimation$rho
# 
# 
# 
# 
# #### Visualisation of the marginal posterior probabilities of variable inclusion
# prob <- fit$estimation$mean.marginal.probabilities
# plot(prob, ylim = c(0, 1))
# 
# #### Visualisation of the estimated dynamic effects over time
# plot_functional_effects(fit, plot=c("Y", "mu", "env", "beta"), mfrow = c(3, 3), add = c("matplot", "quantile"))





