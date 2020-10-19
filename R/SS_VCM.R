
#____________________________________________________________________________________________________________________
# 
#  Varying Coefficient Model (VCM) code
#  B. Heuclin, F. Mortier, C. Trottier, M. Denis
#  22/09/2020
# 
#____________________________________________________________________________________________________________________


library(mvnfast)
library(splines)
library(doParallel)
library(coda)
library(stringr)
library(rlang)
library(utils)



legendre_polynome <- function(n, t){
  res <- 0
  for(k in 0:n) res <- res + choose(n, k)^2 *(t-1)^(n-k) * (t+1)^k
  return(res* 1/2^n)
}


mode <- function(x) stats::density(x)$x[which.max(stats::density(x)$y)]


prob_jointe <- function(tab, start, end){
  proba.g.jointe <- sort(table(apply(tab[, start:end], 1,
                                     function(l) paste(l, collapse=""))) / 100, decreasing = TRUE )#/ nrow(g.est))
  for(i in 1:length(proba.g.jointe)){
    st <- names(proba.g.jointe)[i]
    if(sum(as.numeric(st)) == 0) names(proba.g.jointe)[i] <- "nothing"
    else names(proba.g.jointe)[i] <- paste0( paste(start - 1 + stringr::str_locate_all(st, "1")[[1]][, 1], collapse = "/"))
  }
  return(proba.g.jointe)
}







ss_vcm <- function(Y, X, settings, init, selection = TRUE, var_gp = "global")
{
  if(selection){
    ss_vcm_cpp(as.matrix(Y), as.matrix(X), as.list(settings), as.list(init), var_gp = var_gp)
  }else{
    vcm_cpp(as.matrix(Y), as.matrix(X), as.list(settings), as.list(init), var_gp = var_gp)
  }
}




#' Main function to run the Bayesian varying coefficient model using spike-and-slab prior for variable selection
#'
#' @param Y n x T matrix of phenotypic observations (n is the number of individuals, T is the number of observed time points)
#' @param X n x q matrix of variables (markers), q is the number of variables
#' @param ENV n x ne matrix of environmental variables, ne is the number of environmental variables
#' @param selection boolean indicating if selection is requiered, the default is TRUE
#' @param interpolation character, indicating the approach used to estimate curve of the effects. Possible choise is c('Legendre', 'B-spline', 'P-spline', 'RW'), the default is 'P-spline'.
#' @param df integer, only for Legendre, B-spline or P-spline approach. df is the degree of freedom parameter for these approaches, the default is floor(ncol(Y)/3)+3
#' @param order_diff interger in [1, 2, 3], order of the difference penalty, the default is 2.
#' @param niter number of iterations of the MCMC chain, the default is 10000
#' @param burnin number of initial iterations of the MCMC chain which should be discarded, the defaults to 5000
#' @param thin save only every thin-th iteration, the default is 10
#' @param rep number of repetitions, the default is one
#' @param cores number of cores for parallelisation, the default is -1 indicating that all cores minus 1 are used.
#' @param save boolean indicating if the MCMC chains should me saved, the default is TRUE
#' @param path optional (existing) path to save the MCMC chains if save = TRUE, if NULL, path = Fused_HS if selection = TRUE, Fusion_HS ortherwise
#' @param summary.fct name of the function (mean, mode or median) used to calculate the estimation of the parameters from the MCMC posterior samples, the default is mean.
#' @param scale boolean, indicating if X should be scaled, the default is TRUE
#' @param estimation boolean, indicating summary of MCMC chains should be apply to provide estimation of parameters and curve, the default is TRUE
#' @param gelman.diag boolean, indicating gelman.diag should be apply on MCMC chains, the default is TRUE
#' @param gelman.plot boolean indicating if the gelman plot should be plotted, the default is FALSE
#' @param traceplot boolean indicating if the trace plot should be plotted, the default is FALSE
#' @param ... additional arguments 
#'
#' @return Return a list
#' @export
#'
#' @examples # R/exemple_1.R
VCM_fct <- function(Y, X, ENV = NULL, selection = TRUE, interpolation='P-spline', df = floor(ncol(Y)/3)+3,  
                    order_diff = 2, niter = 10000, burnin = 5000, thin = 10,
                    rep = 1, cores = -1, save = TRUE, path = NULL, 
                    summary.fct = mean,  scale = FALSE, estimation = TRUE, 
                    gelman.diag = TRUE, gelman.plot=FALSE, traceplot=FALSE) #,  rm_old_chain = TRUE, ...)
{
  
  require(mvnfast)
  require(splines)
  require(doParallel)
  require(coda)
  require(stringr)
  require(rlang)
  require(utils)
  
  if(!interpolation %in% c('Legendre', 'B-spline', 'P-spline', 'RW')) stop("interpolation must be in the list : 'Legendre', 'B-spline', 'P-spline', 'RW'!")
  
  # print("--------------------------------")
  print("Bayesian Varying Coefficient Model using spike-and-slab prior")
  
  
  n <- nrow(Y); T <- ncol(Y); q <- ncol(X); n; q; T
  print(paste0("Number of individuals: ", n))
  print(paste0("Number of time points: ", T))
  print(paste0("Number of variables: ", q))
  
  if(cores >  parallel::detectCores()-1) warning(paste0("Number of cores can not be upper than ", parallel::detectCores()-1, " on your computer!"), immediate. = TRUE)
  if(cores == -1) cores <- parallel::detectCores()-1
  pars <- expand.grid(rep = 1:rep) #, fold = 1:CV)
  cores <- min(cores, nrow(pars))
  doParallel::registerDoParallel(cores = cores)
  print(paste("Nb cores = ", cores))
  
  
  if(is.null(path) & save){
    path <- "SS_VCM_results" 
    print(paste0("Default path: '", path, "'." ), immediate. = TRUE)
    system(paste0("mkdir ", path)) 
    
  }
  
  if(save){
    files <- system(paste0("ls ", path) , intern = TRUE)
    if(!rlang::is_empty(files)){
      warning(paste0(path, " directory is not empty !"), immediate. = TRUE)
      rm_old_chain <- utils::askYesNo("Do you want to delete all files ?" )
      if(is.na(rm_old_chain)) stop("Function have been stoped !")
      if(rm_old_chain) system(paste0("rm -rf ", path, "/* "))
    }
  }
  
  df_mu <- df_env <- df-1
  if(interpolation == "AR") df_beta <- T else df_beta <- df
  
  if(scale) X <- scale(X)
  output <- list()
  output$data$Y <- Y
  output$data$X <- X
  output$data$ENV <- ENV
  
  output$parameters$n <- n
  output$parameters$q <- q
  output$parameters$T <- T
  output$parameters$interpolation = interpolation
  output$parameters$order_diff = order_diff
  output$parameters$niter = niter
  output$parameters$burnin = burnin
  output$parameters$thin = thin
  output$parameters$cores = cores
  output$parameters$rep = rep
  output$parameters$save = save
  output$parameters$path = path
  output$parameters$gelman.plot = gelman.plot
  output$parameters$traceplot = traceplot
  output$parameters$summary.fct = summary.fct
  output$parameters$df_mu <- output$parameters$df_env <- df_mu
  output$parameters$df_beta <- df_beta
  
  # Settings _____________________________________________________________________________
  print("Settings")
  
  settings <- list()
  settings$l <- df           # p-spline : nombre de noeuds + degre spline  = 3
  settings$n.iter <- niter
  settings$burnin <- burnin
  settings$thin = thin
  settings$k <- 1               # nbr d'iterations pour le Metropolis-Hastings
  settings$Sigma_m0 <- 1e6*diag(settings$l)   # Matrice de la loi a priori sur m
  settings$a <- 1                # parametre de se2
  settings$bb <- 1               # parametre de se2
  settings$s <- 1                # parametre de lambda2
  settings$r <- 1                # parametre de lambda2
  settings$rho_tuning <- 0.05
  # settings$pi <- 0.1
  settings$s2_a <- 10
  tmp_temps <- 2*((1:T)-1)/(T-1) -1
  if(interpolation == 'Legendre'){
    settings$Bt <- matrix(1, T, settings$l); 
    for(t in 1:T) for(j in 1:settings$l) settings$Bt[t, j] <- legendre_polynome(j-1, tmp_temps[t])
  }else{
    settings$Bt <- as.matrix(splines::bs(tmp_temps, degree = 3, df = settings$l, intercept = TRUE)) # Base des b-spline
  }
  if(interpolation == 'RW') settings$B <- diag(T) 
  if(interpolation == 'B-spline' | interpolation == 'P-spline') settings$B <- as.matrix(splines::bs(tmp_temps, degree = 3, df = settings$l, intercept = TRUE)) # Base des b-spline
  if(interpolation == 'Legendre'){ settings$B <- matrix(1, T, settings$l); for(t in 1:T) for(j in 1:settings$l) settings$B[t, j] <- legendre_polynome(j-1, tmp_temps[t]) }
  
  if(!is.null(ENV)){
    settings$n_env <- ncol(ENV)
    settings$Be <- array(1, c(T, settings$l, settings$n_env))
    settings$Be_tmp <- array(NA, c(T, settings$l-1, settings$n_env))
    for(i in 1:settings$n_env){
      tmp_env <- (ENV[, i] - min(ENV[, i]))/(max(ENV[, i]) - min(ENV[, i])) *2 -1 
      if(interpolation == 'Legendre'){
        for(t in 1:T) for(j in 1:settings$l) settings$Be[t, j, i] <- legendre_polynome(j-1, tmp_env[t])
      }else{
        settings$Be[, , i] <- as.matrix(splines::bs(tmp_env, degree = 3, df = settings$l, intercept = TRUE)) # Base des b-spline
      }
    }
    Z=list()
    B_tilde <- settings$Be
    for ( j in (1:dim(B_tilde)[3])){
      QR = qr(t(t(rep(1, T))%*%B_tilde[ , , j]))
      QR_Q = qr.Q(QR,complete =TRUE)
      Z[[j]] = QR_Q[, 2:ncol(QR_Q)]
      settings$Be_tmp[, , j] <- settings$Be[, , j] %*% Z[[j]]
    }
    settings$Be <- settings$Be_tmp
  }else{
    settings$Be <- array(0, c(T, settings$l-1, 1))
    settings$n_env <- 1
  }
  
  if(interpolation == 'B-spline' | interpolation == 'Legendre') {
    settings$K <- diag(ncol(settings$B))
    settings$D <- diag(ncol(settings$B))
  }else{
    settings$D <- diff(diag(ncol(settings$B)), differences = order_diff)
    settings$K <- t(settings$D)%*%settings$D
    settings$K <- settings$K + 0.001 * diag(ncol(settings$B))
  }
  settings$K2 <- t(diff(diag(settings$l-1), differences = 2))%*%diff(diag(settings$l-1), differences = 2)
  settings$K2 <- settings$K2 + 0.001 * diag(settings$l-1)
  
  Z=list()
  B_tilde <- list(settings$Bt)
  QR = qr(t(t(rep(1, T))%*%B_tilde[[1]]))
  QR_Q = qr.Q(QR,complete =TRUE)
  Z[[1]] = QR_Q[, 2:ncol(QR_Q)]
  settings$Bt <- settings$Bt %*% Z[[1]]
  
  
  settings$epsilon = 1e-5
  
  output$settings <- settings
  
  # MCMC _____________________________________________________________________________
  print("MCMC sampler")
  
  
  
  list_chain <- foreach::foreach(k = 1:nrow(pars), .verbose = FALSE) %dopar% {
    init <- list()
    init$alpha <- stats::rnorm(1, 0, 3)
    init$pi = stats::runif(1, 0.001, 0.99)
    init$m <- stats::rnorm(settings$l-1, 0, 1)
    init$e <- matrix(0, settings$l-1, settings$n_env)
    init$b <- matrix(stats::rnorm( ncol(settings$B)*q, 0, 1), ncol(settings$B), q)
    init$rho <- stats::runif(1, 0.001, 0.99)            # rho, parametre auto-regressif sur la matrice de variance residuelle
    init$se2 <- abs(stats::rnorm(1, 1, 3))		           # sigma^2, variance residuelle, scalaire
    init$g <- rep(1, q) # sample(0:1, q, replace = TRUE);       # parametre gamma pour le Spike and Slab
    init$tau2 <- rep(100, q)    # tau2, parametres de groupe lasso, vecteur de longueur q
    init$tau0 <- 100
    init$tau0_e <- rep(100, settings$n_env)
    init$xi <- matrix(1, nrow(settings$D), q)
    
    # init <- list()
    # init$alpha <- 0
    # init$pi = 0.5
    # init$m <- rep(0, settings$l-1)
    # init$e <- matrix(0, settings$l-1, settings$n_env)
    # init$b <- matrix(0, ncol(settings$B), q)
    # init$tau2 <- rep(100, q)    # tau2, parametres de groupe lasso, vecteur de longueur q
    # init$rho <- 0.5            # rho, parametre auto-regressif sur la matrice de variance residuelle
    # init$se2 <- 5		           # sigma^2, variance residuelle, scalaire
    # init$g <- rep(1, q);       # parametre gamma pour le Spike and Slab
    # init$tau0 <- 100
    # init$tau0_e <- rep(100, settings$n_env)
    # init$xi <- matrix(1, nrow(settings$D), q)
    
    # chain <- ss_vcm(Y = Y, X = X, settings = settings, init = init, selection = selection)
    
    if(selection){
      chain <- ss_vcm_cpp(as.matrix(Y), as.matrix(X), as.list(settings), as.list(init), var_gp = "global")
    }else{
      chain <- vcm_cpp(as.matrix(Y), as.matrix(X), as.list(settings), as.list(init), var_gp = "global")
    }
    
    if(save) save(chain, settings, init, file = paste0(path, "/chain_rep_", k, ".Rdata"))
    
    if(save) return() else return(chain)
  }
  
  if(!save) {output$list_chain <- list_chain ; rm(list_chain)}
  if(save) save(output, file = paste0(path, "/output.Rdata"))
  
  # Estimation _____________________________________________________________________________
  if(estimation){
    output <- estimation.VCM(output)
    if(save) save(output, file = paste0(path, "/output.Rdata"))
  }
  
  # Gelman diag _____________________________________________________________________________
  if(rep>1 & gelman.diag) output <- gelman.diag.VCM(output, gelman.plot = gelman.plot, traceplot = traceplot) else output$gelman.diag <- NULL
  if(save) save(output, file = paste0(path, "/output.Rdata"))
  
  
  print("End")
  return(output)
  
}



#___________________________________________________________________________________________________________________







#' Gelman-Rubin's convergence diagnostics of the MCMC chain generated by the \code{VCM_fct} function
#'
#' This function can be apply only if the number of repetition applied in the \code{VCM_fct} function is upper than one.
#'
#' @param object list generated by the \code{VCM_fct} function 
#' @param gelman.plot boolean indicating if the gelman plot should be plotted 
#' @param traceplot boolean indicating if the trace plot should be plotted 
#'
#' @return return input \code{object} including Gelman Rubin's convergence diagnostics
#' @export
#'
#' @examples
gelman.diag.VCM <- function(object, gelman.plot = FALSE, traceplot = FALSE){
  require(coda)
  
  rep <- object$parameters$rep
  df_beta <- object$parameters$df_beta
  df_mu <- df_env <- object$parameters$df_mu
  q <- object$parameters$q
  ENV <- object$data$ENV
  path <- object$parameters$path
  save <- object$parameters$save
  
  print("Gelman diagnostic")
  mcmc_list_chain <- b <- m <- e <-  list()
  for(i in 1:q){
    b[[i]] <- list()
    for(j in 1:df_beta) b[[i]][[j]] <- list()
  }
  g <- NULL
  k=1
  for(k in 1:rep){
    if(save) load(paste0(path, "/chain_rep_", k, ".Rdata")) else chain <- object$list_chain[[k]]
    prob <- colMeans(chain$g)
    g <- cbind(g, prob)
    
    for(i in 1:q) for(j in 1:df_beta) if(prob[i]>0.5) b[[i]][[j]][[k]] <- coda::mcmc(chain$b[, j, i])
    
    tmp2 <-tmp3 <-  NULL
    tmp <- do.call(cbind, chain[c("alpha", "m")]); colnames(tmp) <- c("alpha", paste0("m", 1:(df_mu)))
    if(!is.null(ENV)) {for(i in 1:ncol(ENV)) tmp2 <- cbind(tmp2, chain$e[, , i]); colnames(tmp2) <-  paste0("e", 1:((df_env)*ncol(ENV)))}
    tmp3 <- do.call(cbind, chain[c("pi", "rho", "se2")]); colnames(tmp3) <- c("pi", "rho", "se2")
    mcmc_list_chain[[k]] <- coda::mcmc(cbind(tmp, tmp2, tmp3)) #, thin = thinin, start = burnin+1, end = niter)
  }
  object$mcmc_list <- coda::mcmc.list(mcmc_list_chain); rm(mcmc_list_chain)
  object$gelman.diag <- coda::gelman.diag(object$mcmc_list)
  print((object$gelman.diag))
  # graphics::plot(object$gelman.diag$psrf[, 1], ylab="psrf", xlab= "parameters")
  print("Summary of the Gelman diagnostic:")
  print(summary(object$gelman.diag$psrf))
  print("mpsrf:")
  print(object$gelman.diag$mpsrf)
  if(gelman.plot) coda::gelman.plot(object$mcmc_list)
  if(traceplot) plot(object$mcmc_list)
  
  # graphics::par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))
  # graphics::matplot(apply(g, 1, cumsum)/ 1:rep, t="l", lty = 1, lwd = 2, xlab = "number of repetitions", ylab = "average marg. prob.", ylim = c(0, 1))
  # # dim(g)
  prob <- rowMeans(g)
  
  
  gelman_diag_b <- matrix(NA, q, df_beta)
  for(i in 1:q){
    # if(prob[i] > 0.1)
    for(j in 1:df_beta){
      if(sum(sapply(b[[i]][[j]], length) > 0) > 1) gelman_diag_b[i, j] <- coda::gelman.diag(coda::mcmc.list(b[[i]][[j]][which(sapply(b[[i]][[j]], length) > 0)]))$psrf[1]
    }
  }
  object$gelman.diag.b.psrf <- round(gelman_diag_b, 2)
  object$gelman.diag.b.psrf.median <- round(apply(gelman_diag_b , 1, stats::median, na.rm = TRUE), 2)
  tmp <- object$gelman.diag.b.psrf.median
  if(sum(!is.na(tmp))) graphics::plot(object$gelman.diag.b.psrf.median, ylab = "median of psrf", xlab = "variables", main = "psrf of parameters 'b'"); graphics::abline(1.1, 0, col = 2)
  return(object)
}


#___________________________________________________________________________________________________________________





#' Function to obtain parameter estimations from many MCMC chains generated by the \code{VCM_fct} function 
#'
#' @param object list generated by the \code{VCM_fct} function 
#'
#' @return return the input \code{object} including the parameter estimations
#' @export
#'
#' @examples
estimation.VCM <- function(object){
  require(stringr)
  
  print("Summary")
  
  rep <- object$parameters$rep
  q <- object$parameters$q
  T <- object$parameters$T
  ENV <- object$data$ENV
  path <- object$parameters$path
  summary.fct <- object$parameters$summary.fct
  
  res_selec_mar <- res_selec_joint <- res_selec_prob <- matrix(0, rep, q)
  prob_joint <- rep(NA, rep)
  tmp <- matrix(NA, rep*q, T)
  curve_beta <- cbind(expand.grid(chain = 1:rep, pos = 1:q), tmp)
  if(!is.null(ENV)) {
    tmp <- matrix(NA, rep*ncol(ENV), T)
    curve_env <- cbind(expand.grid(chain = 1:rep, pos = 1:ncol(ENV)), tmp)
  }
  tmp <- matrix(NA, rep, T)
  curve_mu <- cbind(expand.grid(chain = 1:rep, pos = 1), tmp)
  
  full_chain <- list()
  full_chain$alpha <- NULL
  full_chain$m <- NULL
  full_chain$e <- NULL
  full_chain$b <- NULL
  full_chain$g <- NULL
  full_chain$rho <- NULL
  full_chain$se2 <- NULL
  g <- NULL
  i=1
  for(i in 1:rep){ # length(list_chain)){
    # print(i)
    if(!object$parameters$save ){
      chain <- object$list_chain[[i]]
    }else{
      load(paste0(path, "/chain_rep_", i, ".Rdata"))
    }
    g <- cbind(g, colMeans(chain$g))
    res_selec_mar[i, which(colMeans(chain$g[, ])>0.5)] <- 1
    proba.g.jointe <- sort(table(apply(chain$g[, 1:q], 1,
                                       function(l) paste(l, collapse=""))) / nrow(chain$g), decreasing = TRUE )[1]#/ nrow(g.est))
    prob_joint[i] <- proba.g.jointe
    res_selec_joint[i, stringr::str_locate_all(names(proba.g.jointe), "1")[[1]][, 1]] <- 1
    res_selec_prob[i, ] <-  colMeans(chain$g[, ])
    for(j in 1:q){
      if(res_selec_prob[i, j] >= 0.5){
        curve_beta[which(curve_beta$chain == i & curve_beta$pos == j), -(1:2)] <- object$settings$B %*% apply(chain$b[chain$g[, j]==1, , j], 2, summary.fct) #/ res_selec_prob[i, j]
      }
    }
    if(!is.null(ENV)) for(j in 1:ncol(ENV)) curve_env[which(curve_env$chain == i & curve_env$pos == j), -(1:2)] <- object$settings$Be[, , j] %*% apply(chain$e[, , j], 2, summary.fct) 
    curve_mu[which(curve_mu$chain == i), -(1:2)] <- object$settings$Bt %*% apply(chain$m[, ], 2, summary.fct) 
    
    full_chain$pi <- c(full_chain$pi, chain$pi)
    full_chain$alpha <- c(full_chain$alpha, chain$alpha)
    full_chain$m <- rbind(full_chain$m, chain$m)
    full_chain$e <- DescTools::Abind(full_chain$e, chain$e, along = 1)
    full_chain$b <- DescTools::Abind(full_chain$b, chain$b, along = 1)
    full_chain$g <- rbind(full_chain$g, chain$g)
    full_chain$rho <- c(full_chain$rho, chain$rho)
    full_chain$se2 <- c(full_chain$se2, chain$se2)
  }
  object$full_chain <- full_chain
  
  
  object$estimation$marginal.probabilities <- g
  object$estimation$mean.marginal.probabilities <-  colMeans(full_chain$g)
  object$estimation$res_selec_mar <- res_selec_mar
  
  object$estimation$joint_probabilities <- prob_joint
  object$estimation$res_selec_joint <- res_selec_joint
  object$estimation$res_selec_prob <- res_selec_prob
  
  object$estimation$curve_mu <- curve_mu
  if(!is.null(ENV)) object$estimation$curve_env <- curve_env else object$estimation$curve_env <- NULL
  object$estimation$curve_beta <- curve_beta
  object$estimation$pi <- summary.fct(full_chain$pi)
  object$estimation$rho <- summary.fct(full_chain$rho)
  object$estimation$alpha <- summary.fct(full_chain$alpha)
  object$estimation$se2 <- summary.fct(full_chain$se2)
  return(object)
}




#___________________________________________________________________________________________________________________






#' Function to plot the estimated functional effects
#'
#' @param object list generated by the \code{VCM_fct} function 
#' @param mfrow a numerical vector of the form c(nr, nc). Subsequent figures will be drawn in an nr-by-nc array on the device by rows. The default is c(6, 7).
#' @param mar a numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(1, 1, 4, 1).
#' @param plot caracter list ("Y", "mu", "env", "beta") containing object names for which you want the plot curves. The default is ("Y", "mu", "env", "beta").
#' @param id variable indinces for which you want to plot the functional effects. By defaults the argument id contains indicies of variables with posterior marginal probabilities upper than 0.5.
#' @param add caracter list ("quantile", "matplot") to indicate if you want to add credible intervals at 95\% or matpot of estimated curve of each repetition.
#' @param name optional caracter vector of variable names
#'
#' @return 
#' @export
#'
#' @example R/exemple_1.R
#' 
plot_functional_effects <- function(object, mfrow = c(6, 7), mar = c(1, 1, 4, 1), plot = c("Y", "mu", "env", "beta"), 
                       id = which(object$estimation$mean.marginal.probabilities > 0.5), add = c("quantile", "matplot"),
                       name = NULL)
{
  if(is.null(name)) name <- 1:object$parameters$q
  graphics::par(mfrow = mfrow, mar = mar)
  ENV <- object$data$ENV
  
  if("Y" %in% plot){
    Y <- object$data$Y
    graphics::matplot(t(Y), t='l', lty = 1, col = "gray65", main = "Y")
    graphics::lines(colMeans(Y), lwd = 2)
  }
  
  if("mu" %in% plot){
    if("matplot" %in% add){
      graphics::matplot(t(object$estimation$curve_mu[, -c(1:2)]), t='l', col = "gray65", lty = 1, ylab = "estimation", xlab = "Time", main = "mu"); graphics::abline(0, 0)
      graphics::lines(colMeans(object$estimation$curve_mu[, -c(1:2)]), lwd = 2, t='l')
    }else{
      graphics::plot(colMeans(object$estimation$curve_mu[, -c(1:2)]), lwd = 2, ylab = "estimation", xlab = "Time", main = "mu", t='l'); graphics::abline(0, 0)
    }
    if("quantile" %in% add) graphics::lines(apply(object$estimation$curve_mu[, -c(1:2)], 2, stats::quantile, 0.025), lty = 3, t="l", lwd = 2)
    if("quantile" %in% add) graphics::lines(apply(object$estimation$curve_mu[, -c(1:2)], 2, stats::quantile, 0.975), lty = 3, t="l", lwd = 2)
  }
  
  if("env" %in% plot){
    if(!is.null(ENV)){
      if(is.null(colnames(ENV))) colnames(ENV) <- 1:ncol(ENV)
      for(i in 1:ncol(ENV)) {
        if("matplot" %in% add){
          graphics::matplot(t(object$estimation$curve_env[which(object$estimation$curve_env$pos == i), -c(1:2)]), t='l', col = "gray65", lty = 1, ylab = "estimation", xlab = "Time", main = paste0("env ", colnames(ENV)[i])); graphics::abline(0, 0)
          graphics::lines(colMeans(object$estimation$curve_env[which(object$estimation$curve_env$pos == i), -c(1:2)]), lwd = 2, t='l')
        }else{
          graphics::plot(colMeans(object$estimation$curve_env[which(object$estimation$curve_env$pos == i), -c(1:2)]), lwd = 2, ylab = "estimation", xlab = "Time", main = paste0("env ", colnames(ENV)[i]), t='l'); graphics::abline(0, 0)
        }
        if("quantile" %in% add) graphics::lines(apply(object$estimation$curve_env[which(object$estimation$curve_env$pos == i), -c(1:2)], 2, stats::quantile, 0.025), lty = 3, t="l", lwd = 2)
        if("quantile" %in% add) graphics::lines(apply(object$estimation$curve_env[which(object$estimation$curve_env$pos == i), -c(1:2)], 2, stats::quantile, 0.975), lty = 3, t="l", lwd = 2)
      }
    }
  }
  
  if("beta" %in% plot){
    for(i in id) {
      if("matplot" %in% add){
        graphics::matplot(t(object$estimation$curve_beta[which(object$estimation$curve_beta$pos == i), -c(1:2)]), t='l', col = "gray65", lty = 1, ylab = "estimation", xlab = "Time", main = paste0("beta ", name[i])); graphics::abline(0, 0)
        graphics::lines(colMeans(object$estimation$curve_beta[which(object$estimation$curve_beta$pos == i), -c(1:2)], na.rm = TRUE), lwd = 2, t='l')
      }else{
        graphics::plot(colMeans(object$estimation$curve_beta[which(object$estimation$curve_beta$pos == i), -c(1:2)], na.rm = TRUE), lwd = 2, t='l', ylab = "estimation", xlab = "Time", main = paste0("beta ", name[i])); graphics::abline(0, 0)
      }
      if("quantile" %in% add) graphics::lines(apply(object$estimation$curve_beta[which(object$estimation$curve_beta$pos == i), -c(1:2)], 2, stats::quantile, 0.025, na.rm = TRUE), lty = 3, t="l", lwd = 2)
      if("quantile" %in% add) graphics::lines(apply(object$estimation$curve_beta[which(object$estimation$curve_beta$pos == i), -c(1:2)], 2, stats::quantile, 0.975, na.rm = TRUE), lty = 3, t="l", lwd = 2)
    }
  }
}



#___________________________________________________________________________________________________________________






#' Function to plot convergence diagnostic of variable marginal posterior inclusion probabilities
#'
#' @param object  list generated by the \code{VCM_fct} function 
#'
#' @return none
#' @export
#'
#' @example R/exemple_1.R
#' 
plot_diagnostic_gamma <- function(object){
  rep <- object$parameters$rep
  g <- object$estimation$marginal.probabilities
  graphics::par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
  graphics::matplot(apply(g, 1, cumsum)/ 1:rep, t="l", lty = 1, lwd = 2, xlab = "number of repetitions", ylab = "average marg. prob.", ylim = c(0, 1))
}






