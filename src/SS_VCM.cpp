#include <RcppArmadillo.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerVector my_sample(IntegerVector A, int n){
  IntegerVector r (n);
  r = RcppArmadillo::sample(A, n, FALSE, NumericVector::create());
  return(r);
}

//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
double rInvGauss( double mu, double l){
  double v = R::rnorm(0, 1);
  double y = v*v;
  double x = mu + mu * mu * y /(2*l) - mu * sqrt(4*mu * l * y + mu * mu * y * y) / (2*l);
  double z = R::runif(0, 1);
  if(z<= mu/(mu+x)) return(x); else return(mu * mu / x);
  return(x);
}

//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
double mmax(double a, double b){
  if(a>b) 
    return(a); 
  else 
    return(b);
}

//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
double mmin(double a, double b){
  if(a<b) 
    return(a); 
  else 
    return(b);
}

//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec rmnorm(arma::vec mean, arma::mat Sigma){
  int n = mean.n_elem;
  arma::mat s(n, n); //, d(n, n);
  // arma::vec v(n);
  // arma::svd(s, v, d, Sigma);
  // arma::vec r = s* diagmat(sqrt(v)) *arma::randn(n, 1) + mean;
  arma::chol(s, Sigma);
  arma::vec r = s.t() *arma::randn(n, 1) + mean;
  return(r);
}

//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec rmnorm_svd(arma::vec r, arma::vec mean, arma::mat V, arma::vec s){
  r = V* diagmat(sqrt(s)) *(arma::randn(mean.n_elem, 1) + mean);
  return(r);
}

//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec rmnorm_svd2(arma::vec mean, arma::mat V, arma::vec s){
  arma::vec r = V* diagmat(sqrt(s)) *(arma::randn(mean.n_elem, 1) + mean);
  return(r);
}






//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
List gl_vcm_cpp(arma::mat Y, arma::mat X, List settings, List init) {
  //
  //
  //
  //  Input :
  //    Y : matrice (n x T), observation dynamique
  //    X : matrice (n x q), marqueurs moleculaire
  //    init : list, valeurs initiales pour les parametres
  //    settings : list, reglages
  //    sim : liste, parametres simules pour les graphs
  //
  //  Output :
  //    Chain : liste contenant six matrices, chaques matrices contient les echantillonges selont les
  //            lois conditionnelles completes de chaques parametres (m, b, to, rho, se2, lambda2)
  
  int n = Y.n_rows, T = Y.n_cols, q = X.n_cols, j;
  
  // Rcout << " 0,1 \n";
  // parametre :
  int nIter = settings["n.iter"];             // nombre d'iteration de gibbs
  int k = settings["k"];                      // Nombre d'iteration de M-H pour gamma.head(dg)
  int l = settings["l"];                      // degre des spline + nombre de noeuds
  double a = settings["a"];                   // parametre de sigma^2
  double bb = settings["bb"];                 // parametre de sigma^2
  double s = settings["s"];                   // parametre shape de tau2
  double r = settings["r"];                   // parametre rate de tau2
  double rho_tuning = settings["rho_tuning"]; // Parametre pour le MH de rho
  arma::mat B = settings["B"];                // Matrice de base temp
  arma::mat Bt = settings["Bt"];              // Matrice de base temps (QR decomposition)
  arma::mat Be = settings["Be"];              // Matrice de base environnement (QR decomposition)
  arma::mat K = settings["K"];                // matrice des differences (p-spline)
  // double s2_a = settings["s2_a"];              // variance de alpha
  
  // Rcout << " 0,2 \n";
  // Initialisation :
  arma::vec m = init["m"];			              //  m, intercept, vecteur de longueur v (degre de l'approximation par polynome de legendre
  arma::vec e = init["e"];			              //  e, environment, vecteur de longueur v (degre de l'approximation par polynome de legendre
  arma::mat b = init["b"];			              // b, coefficient dynamique, matrice de taile (T x q)
  double rho = init["rho"];                   // rho, parametre auto-regressif sur la matrice de variance residuelle
  double se2 = init["se2"];	                  // sigma^2, variance residuelle, scalaire
  arma::vec tau2 = init["tau2"];              // variance spike and slab sur beta
  double lambda2 = init["lambda2"]; // lambda^2, parametre de penalisation
  double tau0 = init["tau0"];                 // variance de la loi a priori sur m
  double tau0_e = init["tau0_e"];             // variance de la loi a priori sur e
  double alpha = init["alpha"];
  
  // Rcout << " 0,3 \n";
  // declaration des matrices de sortie
  arma::mat chain_m (nIter/2+1, l); //colnames(chain.m) = paste0("m.", 1:v);
  arma::mat chain_e (nIter/2+1, l);
  arma::cube chain_b (nIter/2+1, l, q);
  arma::mat chain_lambda1 (nIter/2+1, q);
  arma::mat chain_tau2 (nIter/2+1,q);
  arma::vec chain_rho (nIter/2+1);
  arma::vec chain_se2 (nIter/2+1);
  arma::vec chain_tau0 (nIter/2+1);
  arma::vec chain_tau0_e (nIter/2+1);
  arma::vec chain_alpha (nIter/2+1);
  arma::vec chain_lambda2 (nIter/2+1);
  
  chain_m.row(0) = m.t();
  chain_e.row(0) = e.t();
  for(j=0; j<=q-1; j++){
    chain_b.slice(j).row(0) = b.col(j).t();
  }
  chain_tau2.row(0) = tau2.t();
  chain_rho(0) = rho;
  chain_se2(0) = se2;
  chain_tau0(0) = tau0;
  chain_tau0_e(0) = tau0_e;
  chain_alpha(0) = alpha;
  chain_lambda2(0) = lambda2;
  
  // Declaration des variables internes :
  IntegerVector ind = Range(0, q-1);
  int iter = 1, i, ii, perc=2, iter_tmp = 1;
  double rho_i, rho_star, tmp_i, tmp_star, log_Qnew, log_Qold, rate_se2, log_ratio, a0, a1, a2, var_alpha, mean_alpha;
  arma::vec tmp_m(T), tmp_a(T), mean_m(l), tmp_e(T), mean_e(l), mean_bj(l), tmp_mean_bj(T), tmp(T), tmp_mean_cj(T), diag_0_Gamma(T), diag_1_Gamma(T);
  arma::vec tmp2(T), s_tmp(l), X2 (q), sumX(q), sumY(T), b2(q), one = arma::ones(T);
  arma::mat Gamma(T, T),Gamma_i(T, T), Gamma_star(T, T), Gamma_inv(T, T), Gamma_i_inv(T, T), Sigma_m(l, l), Sigma_e(l, l), tmp_sig_bj(l, l), Sigma_bj(l, l);
  arma::mat V1 (l, l), V2 (l,l);
  arma::mat Ytild(T, n), BGinv(T, l), YX(T, q), XX(q, q);
  
  X2 = sum(pow(X, 2)).t();
  sumX = sum(X).t();
  sumY = sum(Y).t();
  YX = Y.t()*X;
  XX = X.t()*X;
  
  
  // Pour afficher l'avancement du Gibbs
  Rprintf("\nSampler Progress...\n| ");
  for(i=0; i<100; i+=20) Rprintf("%d%%      ",i);
  Rprintf("%d%% |\n|",i);
  
  while(iter<= nIter){          // GIBBS <----------------------------------
    // Rcout << iter << "\n";
    diag_0_Gamma =  arma::ones(T, 1) * (1+rho*rho); diag_0_Gamma(0) = diag_0_Gamma(T-1) = 1;  //arma::repmat( 1+rho*rho, T, 1);
    diag_1_Gamma = -rho* arma::ones(T-1, 1); //arma::repmat( -rho, T-1, 1);
    Gamma_inv = (arma::diagmat(diag_1_Gamma, -1) + arma::diagmat(diag_0_Gamma, 0) + arma::diagmat(diag_1_Gamma, 1))/(1-rho*rho);
    BGinv = B.t() * Gamma_inv;
    
    // Rcout << " 1.0 \n";
    // 1.0_ update alpha:
    var_alpha = as<double>(wrap(n * one.t() * Gamma_inv * one / se2));
    tmp_a = sumY - n*Bt*m - n*Be*e - B*b*sumX;
    mean_alpha = as<double>(wrap( one.t() * Gamma_inv/se2 * tmp_a / var_alpha ));
    alpha = R::rnorm(mean_alpha, sqrt(1/var_alpha));
    
    // Rcout << " 1.1 \n";
    // 1_ update de m:
    Sigma_m = inv(K /tau0 + Bt.t() * Gamma_inv *Bt * n / se2);
    tmp_m = sumY - n*alpha - n * Be * e - B*b*sumX;
    mean_m = Sigma_m* Bt.t() * Gamma_inv *tmp_m / se2;
    m = rmnorm(mean_m, Sigma_m);
    tau0 = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(m.t()*K*m /2 + r)));
    // tau0 = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(m.t()*K*m /se2/2 + r)));
    
    // Rcout << " 1.2 \n";
    // 1_ update e (environnement):
    Sigma_e = inv(K /tau0_e + Be.t() * Gamma_inv * Be * n / se2);
    tmp_e = sumY - n*alpha -n*Bt*m - B*b*sumX;
    mean_e = Sigma_e * Be.t() * Gamma_inv * tmp_e / se2;
    e = rmnorm(mean_e, Sigma_e);
    tau0_e = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(e.t()*K*e /2 + r)));
    
    // Rcout << " 2 \n";
    // 2 update beta and tau2
    //---------------- Dirac on the Slab
    for(j=0; j<q; j++){
      Sigma_bj = inv(K / (se2 * tau2(j)) + as<double>(wrap(X2(j))) *BGinv * B /se2 );
      tmp_mean_bj = YX.col(j) - sumX(j)*one*alpha - sumX(j)*Bt*m - sumX(j)*Be*e - B*b* XX.col(j) + B*b.col(j)* X2(j);
      mean_bj = 1/se2  * Sigma_bj * BGinv * tmp_mean_bj;
      b.col(j) = rmnorm(mean_bj, Sigma_bj);
      tau2(j) = 1/rInvGauss(as<double>(wrap(sqrt(lambda2*se2/(b.col(j).t()*K*b.col(j))))), lambda2);
      
    }
    
    // Rcout << " 3 \n";
    // 3. update lambda2
    lambda2 = R::rgamma((double)(q*l + q)/2+s, (double)(sum(tau2) /2 + r));
    // lambda2 = R::rgamma((double)(q*l + q)/2+s, (double)(l*sum(tau2) /2 + r));
    
    for(i=0; i<n; i++) Ytild.col(i) = (Y.row(i).t() - alpha - Bt*m - Be*e - B*b*X.row(i).t());
    a0 = as<double>(wrap(accu(pow( Ytild.rows(1, T-2), 2))));
    a1 = as<double>(wrap(accu(Ytild.rows(0, T-2) % Ytild.rows(1, T-1))));
    a2 = a0 + as<double>(wrap(accu(pow(Ytild.row(0), 2)))) + as<double>(wrap(accu(pow(Ytild.row(T-1), 2))));
    
    // Rcout << " 4 \n";
    // 4. update rho (Metropolis-Hastings)
    rho_i = rho;
    tmp_i = (rho_i*rho_i * a0 - 2 * rho_i * a1 + a2)/(1-rho_i*rho_i);
    for(ii=0; ii<k; ii++){
      rho_star = R::runif(mmax(-1,rho-rho_tuning), mmin(1,rho+rho_tuning));
      log_Qnew=R::dunif(rho_star, mmax(-1,rho_i-rho_tuning), mmin(1,rho_i+rho_tuning), 1); // proposal distribution Q
      log_Qold=R::dunif(rho_i, mmax(-1,rho_star-rho_tuning), mmin(1,rho_star+rho_tuning), 1);
      
      tmp_star = (rho_star*rho_star * a0 - 2 * rho_star * a1 + a2)/(1-rho_star*rho_star);
      log_ratio = -n/2 *(T-1)* log(1-rho_star*rho_star) - 1/(2*se2)*tmp_star + log_Qold + n/2 * (T-1) * log(1-rho_i*rho_i) + 1/(2*se2)*  tmp_i - log_Qnew;
      
      if(R::runif(0, 1)< exp(log_ratio)) {
        rho_i = rho_star;
        tmp_i = tmp_star;
      }
    }
    rho = rho_i;
    
    // Rcout << " 5 \n";
    // 5. update de sigma^2:
    rate_se2 = tmp_i/2;
    for(j=0; j<q; j++) rate_se2 += as<double>(wrap(b.col(j).t()*K*b.col(j)/(tau2[j]*2)));
    rate_se2 += bb;
    se2 = 1 / R::rgamma(a + (double)(n*T)/2 + (double)(l*q)/2, 1/rate_se2 );
    
    
    //---------------------------------------------------------------------------------------------------
    // Rcout << iter << "\n";
    if(iter>nIter/2){
      chain_m.row(iter_tmp) = m.t();
      chain_e.row(iter_tmp) = e.t();
      for(j=0; j<=q-1; j++){
        chain_b.slice(j).row(iter_tmp) = b.col(j).t();
      }
      chain_tau2.row(iter_tmp) = tau2.t();
      chain_lambda2(iter_tmp) = lambda2;
      chain_tau0(iter_tmp) = tau0;
      chain_tau0_e(iter_tmp) = tau0_e;
      chain_alpha(iter_tmp) = alpha;
      chain_rho(iter_tmp) = rho;
      chain_se2(iter_tmp) = se2;
      iter_tmp ++;
    }
    iter ++;
    
    if(iter/(double)(nIter) >= perc/100.0){
      Rprintf("*");
      perc +=2;
    }
    
  }
  Rprintf("|");
  List chain = List::create( Named("m") = chain_m, Named("e") = chain_e, Named("b") = chain_b, Named("rho") = chain_rho,
                             Named("se2") = chain_se2, Named("tau2") = chain_tau2, Named("tau0") = chain_tau0, 
                                   Named("tau0_e") = chain_tau0_e, Named("alpha") = chain_alpha, Named("lambda2") = chain_lambda2);
  return(chain);
  
}

// [[Rcpp::export]]
Rcpp::IntegerVector sample_int(int n, int min, int max) {
  Rcpp::IntegerVector pool = Rcpp::seq(min, max);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[Rcpp::Range(0, n - 1)];
}

// [[Rcpp::export]]
double sumlogs(double a, double b){
  double M = mmax(a, b);
  return(  M + log(exp(a-M) + exp(b-M))  );
}


//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
List ss_vcm_cpp(arma::mat Y, arma::mat X, List settings, List init, std::string var_gp) {
  //
  //
  //
  //  Input :
  //    Y : matrice (n x T), observation dynamique
  //    X : matrice (n x q), marqueurs moleculaire
  //    init : list, valeurs initiales pour les parametres
  //    settings : list, reglages
  //    sim : liste, parametres simules pour les graphs
  //
  //  Output :
  //    Chain : liste contenant six matrices, chaques matrices contient les echantillonges selont les
  //            lois conditionnelles completes de chaques parametres (m, b, to, rho, se2, lambda2)
  
  int n = Y.n_rows, T = Y.n_cols, q = X.n_cols, j;
  
  // Rcout << " Settings \n";
  // parametre :
  // std::string var_gp = settings["GP_var"];    // 'local' or 'global' 
  Rcout << var_gp << "\n";
  
  int nIter = settings["n.iter"];             // nombre d'iteration de gibbs
  int burnin = settings["burnin"];
  int thin = settings["thin"];
  int k = settings["k"];                      // Nombre d'iteration de M-H pour gamma.head(dg)
  int l = settings["l"];                      // degre des spline + nombre de noeuds
  double a = settings["a"];                   // parametre de sigma^2
  double bb = settings["bb"];                 // parametre de sigma^2
  double s = settings["s"];                   // parametre shape de tau2
  double r = settings["r"];                   // parametre rate de tau2
  double rho_tuning = settings["rho_tuning"]; // Parametre pour le MH de rho
  arma::mat B = settings["B"];                // Matrice de base temp
  arma::mat Bt = settings["Bt"];              // Matrice de base temps (QR decomposition)
  int n_env = settings["n_env"];
  // arma::mat Be = settings["Be"];              // Matrice de base environnement (QR decomposition)
  arma::cube Be  = settings["Be"];
  // arma::mat K = settings["K"]; // matrice des differences (p-spline)
  arma::mat D = settings["D"]; // matrix of finite difference operator
  arma::mat K2 = settings["K2"];                // matrice des differences (p-spline)
  // double s2_a = settings["s2_a"];              // variance de alpha
  double epsilon = settings["epsilon"];
  
  // Rcout << " Initialisation \n";
  // Initialisation :
  arma::vec m = init["m"];			              //  m, intercept, vecteur de longueur v (degre de l'approximation par polynome de legendre
  arma::mat e = init["e"];			              //  e, environment, vecteur de longueur v (degre de l'approximation par polynome de legendre
  arma::mat b = init["b"];			              // b, coefficient dynamique, matrice de taile (T x q)
  double rho  = init["rho"];                   // rho, parametre auto-regressif sur la matrice de variance residuelle
  double se2  = init["se2"];	                  // sigma^2, variance residuelle, scalaire
  double pi   = init["pi"];                 // proportion de variables selectionnees a priori
  arma::vec g = init["g"];                // variables indicatrices gamma
  double tau0 = init["tau0"];                 // variance de la loi a priori sur m
  arma::vec tau0_e = init["tau0_e"];             // variance de la loi a priori sur e
  double alpha = init["alpha"];
  
  arma::vec tau2 = arma::ones(q);              // variance spike and slab sur beta
  arma::vec phi_tau2 = 2*arma::ones(arma::size(tau2));
  arma::mat xi = arma::ones(arma::size(D)(0), q);   // local variance parameter on the gaussian process prior of the effect b
  arma::mat phi_xi = 2*arma::ones(arma::size(xi));
  if(var_gp == "global"){
    arma::vec tau2 = init["tau2"];              // variance spike and slab sur beta
  }
  if(var_gp == "local"){
    arma::mat xi = init["xi"];   // local variance parameter on the gaussian process prior of the effect b
  }
  
  // declaration des matrices de sorties
  // arma::mat chain_m (nIter/2+1, l-1); //colnames(chain.m) = paste0("m.", 1:v);
  // arma::cube chain_e (nIter/2+1, l-1, n_env);
  // arma::cube chain_b (nIter/2+1, l, q);
  // arma::mat chain_lambda1 (nIter/2+1, q);
  // arma::mat chain_tau2 (nIter/2+1,q);
  // arma::vec chain_rho (nIter/2+1);
  // arma::vec chain_se2 (nIter/2+1);
  // arma::mat chain_g (nIter/2+1, q);
  // arma::vec chain_tau0 (nIter/2+1);
  // arma::mat chain_tau0_e (nIter/2+1, n_env);
  // arma::vec chain_alpha (nIter/2+1);
  
  arma::mat chain_m ((nIter-burnin)/thin, l-1); //colnames(chain.m) = paste0("m.", 1:v);
  arma::cube chain_e ((nIter-burnin)/thin, l-1, n_env);
  arma::cube chain_b ((nIter-burnin)/thin, arma::size(B)(1), q);
  arma::mat chain_lambda1 ((nIter-burnin)/thin, q);
  arma::mat chain_tau2 ((nIter-burnin)/thin, q);
  arma::vec chain_rho ((nIter-burnin)/thin);
  arma::vec chain_se2 ((nIter-burnin)/thin);
  arma::vec chain_pi ((nIter-burnin)/thin);
  arma::mat chain_g ((nIter-burnin)/thin, q);
  arma::vec chain_tau0 ((nIter-burnin)/thin);
  arma::mat chain_tau0_e ((nIter-burnin)/thin, n_env);
  arma::vec chain_alpha ((nIter-burnin)/thin);
  arma::cube chain_xi ((nIter-burnin)/thin, arma::size(D)(0), q);
  
  // Rcout << "0.1 \n";
  // chain_m.row(0) = m.t();
  // for(j=0; j < n_env; j++){
  //   chain_e.slice(j).row(0) = e.col(j).t();
  // }
  // for(j=0; j<=q-1; j++){
  //   chain_b.slice(j).row(0) = b.col(j).t();
  // }
  // chain_tau2.row(0) = tau2.t();
  // chain_rho(0) = rho;
  // chain_se2(0) = se2;
  // chain_g.row(0) = g.t();
  // chain_tau0(0) = tau0;
  // chain_tau0_e.row(0) = tau0_e.t();
  // chain_alpha(0) = alpha;
  
  // Declaration des variables internes :
  IntegerVector ind = Range(0, q-1);
  int iter = 0, i, ii, perc=2, iter_tmp = 0;
  double rho_i, rho_star, tmp_i, tmp_star, log_Qnew, log_Qold, rate_se2, log_ratio, log_R, ratio, a0, a1, a2, var_alpha, mean_alpha;
  arma::vec tmp_m(T), tmp_a(T), mean_m(l-1), tmp_e(T), mean_e(l-1), mean_bj(l), tmp_mean_bj(T), tmp(T), tmp_mean_cj(T), diag_0_Gamma(T), diag_1_Gamma(T);
  arma::vec tmp2(T), s_tmp(l), X2 (q), sumX(q), sumY(T), b2(q), one = arma::ones(T), Bee (T);
  arma::mat Gamma(T, T),Gamma_i(T, T), Gamma_star(T, T), Gamma_inv(T, T), Gamma_i_inv(T, T), Sigma_m(l-1, l-1), Sigma_e(l-1, l-1), tmp_sig_bj(l, l), Sigma_bj(l, l);
  arma::mat V1 (l, l), V2 (l,l);
  arma::mat Ytild(T, n), BGinv(T, l), YX(T, q), XX(q, q);
  Rcpp::IntegerVector rand(q);
  arma::vec one2 = arma::ones(arma::size(D)(1));
  // Rcout << "calculus \n" ;
  X2 = sum(pow(X, 2)).t();
  sumX = sum(X).t();
  sumY = sum(Y).t();
  YX = Y.t()*X;
  XX = X.t()*X;
  
  Bee = Be.slice(0) * e.col(0);
  for(i=1; i<n_env; i++) Bee += Be.slice(i) * e.col(i);
  
  
  // Pour afficher l'avancement du Gibbs
  Rprintf("\nSampler Progress...\n| ");
  for(i=0; i<100; i+=20) Rprintf("%d%%      ",i);
  Rprintf("%d%% |\n|",i);
  
  while(iter< nIter){          // GIBBS <----------------------------------
    // Rcout << iter << "\n";
    diag_0_Gamma =  arma::ones(T, 1) * (1+rho*rho); diag_0_Gamma(0) = diag_0_Gamma(T-1) = 1;  //arma::repmat( 1+rho*rho, T, 1);
    diag_1_Gamma = -rho* arma::ones(T-1, 1); //arma::repmat( -rho, T-1, 1);
    Gamma_inv = (arma::diagmat(diag_1_Gamma, -1) + arma::diagmat(diag_0_Gamma, 0) + arma::diagmat(diag_1_Gamma, 1))/(1-rho*rho);
    BGinv = B.t() * Gamma_inv;
    
    // Rcout << " 1.0 \n";
    // 1.0_ update alpha:
    var_alpha = as<double>(wrap(n * one.t() * Gamma_inv * one / se2));
    tmp_a = sumY - n*Bt*m - n*Bee - B*b*sumX;
    mean_alpha = as<double>(wrap( one.t() * Gamma_inv/se2 * tmp_a / var_alpha ));
    alpha = R::rnorm(mean_alpha, sqrt(1/var_alpha));
    
    // Rcout << " 1.1 \n";
    // 1_ update m:
    Sigma_m = inv(K2 /tau0 + Bt.t() * Gamma_inv *Bt * n / se2);
    tmp_m = sumY - n*alpha - n * Bee - B*b*sumX;
    mean_m = Sigma_m * Bt.t() * Gamma_inv * tmp_m / se2;
    m = rmnorm(mean_m, Sigma_m);
    tau0 = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(m.t()*K2*m /2 + r)));
    // tau0 = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(m.t()*K2*m /se2/2 + r)));
    
    // Rcout << " 1.2 \n";
    // 1_ update e (environnement):
    for(j=0; j< n_env; j++){
      Sigma_e = inv(K2 /tau0_e(j) + Be.slice(j).t() * Gamma_inv * Be.slice(j) * n / se2);
      tmp_e = sumY - n*alpha -n*Bt*m - n*(Bee - Be.slice(j) * e.col(j)) - B*b*sumX;
      mean_e = Sigma_e * Be.slice(j).t() * Gamma_inv * tmp_e / se2;
      e.col(j) = rmnorm(mean_e, Sigma_e);
      tau0_e(j) = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(e.col(j).t()*K2*e.col(j) /2 + r)));
      // tau0_e = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(e.t()*K2*e /se2/2 + r)));
      Bee = Be.slice(0) * e.col(0);
      for(i=1; i<n_env; i++) Bee += Be.slice(i) * e.col(i);
    }
    
    
    // Rcout << " 2 \n";
    // 2 update gamma, beta et tau2
    //---------------- Dirac sur le Slab
    rand = sample_int(q, 0, q-1);
    for(i=0; i<q; i++){
      // Rcout << i << "  \n";
      
      j = rand(i);
      // Sigma_bj = inv(K / (se2 * tau2(j)) + as<double>(wrap(X2(j))) *BGinv * B /se2 );
      Sigma_bj = inv((D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2)) / (se2 * tau2(j)) + as<double>(wrap(X2(j))) *BGinv * B /se2 );
      tmp_mean_bj = YX.col(j) - sumX(j)*one*alpha - sumX(j)*Bt*m - sumX(j)*Bee - B*b* XX.col(j) + B*b.col(j)* X2(j);
      // Rcout << " 2.2 \n";
      
      // R = as<double>(wrap( pi/(1-pi)* sqrt(det(Sigma_bj)) / pow(se2*tau2(j), l/2) * exp(tmp_mean_bj.t() * Gamma_inv/se2  * B * Sigma_bj * B.t() * Gamma_inv/se2 * tmp_mean_bj / 2)));
      // if(isinf(R)) ratio = 1; else ratio = R/(1+R);
      log_R = as<double>(wrap(  log( pi/(1-pi) * sqrt(det(Sigma_bj)) / pow(se2*tau2(j), l/2) ) + tmp_mean_bj.t() * Gamma_inv/se2  * B * Sigma_bj * B.t() * Gamma_inv/se2 * tmp_mean_bj / 2  ));
      // log_R = as<double>(wrap(  log( pi/(1-pi)) + log(det(D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2)))/2 + log(det(Sigma_bj))/2 -  l/2* log(tau2(j)*se2)  + tmp_mean_bj.t() * Gamma_inv/se2  * B * Sigma_bj * B.t() * Gamma_inv/se2 * tmp_mean_bj / 2  ));
      ratio = exp(log_R - sumlogs(1, log_R));
      // Rcout << " 2.3 \n";
      
      if(arma::randu()< ratio){
        g(j) = 1;
        mean_bj = 1/se2  * Sigma_bj * BGinv * tmp_mean_bj;
        // Rcout << " 2.4 \n";
        
        b.col(j) = rmnorm(mean_bj, Sigma_bj);
        if(var_gp == "global"){
          tau2(j) = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(b.col(j).t()*(D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2))*b.col(j) /se2/2 + r)));
          // phi_tau2(j) = 1/R::rgamma(1, 1/as<double>(wrap(1+1/tau2(j))));
        }
        if(var_gp == "local"){
          for(ii=0; ii<arma::size(D)(0); ii++){
            xi(ii, j) = R::rgamma((double)(1)/2 + s, 1/as<double>(wrap(arma::pow(D.row(ii)*b.col(j), 2) /tau2(j)/se2/2 + r)));
            // phi_xi(ii, j) = 1/R::rgamma(1, 1/as<double>(wrap(1+1/xi(ii, j))));
          }
        }
      }
      else{
        g(j) = 0;
        // Rcout << " 2.5 \n";
        
        b.col(j) = arma::zeros(arma::size(B)(1));
        tau2(j) = 1; //1/R::rgamma(s, 1/r);
        xi.col(j) = arma::zeros(arma::size(D)(0));
        // for(ii=0; ii<arma::size(D)(0); ii++){
        //   xi(ii, j) = R::rgamma((double)(1)/2 + s, 1/as<double>(wrap(arma::pow(D.row(ii)*b.col(j), 2) /tau2(j)/se2/2 + r)));
        // }
      }
    }
    
    pi = R::rbeta(1+sum(g), 1+q-sum(g));
    
    for(i=0; i<n; i++) Ytild.col(i) = (Y.row(i).t() - alpha - Bt*m - Bee - B*b*X.row(i).t());
    a0 = as<double>(wrap(accu(pow( Ytild.rows(1, T-2), 2))));
    a1 = as<double>(wrap(accu(Ytild.rows(0, T-2) % Ytild.rows(1, T-1))));
    a2 = a0 + as<double>(wrap(accu(pow(Ytild.row(0), 2)))) + as<double>(wrap(accu(pow(Ytild.row(T-1), 2))));
    
    // Rcout << " 3 \n";
    // 3. update rho (Metropolis-Hastings)
    rho_i = rho;
    tmp_i = (rho_i*rho_i * a0 - 2 * rho_i * a1 + a2)/(1-rho_i*rho_i);
    for(ii=0; ii<k; ii++){
      rho_star = R::runif(mmax(-1,rho-rho_tuning), mmin(1,rho+rho_tuning));
      log_Qnew=R::dunif(rho_star, mmax(-1, rho_i-rho_tuning), mmin(1, rho_i+rho_tuning), 1); // proposal distribution Q
      log_Qold=R::dunif(rho_i, mmax(-1, rho_star-rho_tuning), mmin(1, rho_star+rho_tuning), 1);
      
      tmp_star = (rho_star*rho_star * a0 - 2 * rho_star * a1 + a2)/(1-rho_star*rho_star);
      log_ratio = -n/2 *(T-1)* log(1-rho_star*rho_star) - 1/(2*se2)*tmp_star + log_Qold + n/2 * (T-1) * log(1-rho_i*rho_i) + 1/(2*se2)*  tmp_i - log_Qnew;
      
      if(R::runif(0, 1)< exp(log_ratio)) {
        rho_i = rho_star;
        tmp_i = tmp_star;
      }
    }
    rho = rho_i;
    
    // Rcout << " 4 \n";
    // 4. update sigma^2:
    rate_se2 = tmp_i/2;
    for(j=0; j<q; j++) rate_se2 += as<double>(wrap(b.col(j).t()*(D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2))*b.col(j)/(tau2[j]*2)));
    rate_se2 += bb;
    se2 = 1 / R::rgamma(a + (double)(n*T)/2 + (double)(l*sum(g))/2, 1/rate_se2 );
    
    
    //---------------------------------------------------------------------------------------------------
    // Rcout << iter << "\n";
    if(iter>=burnin & (iter - burnin) % thin == 0){
      chain_m.row(iter_tmp) = m.t();
      for(j=0; j< n_env; j++){
        chain_e.slice(j).row(iter_tmp) = e.col(j).t();
      }
      chain_g.row(iter_tmp) = g.t();
      for(j=0; j<=q-1; j++){
        chain_b.slice(j).row(iter_tmp) = b.col(j).t();
      }
      for(j=0; j<=q-1; j++){
        chain_xi.slice(j).row(iter_tmp) = xi.col(j).t();
      }
      chain_tau2.row(iter_tmp) = tau2.t();
      chain_tau0(iter_tmp) = tau0;
      chain_tau0_e.row(iter_tmp) = tau0_e.t();
      chain_alpha(iter_tmp) = alpha;
      chain_rho(iter_tmp) = rho;
      chain_se2(iter_tmp) = se2;
      chain_pi(iter_tmp) = pi;
      iter_tmp ++;
    }
    iter ++;
    
    if(iter/(double)(nIter) >= perc/100.0){
      Rprintf("*");
      perc +=2;
    }
    
  }
  Rprintf("| \n");
  List chain = List::create(Named("g") = chain_g, Named("m") = chain_m, 
                            Named("e") = chain_e, Named("b") = chain_b, 
                            Named("rho") = chain_rho, Named("pi") = chain_pi,
                            Named("se2") = chain_se2, Named("tau2") = chain_tau2, 
                            Named("tau0") = chain_tau0, Named("tau0_e") = chain_tau0_e, 
                            Named("alpha") = chain_alpha, Named("xi") = chain_xi);
  return(chain);
  
}

//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
List ss_vcm_without_rho_cpp(arma::mat Y, arma::mat X, List settings, List init, std::string var_gp) {
  //
  //
  //
  //  Input :
  //    Y : matrice (n x T), observation dynamique
  //    X : matrice (n x q), marqueurs moleculaire
  //    init : list, valeurs initiales pour les parametres
  //    settings : list, reglages
  //    sim : liste, parametres simules pour les graphs
  //
  //  Output :
  //    Chain : liste contenant six matrices, chaques matrices contient les echantillonges selont les
  //            lois conditionnelles completes de chaques parametres (m, b, to, rho, se2, lambda2)
  
  int n = Y.n_rows, T = Y.n_cols, q = X.n_cols, j;
  
  // Rcout << " Settings \n";
  // parametre :
  // std::string var_gp = settings["GP_var"];    // 'local' or 'global' 
  // Rcout << var_gp << "\n";
  
  int nIter = settings["n.iter"];             // nombre d'iteration de gibbs
  int burnin = settings["burnin"];
  int thin = settings["thin"];
  int k = settings["k"];                      // Nombre d'iteration de M-H pour gamma.head(dg)
  int l = settings["l"];                      // degre des spline + nombre de noeuds
  double a = settings["a"];                   // parametre de sigma^2
  double bb = settings["bb"];                 // parametre de sigma^2
  double s = settings["s"];                   // parametre shape de tau2
  double r = settings["r"];                   // parametre rate de tau2
  double rho_tuning = settings["rho_tuning"]; // Parametre pour le MH de rho
  arma::mat B = settings["B"];                // Matrice de base temp
  arma::mat Bt = settings["Bt"];              // Matrice de base temps (QR decomposition)
  int n_env = settings["n_env"];
  // arma::mat Be = settings["Be"];              // Matrice de base environnement (QR decomposition)
  arma::cube Be  = settings["Be"];
  // arma::mat K = settings["K"]; // matrice des differences (p-spline)
  arma::mat D = settings["D"]; // matrix of finite difference operator
  arma::mat K2 = settings["K2"];                // matrice des differences (p-spline)
  // double s2_a = settings["s2_a"];              // variance de alpha
  double epsilon = settings["epsilon"];
  
  // Rcout << " Initialisation \n";
  // Initialisation :
  arma::vec m = init["m"];			              //  m, intercept, vecteur de longueur v (degre de l'approximation par polynome de legendre
  arma::mat e = init["e"];			              //  e, environment, vecteur de longueur v (degre de l'approximation par polynome de legendre
  arma::mat b = init["b"];			              // b, coefficient dynamique, matrice de taile (T x q)
  double rho  = init["rho"];                   // rho, parametre auto-regressif sur la matrice de variance residuelle
  double se2  = init["se2"];	                  // sigma^2, variance residuelle, scalaire
  double pi   = init["pi"];                 // proportion de variables selectionnees a priori
  arma::vec g = init["g"];                // variables indicatrices gamma
  double tau0 = init["tau0"];                 // variance de la loi a priori sur m
  arma::vec tau0_e = init["tau0_e"];             // variance de la loi a priori sur e
  double alpha = init["alpha"];
  
  arma::vec tau2 = arma::ones(q);              // variance spike and slab sur beta
  arma::vec phi_tau2 = 2*arma::ones(arma::size(tau2));
  arma::mat xi = arma::ones(arma::size(D)(0), q);   // local variance parameter on the gaussian process prior of the effect b
  arma::mat phi_xi = 2*arma::ones(arma::size(xi));
  if(var_gp == "global"){
    arma::vec tau2 = init["tau2"];              // variance spike and slab sur beta
  }
  if(var_gp == "local"){
    arma::mat xi = init["xi"];   // local variance parameter on the gaussian process prior of the effect b
  }
  // Rcout << " declaration des matrices de sorties \n";
  // declaration des matrices de sorties
  // arma::mat chain_m (nIter/2+1, l-1); //colnames(chain.m) = paste0("m.", 1:v);
  // arma::cube chain_e (nIter/2+1, l-1, n_env);
  // arma::cube chain_b (nIter/2+1, l, q);
  // arma::mat chain_lambda1 (nIter/2+1, q);
  // arma::mat chain_tau2 (nIter/2+1,q);
  // arma::vec chain_rho (nIter/2+1);
  // arma::vec chain_se2 (nIter/2+1);
  // arma::mat chain_g (nIter/2+1, q);
  // arma::vec chain_tau0 (nIter/2+1);
  // arma::mat chain_tau0_e (nIter/2+1, n_env);
  // arma::vec chain_alpha (nIter/2+1);
  
  arma::mat chain_m ((nIter-burnin)/thin, l-1); //colnames(chain.m) = paste0("m.", 1:v);
  arma::cube chain_e ((nIter-burnin)/thin, l-1, n_env);
  arma::cube chain_b ((nIter-burnin)/thin, arma::size(B)(1), q);
  arma::mat chain_lambda1 ((nIter-burnin)/thin, q);
  arma::mat chain_tau2 ((nIter-burnin)/thin, q);
  arma::vec chain_rho ((nIter-burnin)/thin);
  arma::vec chain_se2 ((nIter-burnin)/thin);
  arma::vec chain_pi ((nIter-burnin)/thin);
  arma::mat chain_g ((nIter-burnin)/thin +1, q);
  arma::vec chain_tau0 ((nIter-burnin)/thin);
  arma::mat chain_tau0_e ((nIter-burnin)/thin, n_env);
  arma::vec chain_alpha ((nIter-burnin)/thin);
  arma::cube chain_xi ((nIter-burnin)/thin, arma::size(D)(0), q);
  
  // Rcout << "0.1 \n";
  // chain_m.row(0) = m.t();
  // for(j=0; j < n_env; j++){
  //   chain_e.slice(j).row(0) = e.col(j).t();
  // }
  // for(j=0; j<=q-1; j++){
  //   chain_b.slice(j).row(0) = b.col(j).t();
  // }
  // chain_tau2.row(0) = tau2.t();
  // chain_rho(0) = rho;
  // chain_se2(0) = se2;
  chain_g.row(0) = g.t();
  // chain_tau0(0) = tau0;
  // chain_tau0_e.row(0) = tau0_e.t();
  // chain_alpha(0) = alpha;
  
  // Rcout << "Declaration des variables internes \n" ;
  // Declaration des variables internes :
  IntegerVector ind = Range(0, q-1);
  int iter = 0, i, ii, perc=2, iter_tmp = 0;
  double rho_i, rho_star, tmp_i, tmp_star, log_Qnew, log_Qold, rate_se2, log_ratio, log_R, ratio, a0, a1, a2, var_alpha, mean_alpha;
  arma::vec tmp_m(T), tmp_a(T), mean_m(l-1), tmp_e(T), mean_e(l-1), mean_bj(l), tmp_mean_bj(T), tmp(T), tmp_mean_cj(T), diag_0_Gamma(T), diag_1_Gamma(T);
  arma::vec tmp2(T), s_tmp(l), X2 (q), sumX(q), sumY(T), b2(q), one = arma::ones(T), Bee (T);
  arma::mat Gamma(T, T),Gamma_i(T, T), Gamma_star(T, T), Gamma_inv(T, T), Gamma_i_inv(T, T), Sigma_m(l-1, l-1), Sigma_e(l-1, l-1), tmp_sig_bj(l, l), Sigma_bj(l, l);
  arma::mat V1 (l, l), V2 (l,l);
  arma::mat Ytild(T, n), BGinv(T, l), YX(T, q), XX(q, q);
  Rcpp::IntegerVector rand(q);
  arma::vec one2 = arma::ones(arma::size(D)(1));
  // Rcout << "calculus \n" ;
  X2 = sum(pow(X, 2)).t();
  sumX = sum(X).t();
  sumY = sum(Y).t();
  YX = Y.t()*X;
  XX = X.t()*X;
  
  Bee = Be.slice(0) * e.col(0);
  for(i=1; i<n_env; i++) Bee += Be.slice(i) * e.col(i);
  
  
  // Pour afficher l'avancement du Gibbs
  Rprintf("\nSampler Progress...\n| ");
  for(i=0; i<100; i+=20) Rprintf("%d%%      ",i);
  Rprintf("%d%% |\n|",i);
  
  while(iter< nIter){          // GIBBS <----------------------------------
    // Rcout << iter << "\n";
    diag_0_Gamma =  arma::ones(T, 1) * (1+rho*rho); diag_0_Gamma(0) = diag_0_Gamma(T-1) = 1;  //arma::repmat( 1+rho*rho, T, 1);
    diag_1_Gamma = -rho* arma::ones(T-1, 1); //arma::repmat( -rho, T-1, 1);
    Gamma_inv = (arma::diagmat(diag_1_Gamma, -1) + arma::diagmat(diag_0_Gamma, 0) + arma::diagmat(diag_1_Gamma, 1))/(1-rho*rho);
    BGinv = B.t() * Gamma_inv;
    
    // Rcout << " 1.0 \n";
    // 1.0_ update alpha:
    var_alpha = as<double>(wrap(n * one.t() * Gamma_inv * one / se2));
    tmp_a = sumY - n*Bt*m - n*Bee - B*b*sumX;
    mean_alpha = as<double>(wrap( one.t() * Gamma_inv/se2 * tmp_a / var_alpha ));
    alpha = R::rnorm(mean_alpha, sqrt(1/var_alpha));
    
    // Rcout << " 1.1 \n";
    // 1_ update m:
    Sigma_m = inv(K2 /tau0 + Bt.t() * Gamma_inv *Bt * n / se2);
    tmp_m = sumY - n*alpha - n * Bee - B*b*sumX;
    mean_m = Sigma_m * Bt.t() * Gamma_inv * tmp_m / se2;
    m = rmnorm(mean_m, Sigma_m);
    tau0 = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(m.t()*K2*m /2 + r)));
    // tau0 = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(m.t()*K2*m /se2/2 + r)));
    
    // Rcout << " 1.2 \n";
    // 1_ update e (environnement):
    for(j=0; j< n_env; j++){
      Sigma_e = inv(K2 /tau0_e(j) + Be.slice(j).t() * Gamma_inv * Be.slice(j) * n / se2);
      tmp_e = sumY - n*alpha -n*Bt*m - n*(Bee - Be.slice(j) * e.col(j)) - B*b*sumX;
      mean_e = Sigma_e * Be.slice(j).t() * Gamma_inv * tmp_e / se2;
      e.col(j) = rmnorm(mean_e, Sigma_e);
      tau0_e(j) = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(e.col(j).t()*K2*e.col(j) /2 + r)));
      // tau0_e = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(e.t()*K2*e /se2/2 + r)));
      Bee = Be.slice(0) * e.col(0);
      for(i=1; i<n_env; i++) Bee += Be.slice(i) * e.col(i);
    }
    
    
    // Rcout << " 2 \n";
    // 2 update gamma, beta et tau2
    //---------------- Dirac sur le Slab
    rand = sample_int(q, 0, q-1);
    for(i=0; i<q; i++){
      // Rcout << i << "  \n";
      
      j = rand(i);
      // Sigma_bj = inv(K / (se2 * tau2(j)) + as<double>(wrap(X2(j))) *BGinv * B /se2 );
      Sigma_bj = inv((D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2)) / (se2 * tau2(j)) + as<double>(wrap(X2(j))) *BGinv * B /se2 );
      tmp_mean_bj = YX.col(j) - sumX(j)*one*alpha - sumX(j)*Bt*m - sumX(j)*Bee - B*b* XX.col(j) + B*b.col(j)* X2(j);
      // Rcout << " 2.2 \n";
      
      // R = as<double>(wrap( pi/(1-pi)* sqrt(det(Sigma_bj)) / pow(se2*tau2(j), l/2) * exp(tmp_mean_bj.t() * Gamma_inv/se2  * B * Sigma_bj * B.t() * Gamma_inv/se2 * tmp_mean_bj / 2)));
      // if(isinf(R)) ratio = 1; else ratio = R/(1+R);
      log_R = as<double>(wrap(  log( pi/(1-pi) * sqrt(det(Sigma_bj)) / pow(se2*tau2(j), l/2) ) + tmp_mean_bj.t() * Gamma_inv/se2  * B * Sigma_bj * B.t() * Gamma_inv/se2 * tmp_mean_bj / 2  ));
      ratio = exp(log_R - sumlogs(1, log_R));
      // Rcout << " 2.3 \n";
      
      if(arma::randu()< ratio){
        g(j) = 1;
        mean_bj = 1/se2  * Sigma_bj * BGinv * tmp_mean_bj;
        // Rcout << " 2.4 \n";
        
        b.col(j) = rmnorm(mean_bj, Sigma_bj);
        if(var_gp == "global"){
          tau2(j) = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(b.col(j).t()*(D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2))*b.col(j) /se2/2 + r)));
          // phi_tau2(j) = 1/R::rgamma(1, 1/as<double>(wrap(1+1/tau2(j))));
        }
        if(var_gp == "local"){
          for(ii=0; ii<arma::size(D)(0); ii++){
            xi(ii, j) = R::rgamma((double)(1)/2 + s, 1/as<double>(wrap(arma::pow(D.row(ii)*b.col(j), 2) /tau2(j)/se2/2 + r)));
            // phi_xi(ii, j) = 1/R::rgamma(1, 1/as<double>(wrap(1+1/xi(ii, j))));
          }
        }
      }
      else{
        g(j) = 0;
        // Rcout << " 2.5 \n";
        
        b.col(j) = arma::zeros(arma::size(B)(1));
        tau2(j) = 1; //1/R::rgamma(s, 1/r);
        xi.col(j) = arma::zeros(arma::size(D)(0));
        // for(ii=0; ii<arma::size(D)(0); ii++){
        //   xi(ii, j) = R::rgamma((double)(1)/2 + s, 1/as<double>(wrap(arma::pow(D.row(ii)*b.col(j), 2) /tau2(j)/se2/2 + r)));
        // }
      }
    }
    
    pi = R::rbeta(1+sum(g), 1+q-sum(g));
    
    for(i=0; i<n; i++) Ytild.col(i) = (Y.row(i).t() - alpha - Bt*m - Bee - B*b*X.row(i).t());
    a0 = as<double>(wrap(accu(pow( Ytild.rows(1, T-2), 2))));
    a1 = as<double>(wrap(accu(Ytild.rows(0, T-2) % Ytild.rows(1, T-1))));
    a2 = a0 + as<double>(wrap(accu(pow(Ytild.row(0), 2)))) + as<double>(wrap(accu(pow(Ytild.row(T-1), 2))));
    
    // Rcout << " 3 \n";
    // 3. update rho (Metropolis-Hastings)
    rho_i = rho;
    tmp_i = (rho_i*rho_i * a0 - 2 * rho_i * a1 + a2)/(1-rho_i*rho_i);
    for(ii=0; ii<k; ii++){
      rho_star = 0; // R::runif(mmax(-1,rho-rho_tuning), mmin(1,rho+rho_tuning));
      log_Qnew=R::dunif(rho_star, mmax(-1, rho_i-rho_tuning), mmin(1, rho_i+rho_tuning), 1); // proposal distribution Q
      log_Qold=R::dunif(rho_i, mmax(-1, rho_star-rho_tuning), mmin(1, rho_star+rho_tuning), 1);
      
      tmp_star = (rho_star*rho_star * a0 - 2 * rho_star * a1 + a2)/(1-rho_star*rho_star);
      log_ratio = -n/2 *(T-1)* log(1-rho_star*rho_star) - 1/(2*se2)*tmp_star + log_Qold + n/2 * (T-1) * log(1-rho_i*rho_i) + 1/(2*se2)*  tmp_i - log_Qnew;
      
      if(R::runif(0, 1)< exp(log_ratio)) {
        rho_i = rho_star;
        tmp_i = tmp_star;
      }
    }
    // rho = rho_i;
    rho = 0;
    
    
    
    // Rcout << " 4 \n";
    // 4. update sigma^2:
    rate_se2 = tmp_i/2;
    for(j=0; j<q; j++) rate_se2 += as<double>(wrap(b.col(j).t()*(D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2))*b.col(j)/(tau2[j]*2)));
    rate_se2 += bb;
    se2 = 1 / R::rgamma(a + (double)(n*T)/2 + (double)(l*sum(g))/2, 1/rate_se2 );
    
    
    //---------------------------------------------------------------------------------------------------
    // Rcout << iter << "\n";
    if(iter>=burnin & (iter - burnin) % thin == 0){
      chain_m.row(iter_tmp) = m.t();
      for(j=0; j< n_env; j++){
        chain_e.slice(j).row(iter_tmp) = e.col(j).t();
      }
      chain_g.row(iter_tmp+1) = g.t();
      for(j=0; j<=q-1; j++){
        chain_b.slice(j).row(iter_tmp) = b.col(j).t();
      }
      for(j=0; j<=q-1; j++){
        chain_xi.slice(j).row(iter_tmp) = xi.col(j).t();
      }
      chain_tau2.row(iter_tmp) = tau2.t();
      chain_tau0(iter_tmp) = tau0;
      chain_tau0_e.row(iter_tmp) = tau0_e.t();
      chain_alpha(iter_tmp) = alpha;
      chain_rho(iter_tmp) = rho;
      chain_se2(iter_tmp) = se2;
      chain_pi(iter_tmp) = pi;
      iter_tmp ++;
    }
    iter ++;
    
    if(iter/(double)(nIter) >= perc/100.0){
      Rprintf("*");
      perc +=2;
    }
    
  }
  Rprintf("| \n");
  List chain = List::create(Named("g") = chain_g, Named("m") = chain_m, 
                            Named("e") = chain_e, Named("b") = chain_b, 
                            Named("rho") = chain_rho, Named("pi") = chain_pi,
                            Named("se2") = chain_se2, Named("tau2") = chain_tau2, 
                            Named("tau0") = chain_tau0, Named("tau0_e") = chain_tau0_e, 
                            Named("alpha") = chain_alpha, Named("xi") = chain_xi);
  return(chain);
  
}







//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
List vcm_cpp(arma::mat Y, arma::mat X, List settings, List init, std::string var_gp) {
  //
  //
  //
  //  Input :
  //    Y : matrice (n x T), observation dynamique
  //    X : matrice (n x q), marqueurs moleculaire
  //    init : list, valeurs initiales pour les parametres
  //    settings : list, reglages
  //    sim : liste, parametres simules pour les graphs
  //
  //  Output :
  //    Chain : liste contenant six matrices, chaques matrices contient les echantillonges selont les
  //            lois conditionnelles completes de chaques parametres (m, b, to, rho, se2, lambda2)
  
  int n = Y.n_rows, T = Y.n_cols, q = X.n_cols, j;
  
  // Rcout << " Settings \n";
  // parametre :
  // std::string var_gp = settings["GP_var"];    // 'local' or 'global' 
  // Rcout << var_gp << "\n";
  
  int nIter = settings["n.iter"];             // nombre d'iteration de gibbs
  int burnin = settings["burnin"];
  int thin = settings["thin"];
  int k = settings["k"];                      // Nombre d'iteration de M-H pour gamma.head(dg)
  int l = settings["l"];                      // degre des spline + nombre de noeuds
  double a = settings["a"];                   // parametre de sigma^2
  double bb = settings["bb"];                 // parametre de sigma^2
  double s = settings["s"];                   // parametre shape de tau2
  double r = settings["r"];                   // parametre rate de tau2
  double rho_tuning = settings["rho_tuning"]; // Parametre pour le MH de rho
  arma::mat B = settings["B"];                // Matrice de base temp
  arma::mat Bt = settings["Bt"];              // Matrice de base temps (QR decomposition)
  int n_env = settings["n_env"];
  // arma::mat Be = settings["Be"];              // Matrice de base environnement (QR decomposition)
  arma::cube Be  = settings["Be"];
  // arma::mat K = settings["K"]; // matrice des differences (p-spline)
  arma::mat D = settings["D"]; // matrix of finite difference operator
  arma::mat K2 = settings["K2"];                // matrice des differences (p-spline)
  // double s2_a = settings["s2_a"];              // variance de alpha
  double epsilon = settings["epsilon"];
  
  // Rcout << " Initialisation \n";
  // Initialisation :
  arma::vec m = init["m"];			              //  m, intercept, vecteur de longueur v (degre de l'approximation par polynome de legendre
  arma::mat e = init["e"];			              //  e, environment, vecteur de longueur v (degre de l'approximation par polynome de legendre
  arma::mat b = init["b"];			              // b, coefficient dynamique, matrice de taile (T x q)
  double rho  = init["rho"];                   // rho, parametre auto-regressif sur la matrice de variance residuelle
  double se2  = init["se2"];	                  // sigma^2, variance residuelle, scalaire
  double pi   = 1;                 // proportion de variables selectionnees a priori
  arma::vec g = init["g"];                // variables indicatrices gamma
  double tau0 = init["tau0"];                 // variance de la loi a priori sur m
  arma::vec tau0_e = init["tau0_e"];             // variance de la loi a priori sur e
  double alpha = init["alpha"];
  
  arma::vec tau2 = arma::ones(q);              // variance spike and slab sur beta
  arma::vec phi_tau2 = 2*arma::ones(arma::size(tau2));
  arma::mat xi = arma::ones(arma::size(D)(0), q);   // local variance parameter on the gaussian process prior of the effect b
  arma::mat phi_xi = 2*arma::ones(arma::size(xi));
  if(var_gp == "global"){
    arma::vec tau2 = init["tau2"];              // variance spike and slab sur beta
  }
  if(var_gp == "local"){
    arma::mat xi = init["xi"];   // local variance parameter on the gaussian process prior of the effect b
  }
  // Rcout << " declaration des matrices de sorties \n";
  // declaration des matrices de sorties
  // arma::mat chain_m (nIter/2+1, l-1); //colnames(chain.m) = paste0("m.", 1:v);
  // arma::cube chain_e (nIter/2+1, l-1, n_env);
  // arma::cube chain_b (nIter/2+1, l, q);
  // arma::mat chain_lambda1 (nIter/2+1, q);
  // arma::mat chain_tau2 (nIter/2+1,q);
  // arma::vec chain_rho (nIter/2+1);
  // arma::vec chain_se2 (nIter/2+1);
  // arma::mat chain_g (nIter/2+1, q);
  // arma::vec chain_tau0 (nIter/2+1);
  // arma::mat chain_tau0_e (nIter/2+1, n_env);
  // arma::vec chain_alpha (nIter/2+1);
  
  arma::mat chain_m ((nIter-burnin)/thin, l-1); //colnames(chain.m) = paste0("m.", 1:v);
  arma::cube chain_e ((nIter-burnin)/thin, l-1, n_env);
  arma::cube chain_b ((nIter-burnin)/thin, arma::size(B)(1), q);
  arma::mat chain_lambda1 ((nIter-burnin)/thin, q);
  arma::mat chain_tau2 ((nIter-burnin)/thin, q);
  arma::vec chain_rho ((nIter-burnin)/thin);
  arma::vec chain_se2 ((nIter-burnin)/thin);
  arma::vec chain_pi ((nIter-burnin)/thin);
  arma::mat chain_g ((nIter-burnin)/thin, q);
  arma::vec chain_tau0 ((nIter-burnin)/thin);
  arma::mat chain_tau0_e ((nIter-burnin)/thin, n_env);
  arma::vec chain_alpha ((nIter-burnin)/thin);
  arma::cube chain_xi ((nIter-burnin)/thin, arma::size(D)(0), q);
  
  // Rcout << "0.1 \n";
  // chain_m.row(0) = m.t();
  // for(j=0; j < n_env; j++){
  //   chain_e.slice(j).row(0) = e.col(j).t();
  // }
  // for(j=0; j<=q-1; j++){
  //   chain_b.slice(j).row(0) = b.col(j).t();
  // }
  // chain_tau2.row(0) = tau2.t();
  // chain_rho(0) = rho;
  // chain_se2(0) = se2;
  // chain_g.row(0) = g.t();
  // chain_tau0(0) = tau0;
  // chain_tau0_e.row(0) = tau0_e.t();
  // chain_alpha(0) = alpha;
  
  // Rcout << "Declaration des variables internes \n" ;
  // Declaration des variables internes :
  IntegerVector ind = Range(0, q-1);
  int iter = 0, i, ii, perc=2, iter_tmp = 0;
  double rho_i, rho_star, tmp_i, tmp_star, log_Qnew, log_Qold, rate_se2, log_ratio, a0, a1, a2, var_alpha, mean_alpha;
  arma::vec tmp_m(T), tmp_a(T), mean_m(l-1), tmp_e(T), mean_e(l-1), mean_bj(l), tmp_mean_bj(T), tmp(T), tmp_mean_cj(T), diag_0_Gamma(T), diag_1_Gamma(T);
  arma::vec tmp2(T), s_tmp(l), X2 (q), sumX(q), sumY(T), b2(q), one = arma::ones(T), Bee (T);
  arma::mat Gamma(T, T),Gamma_i(T, T), Gamma_star(T, T), Gamma_inv(T, T), Gamma_i_inv(T, T), Sigma_m(l-1, l-1), Sigma_e(l-1, l-1), tmp_sig_bj(l, l), Sigma_bj(l, l);
  arma::mat V1 (l, l), V2 (l,l);
  arma::mat Ytild(T, n), BGinv(T, l), YX(T, q), XX(q, q);
  Rcpp::IntegerVector rand(q);
  arma::vec one2 = arma::ones(arma::size(D)(1));
  // Rcout << "calculus \n" ;
  X2 = sum(pow(X, 2)).t();
  sumX = sum(X).t();
  sumY = sum(Y).t();
  YX = Y.t()*X;
  XX = X.t()*X;
  
  Bee = Be.slice(0) * e.col(0);
  for(i=1; i<n_env; i++) Bee += Be.slice(i) * e.col(i);
  
  
  Rprintf("\nSampler Progress...\n| ");
  for(i=0; i<100; i+=20) Rprintf("%d%%      ",i);
  Rprintf("%d%% |\n|",i);
  
  while(iter< nIter){          // GIBBS <----------------------------------
    // Rcout << iter << "\n";
    diag_0_Gamma =  arma::ones(T, 1) * (1+rho*rho); diag_0_Gamma(0) = diag_0_Gamma(T-1) = 1;  //arma::repmat( 1+rho*rho, T, 1);
    diag_1_Gamma = -rho* arma::ones(T-1, 1); //arma::repmat( -rho, T-1, 1);
    Gamma_inv = (arma::diagmat(diag_1_Gamma, -1) + arma::diagmat(diag_0_Gamma, 0) + arma::diagmat(diag_1_Gamma, 1))/(1-rho*rho);
    BGinv = B.t() * Gamma_inv;
    
    // Rcout << " 1.0 \n";
    // 1.0_ update alpha:
    var_alpha = as<double>(wrap(n * one.t() * Gamma_inv * one / se2));
    tmp_a = sumY - n*Bt*m - n*Bee - B*b*sumX;
    mean_alpha = as<double>(wrap( one.t() * Gamma_inv/se2 * tmp_a / var_alpha ));
    alpha = R::rnorm(mean_alpha, sqrt(1/var_alpha));
    
    // Rcout << " 1.1 \n";
    // 1_ update m:
    Sigma_m = inv(K2 /tau0 + Bt.t() * Gamma_inv *Bt * n / se2);
    tmp_m = sumY - n*alpha - n * Bee - B*b*sumX;
    mean_m = Sigma_m * Bt.t() * Gamma_inv * tmp_m / se2;
    m = rmnorm(mean_m, Sigma_m);
    tau0 = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(m.t()*K2*m /2 + r)));
    // tau0 = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(m.t()*K2*m /se2/2 + r)));
    
    // Rcout << " 1.2 \n";
    // 1_ update e (environnement):
    for(j=0; j< n_env; j++){
      Sigma_e = inv(K2 /tau0_e(j) + Be.slice(j).t() * Gamma_inv * Be.slice(j) * n / se2);
      tmp_e = sumY - n*alpha -n*Bt*m - n*(Bee - Be.slice(j) * e.col(j)) - B*b*sumX;
      mean_e = Sigma_e * Be.slice(j).t() * Gamma_inv * tmp_e / se2;
      e.col(j) = rmnorm(mean_e, Sigma_e);
      tau0_e(j) = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(e.col(j).t()*K2*e.col(j) /2 + r)));
      // tau0_e = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(e.t()*K2*e /se2/2 + r)));
      Bee = Be.slice(0) * e.col(0);
      for(i=1; i<n_env; i++) Bee += Be.slice(i) * e.col(i);
    }
    
    
    // Rcout << " 2 \n";
    // 2 update gamma, beta et tau2
    //---------------- Dirac sur le Slab
    rand = sample_int(q, 0, q-1);
    for(i=0; i<q; i++){
      // Rcout << i << "  \n";
      
      j = rand(i);
      // Sigma_bj = inv(K / (se2 * tau2(j)) + as<double>(wrap(X2(j))) *BGinv * B /se2 );
      Sigma_bj = inv((D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2)) / (se2 * tau2(j)) + as<double>(wrap(X2(j))) *BGinv * B /se2 );
      tmp_mean_bj = YX.col(j) - sumX(j)*one*alpha - sumX(j)*Bt*m - sumX(j)*Bee - B*b* XX.col(j) + B*b.col(j)* X2(j);
      // Rcout << " 2.2 \n";
      
      // // R = as<double>(wrap( pi/(1-pi)* sqrt(det(Sigma_bj)) / pow(se2*tau2(j), l/2) * exp(tmp_mean_bj.t() * Gamma_inv/se2  * B * Sigma_bj * B.t() * Gamma_inv/se2 * tmp_mean_bj / 2)));
      // // if(isinf(R)) ratio = 1; else ratio = R/(1+R);
      // log_R = as<double>(wrap(  log( pi/(1-pi) * sqrt(det(Sigma_bj)) / pow(se2*tau2(j), l/2) ) + tmp_mean_bj.t() * Gamma_inv/se2  * B * Sigma_bj * B.t() * Gamma_inv/se2 * tmp_mean_bj / 2  ));
      // // log_R = as<double>(wrap(  log( pi/(1-pi)) + log(det(D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2)))/2 + log(det(Sigma_bj))/2 -  l/2* log(tau2(j)*se2)  + tmp_mean_bj.t() * Gamma_inv/se2  * B * Sigma_bj * B.t() * Gamma_inv/se2 * tmp_mean_bj / 2  ));
      // ratio = exp(log_R - sumlogs(1, log_R));
      // // Rcout << " 2.3 \n";
      
      // if(arma::randu()< ratio){
        g(j) = 1;
        mean_bj = 1/se2  * Sigma_bj * BGinv * tmp_mean_bj;
        // Rcout << " 2.4 \n";
        
        b.col(j) = rmnorm(mean_bj, Sigma_bj);
        if(var_gp == "global"){
          tau2(j) = 1/R::rgamma((double)(l)/2+ s, 1/as<double>(wrap(b.col(j).t()*(D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2))*b.col(j) /se2/2 + r)));
          // phi_tau2(j) = 1/R::rgamma(1, 1/as<double>(wrap(1+1/tau2(j))));
        }
        if(var_gp == "local"){
          for(ii=0; ii<arma::size(D)(0); ii++){
            xi(ii, j) = R::rgamma((double)(1)/2 + s, 1/as<double>(wrap(arma::pow(D.row(ii)*b.col(j), 2) /tau2(j)/se2/2 + r)));
            // phi_xi(ii, j) = 1/R::rgamma(1, 1/as<double>(wrap(1+1/xi(ii, j))));
          }
        }
      // }
      // else{
      //   g(j) = 0;
      //   // Rcout << " 2.5 \n";
      //   
      //   b.col(j) = arma::zeros(arma::size(B)(1));
      //   tau2(j) = 1; //1/R::rgamma(s, 1/r);
      //   xi.col(j) = arma::zeros(arma::size(D)(0));
      //   // for(ii=0; ii<arma::size(D)(0); ii++){
      //   //   xi(ii, j) = R::rgamma((double)(1)/2 + s, 1/as<double>(wrap(arma::pow(D.row(ii)*b.col(j), 2) /tau2(j)/se2/2 + r)));
      //   // }
      // }
    }
    
    // pi = R::rbeta(1+sum(g), 1+q-sum(g));
    
    for(i=0; i<n; i++) Ytild.col(i) = (Y.row(i).t() - alpha - Bt*m - Bee - B*b*X.row(i).t());
    a0 = as<double>(wrap(accu(pow( Ytild.rows(1, T-2), 2))));
    a1 = as<double>(wrap(accu(Ytild.rows(0, T-2) % Ytild.rows(1, T-1))));
    a2 = a0 + as<double>(wrap(accu(pow(Ytild.row(0), 2)))) + as<double>(wrap(accu(pow(Ytild.row(T-1), 2))));
    
    // Rcout << " 3 \n";
    // 3. update rho (Metropolis-Hastings)
    rho_i = rho;
    tmp_i = (rho_i*rho_i * a0 - 2 * rho_i * a1 + a2)/(1-rho_i*rho_i);
    for(ii=0; ii<k; ii++){
      rho_star = R::runif(mmax(-1,rho-rho_tuning), mmin(1,rho+rho_tuning));
      log_Qnew=R::dunif(rho_star, mmax(-1, rho_i-rho_tuning), mmin(1, rho_i+rho_tuning), 1); // proposal distribution Q
      log_Qold=R::dunif(rho_i, mmax(-1, rho_star-rho_tuning), mmin(1, rho_star+rho_tuning), 1);
      
      tmp_star = (rho_star*rho_star * a0 - 2 * rho_star * a1 + a2)/(1-rho_star*rho_star);
      log_ratio = -n/2 *(T-1)* log(1-rho_star*rho_star) - 1/(2*se2)*tmp_star + log_Qold + n/2 * (T-1) * log(1-rho_i*rho_i) + 1/(2*se2)*  tmp_i - log_Qnew;
      
      if(R::runif(0, 1)< exp(log_ratio)) {
        rho_i = rho_star;
        tmp_i = tmp_star;
      }
    }
    rho = rho_i;
    
    // Rcout << " 4 \n";
    // 4. update sigma^2:
    rate_se2 = tmp_i/2;
    for(j=0; j<q; j++) rate_se2 += as<double>(wrap(b.col(j).t()*(D.t()*arma::diagmat(xi.col(j))*D + epsilon*arma::diagmat(one2))*b.col(j)/(tau2[j]*2)));
    rate_se2 += bb;
    se2 = 1 / R::rgamma(a + (double)(n*T)/2 + (double)(l*sum(g))/2, 1/rate_se2 );
    
    
    //---------------------------------------------------------------------------------------------------
    // Rcout << iter << "\n";
    if(iter>=burnin & (iter - burnin) % thin == 0){
      chain_m.row(iter_tmp) = m.t();
      for(j=0; j< n_env; j++){
        chain_e.slice(j).row(iter_tmp) = e.col(j).t();
      }
      chain_g.row(iter_tmp) = g.t();
      for(j=0; j<=q-1; j++){
        chain_b.slice(j).row(iter_tmp) = b.col(j).t();
      }
      for(j=0; j<=q-1; j++){
        chain_xi.slice(j).row(iter_tmp) = xi.col(j).t();
      }
      chain_tau2.row(iter_tmp) = tau2.t();
      chain_tau0(iter_tmp) = tau0;
      chain_tau0_e.row(iter_tmp) = tau0_e.t();
      chain_alpha(iter_tmp) = alpha;
      chain_rho(iter_tmp) = rho;
      chain_se2(iter_tmp) = se2;
      chain_pi(iter_tmp) = pi;
      iter_tmp ++;
    }
    iter ++;
    
    if(iter/(double)(nIter) >= perc/100.0){
      Rprintf("*");
      perc +=2;
    }
    
  }
  Rprintf("| \n");
  List chain = List::create(Named("g") = chain_g, Named("m") = chain_m, 
                            Named("e") = chain_e, Named("b") = chain_b, 
                            Named("rho") = chain_rho, Named("pi") = chain_pi,
                            Named("se2") = chain_se2, Named("tau2") = chain_tau2, 
                            Named("tau0") = chain_tau0, Named("tau0_e") = chain_tau0_e, 
                            Named("alpha") = chain_alpha, Named("xi") = chain_xi);
  return(chain);
  
}





