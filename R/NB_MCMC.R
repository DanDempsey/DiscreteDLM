#' Negative Binomial Regression via MCMC
#'
#' Fits a negative binomial regression model using Gibbs sampling and performs Bayesian variable selection.
#'
#' @param formula Formula object to set the symbolic description of the model to be fitted.
#' @param data Optional dataframe, list or environment containing the model variables.
#' @param nsamp The desired sample size from the posterior. Set to 5000 by default.
#' @param nburn The number of iterations of the MCMC to be discarded as burn-in. Set to 5000 by default.
#' @param thin Thinning factor of the MCMC chain after burn-in. Set to 1 by default (no values discarded after burn-in).
#' @param standardize Logical indicating if the data should be standardised prior to fitting the model to zero mean and unit variance. If there is no intercept, then only scaling is applied (with a warning). The posterior sample of beta is transformed to the original scale.
#' @param prior_beta_mu Mean of the Gaussian prior for the regression coefficients, beta. Either a vector of length equal to the number of predictors or a single numeric to represent a constant vector.
#' @param prior_beta_sigma Covariance matrix of the Gaussian prior for the regression coefficients, beta. Either a square matrix of dimension equal to the number of predictors or a single numeric to represent an isotropic covariance matrix.
#' @param prior_gamma_p Probability parameter of the Bernoulli prior for the predictor inclusion parameter, gamma. Either a vector of length equal to the number of predictors or a single numeric, to represent that all predictors have the same prior probability of inclusion. Note that the model must have one predictor included by default; because of this the first value must be equal to 1 if a vector is given.
#' @param prior_xi_shape Shape parameter of the Gamma distribution prior for the negative binomial stopping parameter, xi.
#' @param prior_xi_scale Scale parameter of the Gamma distribution prior for the negative binomial stopping parameter, xi.
#' @param init_beta Initial MCMC values for the beta parameters. Either a vector of length equal to the number of predictors or a single numeric representing the same starting value for each beta component.
#' @param init_gamma Initial MCMC values for the gamma parameters. Either a vector of length equal to the number of predictors or a single numeric representing the same starting value for each gamma component. Note that the model must have one predictor included by default; because of this the first value must be equal to 1 if a vector is given.
#' @param init_xi Initial MCMC value for xi.
#' @return A list containing the MCMC-derived posterior sample, as well as the data that were used.
#' @details This function fits a negative binomial regression model via MCMC. Latent variable representation allows for Gibbs sampling of the parameters; see Pillow and Scott (2012) and Zhou et. al. (2012). The algorithm also includes predictor inclusion uncertainty inference (inferred via a Metroplis step) adapted from Holmes and Held (2006). The parameters of interest for this model are the regression slopes (beta), the binary predictor inclusion indicator (gamma) and the negative binomial stopping parameter (xi). For further details, see Dempsey and Wyse (2024).
#'
#' As this is a Bayesian model, priors must be specified. Beta has a Gaussian prior, gamma has a Bernoulli prior, and xi has a Gamma distribution prior.
#' @references
#' Jonathan Pillow and James Scott. "Fully Bayesian inference for neural models with negative-binomial spiking." Advances in neural information processing systems 25 (2012).
#'
#' Mingyuan Zhou, Lingbo Li, David Dunson and Lawrence Carin. "Lognormal and gamma mixed negative binomial regression." Proceedings of the... International Conference on Machine Learning. International Conference on Machine Learning. Vol. 2012. NIH Public Access, 2012.
#'
#' Chris C. Holmes and Leonhard Held. "Bayesian auxiliary variable models for binary and multinomial regression." Bayesian Analysis 1(1) (2006): 145-168.
#'
#' Daniel Dempsey and Jason Wyse (2024). Bayesian Generalized Distributed Lag Regression with Variable Selection. arXiv: https://arxiv.org/abs/2403.03646
#' @author Daniel Dempsey (<daniel.dempsey0@gmail.com>)
#' @examples
#' X <- dplyr::select( dlnm::chicagoNMMAPS, c('cvd', 'dow', 'temp', 'dptp', 'o3') )
#' X <- na.omit( X )
#' arglag <- list( fun = 'bs', df = 4 )
#' DLM_dat <- dataframe_DLM( X, lag = 40, dynamic_vars =  c('temp', 'dptp', 'o3'), arglag = arglag )
#' myfit <- NB_MCMC( cvd ~ ., data = DLM_dat, nsamp = 50, nburn = 50 )
#' summary( myfit )
#' @importFrom BayesLogit rpg
#' @import utils
#' @useDynLib DiscreteDLM, .registration = TRUE
#' @export
NB_MCMC <- function( formula, data = NULL, nsamp = 1000, nburn = 1000, thin = 1,
                     standardize = TRUE, prior_beta_mu = 0, prior_beta_sigma = 100,
                     prior_gamma_p = 0.5, prior_xi_shape = 2, prior_xi_scale = 1/50,
                     init_beta = 0, init_gamma = FALSE, init_xi = 1 ) {

  ### Initialize
  cat( "Initializing MCMC algorithm...\n" )
  MCMC_length <- nsamp + nburn

  # Convenient alternative to base sample function; doesn't treat scalars differently to vectors
  sample2 <- function( x, size, replace = FALSE, prob = NULL ) {
    if (missing(size)) {
      size <- length(x)
    }
    x[sample.int(length(x), size, replace, prob)]
  }

  # Utility function for back/forward solving instead of explicit inversion using Cholesky
  backfor <- function( x, y ) {
    backsolve( x, forwardsolve(t(x), y) )
  }

  # Set dataset
  X_full <- model.matrix( formula, data )
  nvar <- ncol( X_full )

  # Standardize
  col_mean <- rep( 0, nvar )
  col_sd <- rep( 1, nvar )
  int_ind <- which( colnames(X_full) == '(Intercept)' )
  if ( standardize ) {
    if ( length(int_ind) == 0 ) {
      warning( 'No intercept; only scaling will be applied.' )
      col_sd <- apply( X_full, 2, sd )
      int_ind <- 1
    }
    else {
      col_mean[-int_ind] <- apply(X_full[, -int_ind], 2, mean )
      col_sd[-int_ind] <- apply( X_full[, -int_ind], 2, sd )
    }
  }
  X_full <- scale( X_full, col_mean, col_sd )

  # Prepare priors and associated statistics used in the algorithm
  len_check <- function( x, nm, gamma_correction = FALSE ) {
    err_msg <- paste0( nm, ' must be of length 1 or length equal to the number of predictors. There are ', nvar, ' predictors.\n' )
    lenx <- length( x )
    if ( lenx == 1 ) {
      res <- rep(x, nvar)
      if ( gamma_correction ) {
        res[1] <- 1
      }
      return( res )
    }
    if ( (lenx < nvar) | (lenx > nvar) ) {
      stop( err_msg )
    }
    x
  }
  single_len_check <- function( x, nm ) {
    if ( length(x) > 1 ) {
      warning( paste0(x, ' must be only a single number.\n') )
      return( x[1] )
    }
    x
  }
  prior_beta_mu <- len_check( prior_beta_mu, 'prior_beta_mu' )
  prior_gamma_p <- len_check( prior_gamma_p, 'prior_gamma_p', TRUE )
  init_beta <- len_check( init_beta, 'init_beta' )
  init_gamma <- len_check( init_gamma, 'init_gamma', TRUE )
  if ( !is.matrix(prior_beta_sigma) ) {
    if ( length(prior_beta_sigma) == 1 ) {
      prior_beta_sigma <- diag( prior_beta_sigma, nvar )
    }
    else {
      stop( paste0('prior_beta_sigma must be a numeric of length one, or a square matrix of dimension equal to the number of predictors. There are ', nvar, ' predictors.\n') )
    }
  }
  if ( !all(nvar %in% dim(prior_beta_sigma)) ) {
    stop( paste0('prior_beta_sigma matrix must be square of dimension equal to the number of predictors. There are ', nvar, ' predictors.\n') )
  }
  single_len_check( prior_xi_shape, 'prior_xi_shape' )
  single_len_check( prior_xi_scale, 'prior_xi_scale' )
  single_len_check( init_xi, 'init_xi' )

  V0i_full <- chol2inv( chol(prior_beta_sigma) )
  V0ib0_full <- V0i_full %*% prior_beta_mu

  if ( !init_gamma[int_ind] ) {
    init_gamma[int_ind] <- TRUE
    warning( 'The gamma value corresponding to the intercept (or the first value if there is no intercept) must be TRUE. This has been corrected.\n' )
  }

  # Parameter initialization
  betares <- matrix( 0, ncol = nvar, nrow = MCMC_length ) # Regression slopes
  betares[1, ] <- init_beta
  gammares <- matrix( FALSE, ncol = nvar, nrow = MCMC_length ) # Variable inclusion indicator (starting value set in the next block of code)
  xires <- numeric( MCMC_length ) # Negative-Binomial stopping parameter
  xires[1] <- xi <- init_xi
  colnames( betares ) <- colnames( gammares ) <- colnames( X_full )

  # Variable index to link the dynamic variable gammas and apply starting value
  var_seq <- 1:nvar
  var_index_unique <- var_seq[-int_ind] # Used when proposing new gamma values
  gammares[1, ][var_seq[init_gamma]] <- TRUE
  gam <- gammares[1, ]

  X <- X_full[, gam, drop = FALSE]
  Xb <- X %*% init_beta[gam]
  V0i <- V0i_full[gam, gam]
  V0i_U <- chol( V0i )
  V0ib0 <- V0ib0_full[gam]

  y <- model.response( model.frame(formula, data = data) )
  y_len <- length( y )
  y_max <- max( y ) + 1

  # Starting value for Poyla-Gamma latent parameter and other associated statistics
  omega <- rpg( y_len, y + xi, Xb )
  lambda <- (y - xi) / (2 * omega)

  XtO <- t( X * omega )
  Vi_U <- chol( V0i + XtO%*%X )
  B <- backfor( Vi_U, V0ib0 + XtO%*%lambda )

  # Initialise matrix of pmfs for psi inference
  R <- matrix( 0, nrow = y_max, ncol = y_max )
  R[1, 1] <- R[2, 2] <- 1

  ### Main Loop
  cat( 'Initialization complete. Running algorithm...\n' )
  pb <- txtProgressBar( min = 2, max = MCMC_length, style = 3 )
  for ( i in 2:MCMC_length ) {

    ### Update xi
    # Poisson latent parameter: psi
    invisible(capture.output( sum_psi <- .C('PSI_FUN', xi, as.integer(y), as.integer(y_len),
                                            as.integer(0:(y_max-1)), as.integer(y_max), R,
                                            as.integer(0), PACKAGE = 'DiscreteDLM' )[[7]] ))

    # Negative Binomial overdispersion parameter: xi
    xires[i] <- xi <- rgamma( 1, prior_xi_shape + sum_psi, prior_xi_scale + sum(log(1 + exp(Xb))) )

    ### Update gamma
    # Propose change to parameter inclusion set
    change_ind <- var_seq[ sample2( var_index_unique, 1 ) ]
    gam_star <- gam
    gam_star[change_ind] <- !gam_star[change_ind]

    # Construct log acceptance ratio for proposed move
    X_star <- X_full[, gam_star, drop = FALSE]
    V0i_star <- V0i_full[gam_star, gam_star]
    V0i_star_U <- chol( V0i_star )
    V0ib0_star <- V0ib0_full[gam_star]
    XtO_star <- t( X_star * omega )

    Vi_star_U <- chol( V0i_star + XtO_star%*%X_star )
    B_star <- backfor( Vi_star_U, V0ib0_star + XtO_star%*%lambda )

    ldet_V <- -sum( log(diag(Vi_U)) )
    ldet_V_star <- -sum( log(diag(Vi_star_U)) )
    ldet_V0i <- sum( log( diag(V0i_U) ) )
    ldet_V0i_star <- sum( log(diag(V0i_star_U)) )

    gam_lprior <- dbinom( gam, 1, prior_gamma_p, log = TRUE )
    gam_lprior_star <- dbinom( gam_star, 1, prior_gamma_p, log = TRUE )

    lkernel <- t(B) %*% t(Vi_U) %*% Vi_U %*% B/2
    lkernel_star <- t(B_star) %*% t(Vi_star_U) %*% Vi_star_U %*% B_star/2

    lnum <- sum( ldet_V_star, ldet_V0i_star, lkernel_star, gam_lprior_star )
    ldenom <- sum( ldet_V, ldet_V0i, lkernel, gam_lprior )

    # Accept proposed update to gamma with M-H acceptance probability
    if ( (lnum - ldenom) > log(runif(1)) ) { # accept
      gammares[i, ] <- gam <- gam_star
      X <- X_star
      V0i <- V0i_star
      V0i_U <- V0i_star_U
      V0ib0 <- V0ib0_star
      Xb <- X %*% betares[i-1, gam]
    } else { # reject
      gammares[i, ] <- gammares[i-1, ]
    }

    ### Update beta
    # Polya-Gamma latent parameter: omega
    omega <- rpg( y_len, y + xi, Xb )
    lambda <- ( y - xi ) / ( 2 * omega )

    # Regression slopes: beta
    XtO <- t( X * omega )
    Vi_U <- chol( V0i + XtO%*%X )
    B <- backfor( Vi_U, V0ib0 + XtO%*%lambda )
    V_U <- chol( chol2inv(Vi_U) )

    betares[i, gam] <- B + t(V_U) %*% rnorm(sum(gam))
    Xb <- X %*% betares[i, gam]

    # Update progress
    setTxtProgressBar( pb, i )

  }

  ### Filter Markov chain and return result
  cat( '\nAlgorithm complete.\n' )

  keep <- seq( nburn + 1, MCMC_length, thin )
  gamma_trunc <- gammares[keep, ]

  # Convert betares back to original scale
  betares <- t( t(betares)/col_sd )
  if( length(int_ind) > 0 ) {
    betares[, int_ind] <- betares[, int_ind] - rowSums(t(t(betares) * col_mean))
  }

  # Return result
  res <- list( beta = betares[keep, ], gamma = gammares[keep, ], xi = xires[keep],
               model_matrix = X_full, data = data )
  class( res ) <- c( 'MCMC_DLM', 'NB_MCMC' )
  res

}

