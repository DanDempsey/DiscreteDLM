#' Binary Quantile Regression via MCMC
#'
#' Fits a quantile regression model to a binary response dataset using Gibbs sampling and performs Bayesian variable selection.
#'
#' @param formula Formula object to set the symbolic description of the model to be fitted.
#' @param data Optional dataframe, list or environment containing the model variables.
#' @param quantile The chosen quantile for regression.
#' @param nsamp The desired sample size from the posterior. Set to 5000 by default.
#' @param nburn The number of iterations of the MCMC to be discarded as burn-in. Set to 5000 by default.
#' @param thin Thinning factor of the MCMC chain after burn-in. Set to 1 by default (no values discarded after burn-in).
#' @param prior_beta_mu Mean of the Gaussian prior for the regression coefficients, beta. Either a vector of length equal to the number of predictors or a single numeric to represent a constant vector.
#' @param prior_beta_sigma Covariance matrix of the Gaussian prior for the regression coefficients, beta. Either a square matrix of dimension equal to the number of predictors or a single numeric to represent an isotropic covariance matrix.
#' @param prior_gamma_p Probability parameter of the Bernoulli prior for the predictor inclusion parameter, gamma. Either a vector of length equal to the number of predictors or a single numeric, to represent that all predictors have the same prior probability of inclusion. Note that the model must have one predictor included by default; because of this the first value must be equal to 1 if a vector is given.
#' @param init_beta Initial MCMC values for the beta parameters. Either a vector of length equal to the number of predictors or a single numeric representing the same starting value for each beta component.
#' @param init_gamma Initial MCMC values for the gamma parameters. Either a vector of length equal to the number of predictors or a single numeric representing the same starting value for each gamma component. Note that the model must have one predictor included by default; because of this the first value must be equal to 1 if a vector is given.
#' @return A list containing the MCMC-derived posterior sample, as well as the data that were used.
#' @details This function fits a quantile binary regression model via MCMC. Latent variable representation allows for Gibbs sampling of the parameters; see Benoit (2017). The algorithm also includes predictor inclusion uncertainty inference (inferred via a Metroplis step) adapted from Holmes and Held (2006). The parameters of interest for this model are the regression slopes (beta) and the binary predictor inclusion indicator (gamma).
#'
#' As this is a Bayesian model, priors must be specified. Beta has a Gaussian prior and gamma has a Bernoulli prior.
#' @references
#' Dries F. Benoit and Dirk Van den Poel. "bayesQR: A Bayesian approach to quantile regression." Journal of Statistical Software 76 (2017): 1-32.
#'
#' Chris C. Holmes and Leonhard Held. "Bayesian auxiliary variable models for binary and multinomial regression." (2006): 145-168.
#' @author Daniel Dempsey (<dempsed6@tcd.ie>)
#' @examples
#' set.seed( 100 )
#' Pima_ex <- MASS::Pima.tr
#' Pima_ex$type <- ifelse( MASS::Pima.tr$type == 'Yes', 1, 0 )
#' pima_mcmc <- QB_MCMC( type ~ ., data = Pima_ex )
#' @import statmod
#' @import mvtnorm
#' @import utils
#' @importFrom MASS area
#' @export

### Load libraries
#library( statmod ) # For Inverse Gaussian
#library( mvtnorm ) # For mu
#source( 'Code/Main_Software/ALD.R' ) # ALD functions

### Main wrapper function for performing MCMC-based inference for the Negative Binomial DLM
QB_MCMC <- function( formula, data = NULL, quantile = 0.5, nsamp = 5000,
                     nburn = 5000, thin = 1, prior_beta_mu = 0, prior_beta_sigma = 100,
                     prior_gamma_p = 0.5, init_beta = 0, init_gamma = FALSE ) {

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

  # Set dataset and groups
  X_full <- model.matrix( formula, data )
  nvar <- ncol( X_full )
  groups <- 1:nvar

  # Prepare priors and associated statistics used in the algorithm
  len_check <- function( x, nm, gamma_correction = FALSE ) {
    err_msg <- paste0( nm, ' must be of length 1 or length equal to the number of predictors. There are ', nvar, ' predictors.' )
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
  prior_beta_mu <- len_check( prior_beta_mu, 'prior_beta_mu' )
  prior_gamma_p <- len_check( prior_gamma_p, 'prior_gamma_p', TRUE )
  init_beta <- len_check( init_beta, 'init_beta' )
  init_gamma <- len_check( init_gamma, 'init_gamma', TRUE )
  if ( !is.matrix(prior_beta_sigma) ) {
    if ( length(prior_beta_sigma) == 1 ) {
      prior_beta_sigma <- diag( prior_beta_sigma, nvar )
    }
    else {
      stop( paste0('prior_beta_sigma must be a numeric of length one, or a square matrix of dimension equal to the number of predictors. There are ', nvar, ' predictors.') )
    }
  }
  if ( !all(nvar %in% dim(prior_beta_sigma)) ) {
    stop( paste0('prior_beta_sigma matrix must be square of dimension equal to the number of predictors. There are ', nvar, ' predictors.') )
  }

  V0i_full <- chol2inv( chol(prior_beta_sigma) )
  V0ib0_full <- V0i_full %*% prior_beta_mu

  if ( !init_gamma[1] ) {
    init_gamma[1] <- TRUE
    warning( 'The first element of the gamma starting value must be TRUE. This has been corrected.' )
  }

  # Parameter initialization
  nvar <- ncol( X_full )
  betares <- matrix( 0, ncol = nvar, nrow = MCMC_length ) # Regression slopes
  betares[1, ] <- init_beta
  gammares <- matrix( FALSE, ncol = nvar, nrow = MCMC_length ) # Variable inclusion indicator (starting value set in the next block of code)
  colnames( betares ) <- colnames( gammares ) <- colnames( X_full )

  # Variable index to link the dynamic variable gammas and apply starting value
  var_split <- split( 1:nvar, groups )
  var_index_unique <- unique( groups )[-1] # Used when proposing new gamma values
  gammares[1, ][unlist(var_split[init_gamma])] <- TRUE
  gam <- gammares[1, ]

  X <- X_full[, gam, drop = FALSE]
  Xb <- X %*% init_beta[gam]
  V0i <- V0i_full[gam, gam]
  V0ib0 <- V0ib0_full[gam]

  y <- data[[all.vars(formula)[1]]]
  y_len <- length( y )
  y_max <- max( y ) + 1
  rtrunc <- ifelse( y, TRUE, FALSE )

  # Starting values for latent variables
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  phi <- 2 / ( quantile * (1 - quantile) )
  delta <- 2 + ( ( psi^2 ) / phi )

  ystar <- rTALD( n = y_len, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
  chi <- (ystar - Xb)^2 / phi
  nu <- 1/rinvgauss( y_len, mean = sqrt( delta/chi ), shape = delta )
  lambda <- ystar - (psi * nu)

  omega <- 1 / (phi * nu)
  XtO <- t(X * omega)

  V <- chol2inv( chol(V0i + XtO%*%X) )
  B <- V%*%( V0ib0 + XtO%*%lambda )

  ### Main Loop
  cat( 'Initialization complete. Running algorithm...\n' )
  pb <- txtProgressBar( min = 2, max = MCMC_length, style = 3 )
  for ( i in 2:MCMC_length ) {

    ### Update gamma
    # Propose change to parameter inclusion set
    change_ind <- var_split[[ sample2( var_index_unique, 1 ) ]]
    gam_star <- gam
    gam_star[change_ind] <- !gam_star[change_ind]

    # Construct log acceptance ratio for proposed move
    X_star <- X_full[, gam_star, drop = FALSE]
    V0i_star <- V0i_full[gam_star, gam_star]
    V0ib0_star <- V0ib0_full[gam_star]
    XtO_star <- t( X_star * omega )

    V_star <- chol2inv( chol(V0i_star + XtO_star%*%X_star) )
    B_star <- V_star%*%( V0ib0_star + XtO_star%*%lambda )

    ldet_V <- sum( log(diag(chol(V))) )
    ldet_V_star <- sum( log(diag(chol(V_star))) )
    ldet_V0i <- sum( log( diag(chol(V0i)) ) )
    ldet_V0i_star <- sum( log(diag(chol(V0i_star))) )

    gam_lprior <- dbinom( gam, 1, prior_gamma_p, log = TRUE )
    gam_lprior_star <- dbinom( gam_star, 1, prior_gamma_p, log = TRUE )

    lkernel <- crossprod(B, chol2inv(chol(V)))%*%B/2
    lkernel_star <- crossprod(B_star, chol2inv(chol(V_star)))%*%B_star/2

    lnum <- sum( ldet_V_star, ldet_V0i_star, lkernel_star, gam_lprior_star )
    ldenom <- sum( ldet_V, ldet_V0i, lkernel, gam_lprior )

    # Accept proposed update to gamma with M-H acceptance probability
    if ( (lnum - ldenom) > log(runif(1)) ) { # accept
      gammares[i, ] <- gam <- gam_star
      X <- X_star
      V0i <- V0i_star
      V0ib0 <- V0ib0_star
      Xb <- X %*% betares[i-1, gam]
    } else { # reject
      gammares[i, ] <- gammares[i-1, ]
    }

    ### Update latent parameters
    # ystar
    ystar <- rTALD( n = y_len, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )

    # nu
    chi <- ( ystar - Xb )^2 / phi
    nu <- 1/rinvgauss( y_len, mean = sqrt(delta/chi), shape = delta )

    ### Update beta
    lambda <- ystar - ( psi * nu )
    omega <- 1 / ( phi * nu )
    XtO <- t( X * omega )
    V <- chol2inv( chol(V0i + XtO%*%X) )
    B <- V%*%( V0ib0 + XtO%*%lambda )
    betares[i, gam] <- B + t(chol(V)) %*% rnorm(sum(gam))
    Xb <- X %*% betares[i, gam]

    # Update progress
    setTxtProgressBar( pb, i )

  }

  ### Filter Markov chain and return result
  cat( 'Algorithm complete. Returning result.\n' )

  keep <- seq( nburn + 1, MCMC_length, thin )
  col_inds <- sapply(var_split, '[', 1)
  gamma_trunc <- gammares[keep, col_inds]

  list( beta = betares[keep, ], gamma = gammares[keep, ], quantile = quantile,
        X = X_full, y = y )

}

