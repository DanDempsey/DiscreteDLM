QB_MCMC <- function( formula, data = NULL, quantile = 0.5, nsamp = 1000,
                     nburn = 1000, thin = 1, standardize = TRUE, prior_beta_mu = 0,
                     prior_beta_sigma = 100, prior_gamma_p = 0.5, init_beta = 0,
                     init_gamma = FALSE, verbose = TRUE ) {

  ### Initialize
  if ( verbose ) {
    cat( "Initializing MCMC algorithm...\n" )
  }
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
    backsolve( x, backsolve(x, y, transpose = TRUE) )
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

  V0i_full <- chol2inv( chol(prior_beta_sigma) )
  V0ib0_full <- V0i_full %*% prior_beta_mu

  if ( !init_gamma[int_ind] ) {
    init_gamma[int_ind] <- TRUE
    warning( 'The gamma value corresponding to the intercept (or the first value if there is no intercept) must be TRUE. This has been corrected.\n' )
  }

  # Parameter initialization
  nvar <- ncol( X_full )
  betares <- matrix( 0, ncol = nvar, nrow = MCMC_length ) # Regression slopes
  betares[1, ] <- init_beta
  gammares <- matrix( FALSE, ncol = nvar, nrow = MCMC_length ) # Variable inclusion indicator (starting value set in the next block of code)
  colnames( betares ) <- colnames( gammares ) <- colnames( X_full )

  # Variable index to link the dynamic variable gammas and apply starting value
  var_seq <- 1:nvar
  var_index_unique <- var_seq[-int_ind] # Used when proposing new gamma values
  gammares[1, ][as.logical(init_gamma)] <- TRUE
  gam <- gammares[1, ]

  X <- X_full[, gam, drop = FALSE]
  Xb <- X %*% init_beta[gam]
  V0i <- V0i_full[gam, gam]
  V0i_U <- chol( V0i )
  V0ib0 <- V0ib0_full[gam]

  y <- model.response( model.frame(formula, data = data) )
  y_len <- length( y )
  y_max <- max( y ) + 1
  upper_tail <- ifelse( y, TRUE, FALSE )

  # Starting values for latent variables
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  phi <- 2 / ( quantile * (1 - quantile) )
  delta <- 2 + ( ( psi^2 ) / phi )

  ystar <- rTALD( n = y_len, upper_tail = upper_tail, mu = Xb, sigma = 1, p = quantile )
  chi <- ( ystar - Xb )^2 / phi
  nu <- 1/rinvgauss( y_len, mean = sqrt(delta/chi), shape = delta )
  lambda <- ystar - ( psi * nu )
  omega <- 1 / (phi * nu)

  XtO <- t(X * omega)
  Vi_U <- chol( V0i + XtO%*%X )
  B <- backfor( Vi_U, V0ib0 + XtO%*%lambda )

  ### Main Loop
  if ( verbose ) {
    cat( 'Initialization complete. Running algorithm...\n' )
    pb <- txtProgressBar( min = 2, max = MCMC_length, style = 3 )
  }
  for ( i in 2:MCMC_length ) {

    ### Update gamma
    # Propose change to parameter inclusion set
    change_ind <- var_seq[ sample2(var_index_unique, 1) ]
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

    lkernel <- 0.5 * sum( (Vi_U %*% B)^2 )
    lkernel_star <- 0.5 * sum( (Vi_star_U %*% B_star)^2 )

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

    ### Update latent parameters
    # ystar
    ystar <- rTALD( n = y_len, upper_tail = upper_tail, mu = Xb, sigma = 1, p = quantile )

    # nu
    chi <- ( ystar - Xb )^2 / phi
    nu <- 1/rinvgauss( y_len, mean = sqrt(delta/chi), shape = delta )

    ### Update beta
    lambda <- ystar - ( psi * nu )
    omega <- 1 / ( phi * nu )

    XtO <- t( X * omega )
    Vi_U <- chol( V0i + XtO%*%X )
    B <- backfor( Vi_U, V0ib0 + XtO%*%lambda )

    betares[i, gam] <- B + backsolve( Vi_U, rnorm(sum(gam)) )
    Xb <- X %*% betares[i, gam]

    # Update progress
    if ( verbose ) {
      setTxtProgressBar( pb, i )
    }

  }

  ### Filter Markov chain and return result
  if ( verbose ) {
    cat( '\nAlgorithm complete.\n' )
  }

  keep <- seq( nburn + 1, MCMC_length, thin )
  gamma_trunc <- gammares[keep, ]

  # Convert betares back to original scale
  betares <- t( t(betares)/col_sd )
  if( length(int_ind) > 0 ) {
    betares[, int_ind] <- betares[, int_ind] - rowSums(t(t(betares) * col_mean))
  }

  # Return result
  res <- list( beta = betares[keep, ], gamma = gammares[keep, ], quantile = quantile,
               model_matrix = X_full, data = data )
  class( res ) <- c( 'MCMC_DLM', 'QB_MCMC' )
  res

}

