#' Print Method for MCMC_DLM Object
#'
#' Print output for MCMC_DLM objects. Displays basic information of regression model and the first few rows of the posterior.
#'
#' @param x An MCMC_DLM object.
#' @param n Number of rows of the posterior to display.
#' @param ... Further arguments; currently unused.
#' @author Daniel Dempsey (<dempsed6@tcd.ie>)
#' @export
print.MCMC_DLM <- function( x, n = 6, ... ) {

  NB_flag <- inherits(x, 'NB_MCMC')
  if( NB_flag ) {
    cat( '--- Negative binomial regression fit ---\n' )
  }
  else {
    cat( paste0('Quantile binary regression fit, with quantile = ', x$quantile, '\n') )
  }

  cat( paste0('First ', n, ' rows of beta posterior sample:\n') )
  head( x$beta, n ) %>% print
  cat( paste0('\nFirst ', n, ' rows of gamma posterior sample:\n') )
  head( x$gamma, n ) %>% print
  if( NB_flag ) {
    cat( paste0('\nFirst ', n, ' values of xi posterior sample:\n') )
    print( x$xi[1:n] )
  }

}

#' Summary Method for MCMC_DLM Object
#'
#' Summary output for MCMC_DLM objects. Displays summary statistics of the posterior.
#'
#' @param object An MCMC_DLM object.
#' @param digits Number of decimal places to round to for reported statistics.
#' @param ... Further arguments; currently unused.
#' @author Daniel Dempsey (<dempsed6@tcd.ie>)
#' @export
summary.MCMC_DLM <- function( object, digits = 2, ... ) {

  NB_flag <- inherits( object, 'NB_MCMC' )

  r2 <- function( z ) { round(z, digits) }

  cat( paste0('----- Summary of ', ifelse(NB_flag, 'Negative Binomial', 'Quantile Binary'), ' Regression Results -----\n\n') )

  if( NB_flag ) {
    cat( paste0('xi: mean = ', r2(mean(object$xi)), ', sd = ', r2(sd(object$xi)), '\n\n') )
  }

  for ( i in 1:ncol(object$X) ) {
    varkeep <- object$gamma[, i]
    var_beta_val <- object$beta[varkeep, i]
    cat( paste0(colnames(object$X)[i], ':\t mean = ', r2(mean(var_beta_val)), ', \tsd = ', r2(sd(var_beta_val)), ', \tmean inclusion = ', r2(mean(varkeep)), '\n') )
  }

}

#' Plotting Method for MCMC_DLM Object
#'
#' Provides some visualisations of the posterior.
#'
#' @param x An MCMC_DLM object.
#' @param type One of "beta", "gamma" or "xi". Partial matching is supported. Each type represents a graph for a different parameter: beta gives a ridgeplot of the regression slopes, gamma plots the mean probability of inclusion, and xi (negative binomial regression only) displays a kernel density estimate of the negative binomial stopping parameter.
#' @param ... Further arguments; currently unused.
#' @author Daniel Dempsey (<dempsed6@tcd.ie>)
#' @import ggplot2
#' @import ggridges
#' @importFrom reshape2 melt
#' @importFrom dplyr select
#' @export
plot.MCMC_DLM <- function( x, type = 'beta', ... ) {

  plot_type <- pmatch( type, c('beta', 'gamma', 'xi') )
  if ( is.na(plot_type) ) {
    stop( "type must be one of: 'beta', 'gamma', or 'xi'.\n" )
  }

  if ( plot_type == 1 ) {

    varkeep_long <- as.vector( x$gamma )
    beta_long <- reshape2::melt( x$beta, value.name = 'Value' )
    beta_long$varkeep <- varkeep_long
    beta_long_filter <- beta_long[beta_long$varkeep, ]

    bet_plot <- ggplot( data = beta_long_filter, aes_string(x = "Value", y = "Var2") ) +
      geom_density_ridges( scale = 0.9, alpha = 0.7, fill = 'blue' ) + theme_ridges() +
      labs( title = "Regression Slopes Ridge Plot", x = "", y = "Variable" )
    print( bet_plot )
    return( bet_plot )

  }

  if ( plot_type == 2 ) {

    gamma_dat <- data.frame( y = apply(x$gamma, 2, mean), x = colnames(x$gamma) )
    gam_plot <- ggplot( gamma_dat, aes_string(x = "x", y = "y") ) +
      geom_segment( aes_string(x = "x", xend = "x", y = 0, yend = "y" ), color = "blue") +  # Draw vertical lines
      geom_point( color = "blue", size = 2 ) +  # Add points on top of the lines
      labs( title = "", x = "Predictor", y = "Probability of Inclusion" )
    print( gam_plot )
    return( gam_plot )

  }

  if ( plot_type == 3 ) {

    if ( inherits(x, 'QB_MCMC') ) {
      stop( "type = xi is only available for negative binomial regression.\n")
    }

    xi_dat <- data.frame( xi = x$xi )
    xi_plot <- ggplot( xi_dat, aes_string(x = "xi") ) +
      geom_density( fill = "blue", alpha = 0.7 ) +
      labs( title = "Posterior of Negative Binomial Stopping Parameter", x = "", y = "" )
    print( xi_plot )
    return( xi_plot )

  }

}

