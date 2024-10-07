#' Visualises Results from MCMC_DLM Fits
#'
#' Provides some ggplot visualisations of the posterior from MCMC_DLM objects.
#'
#' @param x An MCMC_DLM object.
#' @param type One of "beta", "gamma" or "xi". Partial matching is supported. Each type represents a graph for a different parameter: beta gives a ridgeplot of the regression slopes, gamma plots the mean probability of inclusion, and xi (negative binomial regression only) displays a kernel density estimate of the negative binomial stopping parameter.
#' @param include_intercept Logical indicating if the intercept should be included in the graphs.
#' @param print_output Logical indicating whether or not the ggplot objects should be printed to the screen.
#' @author Daniel Dempsey (<dempsed6@tcd.ie>)
#' @return A ggplot object. Will also display the ggplot visual if print_output = TRUE.
#' @import ggplot2
#' @import ggridges
#' @importFrom reshape2 melt
#' @importFrom dplyr select
#' @importFrom dlnm onebasis
#' @export
vis_output <- function( x, type = 'beta', include_intercept = FALSE, print_output = TRUE ) {

  if ( !inherits(x, 'MCMC_DLM') ) {
    stop( "x must be an MCMC_DLM object." )
  }

  plot_type <- pmatch( type, c('beta', 'gamma', 'xi', 'lags') )
  if ( is.na(plot_type) ) {
    stop( "type must be one of: 'beta', 'gamma', 'xi', or 'lags'.\n" )
  }

  if ( plot_type == 1 ) {

    varkeep_long <- as.vector( x$gamma )
    beta_long <- melt( x$beta, value.name = 'Value' )
    beta_long$varkeep <- varkeep_long
    beta_long_filter <- beta_long[beta_long$varkeep, ]

    if ( !include_intercept ) {
      int_inds <- beta_long_filter$Var2 == '(Intercept)'
      beta_long_filter <- beta_long_filter[!int_inds, ]
    }

    bet_plot <- ggplot( data = beta_long_filter, aes_string(x = "Value", y = "Var2") ) +
      geom_density_ridges( scale = 0.9, alpha = 0.7, fill = 'blue' ) + theme_ridges() +
      labs( title = "Regression Slopes Ridge Plot", x = "", y = "Variable" )

    if ( print_output ) {
      print( bet_plot )
    }

    return( bet_plot )

  }

  if ( plot_type == 2 ) {

    gamma_dat <- data.frame( y = apply(x$gamma, 2, mean), x = colnames(x$gamma) )
    if ( !include_intercept ) {
      int_inds <- gamma_dat$x == '(Intercept)'
      gamma_dat <- gamma_dat[!int_inds, ]
    }

    gam_plot <- ggplot( gamma_dat, aes_string(x = "x", y = "y") ) +
      geom_segment( aes_string(x = "x", xend = "x", y = 0, yend = "y" ), color = "blue") +  # Draw vertical lines
      geom_point( color = "blue", size = 2 ) +  # Add points on top of the lines
      labs( title = "", x = "Predictor", y = "Probability of Inclusion" )

    if ( print_output ) {
      print( gam_plot )
    }

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

    if ( print_output ) {
      print( xi_plot )
    }

    return( xi_plot )

  }

  if ( plot_type == 4 ) {

    arglag_full <- c( list(x = 0:x$data$lag), x$data$arglag )
    lag_splines <- do.call( 'onebasis', arglag_full )
    beta_list <- lapply( x$data$dynamic_names, function(z) { x = x$beta[, z] } )

    lag_plot_fun <- function( x, nm ) {

      fitted_weights_raw <- lag_splines %*% t( x )
      fitted_weights_raw_mean <- apply( fitted_weights_raw, 1, mean )
      mean_effect <- sum( fitted_weights_raw_mean )
      fitted_weights_mean <- fitted_weights_raw_mean / mean_effect

      graph_dat <- data.frame( x = seq_along(fitted_weights_mean) - 1, y = fitted_weights_mean )

      ggplot( graph_dat, aes_string(x = 'x', y = 'y') ) + geom_line() +
        labs( title = paste0(nm, ' lag response'), x = "Lag" )

    }

    lag_plots <- Map( lag_plot_fun, x = beta_list, nm = names(beta_list) )

    if ( print_output ) {
      for ( i in 1:length(lag_plots) ) {
        print( lag_plots[[i]] )
      }
    }

    return( lag_plots )

  }

}

