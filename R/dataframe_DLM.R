dataframe_DLM <- function( X, lag, dynamic_vars = NULL, arglag = list(fun = 'bs'), ... ) {

  if ( is.null(dynamic_vars) ) {
    warning( 'No dynamic variables listed. Returning given data.' )
    return( X )
  }

  if ( missing(lag) ) {
    stop( "The 'lag' argument must be supplied." )
  }

  X <- as.data.frame( X )
  static_vars <- setdiff( colnames(X), dynamic_vars )
  X_static <- select( X, static_vars )
  X_dynamic_raw <- select( X, dynamic_vars )
  X_dynamic <- lapply( X_dynamic_raw, crossbasis, lag = lag, arglag = arglag, ... ) %>%
    do.call( what = 'cbind' )

  l <- ncol( X_dynamic ) / length( dynamic_vars )
  dynamic_names <- paste0( rep( dynamic_vars, each = l ), paste0('.l', 1:l) )
  colnames( X_dynamic ) <- dynamic_names
  dynamic_names_list <- split( dynamic_names, rep(1:ncol(X_dynamic_raw), each = l) )
  names( dynamic_names_list ) <- dynamic_vars

  res <- list( data = cbind(X_static, X_dynamic) %>% na.omit,
               dynamic_names = dynamic_names_list, lag = lag, arglag = arglag )
  class( res ) <- 'dataframe_DLM'
  res

}

as.data.frame.dataframe_DLM <- function( x, ... ) {
  as.data.frame( x$data, ... )
}

print.dataframe_DLM <- function( x, ... ) {
  cat( 'Data:\n')
  print( x$data, ... )
  cat( paste0('\nLag: ', x$lag, '\n\n') )
  cat( 'Arglag:\n' )
  print( x$arglag, ... )
  cat( '\n' )
}

summary.dataframe_DLM <- function( object, ... ) {
  print( summary(as.data.frame(object), ...) )
  cat( paste0('\nLag: ', object$lag, '\n\n') )
}

plot.dataframe_DLM <- function( x, ... ) {
  plot( as.data.frame(x), ... )
}

