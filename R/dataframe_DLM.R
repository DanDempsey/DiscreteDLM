#' Creating a Lagged Dataframe
#'
#' Transforms a dataframe by smoothing over lags of dynamic variables to prepare it for distributed lag modelling.
#'
#' @param X A dataframe.
#' @param lag Lag length. See the help file for the "crossbasis" function from the dlnm package.
#' @param dynamic_vars The names of the columns that should be lagged. If not supplied, the dataset is not altered in any way (except possibly for scaling) and a warning is supplied.
#' @param arglag A list that is passed into onebasis for generating a basis matrix.
#' @param ... Further arguments to be passed into the crossbasis function.
#' @return A dataframe where the listed dynamic variables have been appropriately lagged. If no dynamic variables are given, the input dataframe is returned unaltered with a warning.
#' @author Daniel Dempsey (<daniel.dempsey0@gmail.com>)
#' @examples
#' X <- dplyr::select( dlnm::chicagoNMMAPS, c('cvd', 'dow', 'temp', 'dptp', 'o3') )
#' X <- na.omit( X )
#' arglag <- list( fun = 'bs', df = 4 )
#' DLM_dat <- dataframe_DLM( X, lag = 40, dynamic_vars =  c('temp', 'dptp', 'o3'), arglag = arglag )
#' @importFrom dplyr select %>%
#' @import splines
#' @import dlnm
#' @export
dataframe_DLM <- function( X, lag, dynamic_vars = NULL, arglag = list(fun = 'bs'), ... ) {

  if ( is.null(dynamic_vars) ) {
    warning( 'No dynamic variables listed. Returning given data.' )
    return( X )
  }

  if ( missing(lag) ) {
    stop( "The 'lag' argument must be supplied." )
  }

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

#' @export
as.data.frame.dataframe_DLM <- function( x, ... ) {
  as.data.frame( x$data, ... )
}

#' @export
print.dataframe_DLM <- function( x, ... ) {
  cat( 'Data:\n')
  print( x$data, ... )
  cat( paste0('\nLag: ', x$lag, '\n\n') )
  cat( 'Arglag:\n' )
  print( x$arglag, ... )
  cat( '\n' )
}

#' @export
summary.dataframe_DLM <- function( object, ... ) {
  print( summary(as.data.frame(object), ...) )
  cat( paste0('\nLag: ', object$lag, '\n\n') )
}

#' @export
plot.dataframe_DLM <- function( x, ... ) {
  plot( as.data.frame(x), ... )
}

