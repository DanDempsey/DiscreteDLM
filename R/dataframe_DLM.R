#' Creating a Lagged Dataframe
#'
#' Transforms a dataframe by lagging dynamic variables (using b-splines), to prepare it for distributed lag modelling. This is essentially a convenience function for quickly creating a dataframe that can be used for distributed lag modelling. This function interpolates the lag window of the dynamic variables using b-splines to create the cross-basis matrix; using the crossbasis function directly gives the user more flexibility.
#'
#' @param X A dataframe.
#' @param lag Lag length. See the help file for the "crossbasis" function from the dlnm package.
#' @param df Degrees of freedom for the b-splines function. See the bs documentation.
#' @param dynamic_vars The names of the columns that should be lagged. If not supplied, the dataset is not altered in any way.
#' @param scale Logical indicating if the user wishes the column of the lagged columns should be scaled to have zero mean and unit variance.
#' @param ... Further arguments to be passed into the crossbasis function.
#' @return A dataframe where the listed dynamic variables have been appropriately lagged. If no dynamic variables are given, the input dataframe is returned unaltered.
#' @author Daniel Dempsey (<dempsed6@tcd.ie>)
#' @examples
#' X <- dplyr::select( dlnm::chicagoNMMAPS, c('cvd', 'dow', 'temp', 'dptp', 'o3') )
#' X <- na.omit( X )
#' X_lagged <- dataframe_DLM( X, lag = 40, df = 7, dynamic_vars =  c('temp', 'dptp', 'o3') )
#' head( X_lagged )
#' @import splines
#' @import dlnm
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export

dataframe_DLM <- function( X, lag, df, dynamic_vars = NULL, scale = TRUE, ... ) {

  if ( is.null(dynamic_vars) ) {
    return( X )
  }
  varnames <- colnames( X )
  static_vars <- setdiff( varnames, dynamic_vars )
  X_static <- select( X, static_vars )
  X_dynamic_raw <- select( X, dynamic_vars )
  X_dynamic <- lapply( X_dynamic_raw, crossbasis, lag = lag, arglag = list(fun = 'bs', df = df), ... ) %>%
    do.call( what = 'cbind' )
  if ( scale ) {
    X_dynamic <- apply( X_dynamic, 2, scale )
  }
  colnames( X_dynamic ) <- paste0( rep( colnames(X_dynamic_raw), each = df ), c('', as.character(2:df)) )

  cbind( X_static, X_dynamic ) %>% na.omit

}

