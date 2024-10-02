#' Creating a Lagged Dataframe
#'
#' Transforms a dataframe by smoothing over lags of dynamic variables to prepare it for distributed lag modelling.
#'
#' @param X A dataframe.
#' @param lag Lag length. See the help file for the "crossbasis" function from the dlnm package.
#' @param dynamic_vars The names of the columns that should be lagged. If not supplied, the dataset is not altered in any way (except possibly for scaling) and a warning is supplied.
#' @param arglag A list that is passed into onebasis for generating a basis matrix.
#' @param scale Logical indicating if the user wishes the columns of the resulting dataframe should be scaled to have zero mean and unit variance.
#' @param ... Further arguments to be passed into the crossbasis function.
#' @return A dataframe where the listed dynamic variables have been appropriately lagged. If no dynamic variables are given, the input dataframe is returned unaltered (besides standardisation if scale = TRUE) with a warning.
#' @author Daniel Dempsey (<dempsed6@tcd.ie>)
#' @examples
#' X <- dplyr::select( dlnm::chicagoNMMAPS, c('cvd', 'dow', 'temp', 'dptp', 'o3') )
#' X <- na.omit( X )
#' arglag <- list( fun = 'bs', df = 4 )
#' DLM_dat <- dataframe_DLM( X, lag = 40, dynamic_vars =  c('temp', 'dptp', 'o3'), arglag = arglag )
#' head( DLM_dat )
#' @import splines
#' @import dlnm
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export

dataframe_DLM <- function( X, lag, dynamic_vars = NULL, arglag = list(fun = 'bs'), scale = TRUE, ... ) {

  if ( is.null(dynamic_vars) ) {
    if ( scale ) {
      X <- apply( X, 2, scale )
    }
    warning( 'No dynamic variables listed.' )
    return( X )
  }

  static_vars <- setdiff( colnames(X), dynamic_vars )
  X_static <- select( X, static_vars )
  X_dynamic_raw <- select( X, dynamic_vars )
  X_dynamic <- lapply( X_dynamic_raw, crossbasis, lag = lag, arglag = arglag, ... ) %>%
    do.call( what = 'cbind' )

  if ( scale ) {
    X_dynamic <- apply( X_dynamic, 2, scale )
  }

  l <- ncol( X_dynamic ) / length( dynamic_vars )
  colnames( X_dynamic ) <- paste0( rep( colnames(X_dynamic_raw), each = l ), paste0('.l', 1:l) )

  cbind( X_static, X_dynamic ) %>% na.omit

}

