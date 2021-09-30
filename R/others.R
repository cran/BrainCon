#' tool functions
#'
#' @param u numeric value.
#'
#' @return Intermediate results.
#'
#' @name tool
#' @keywords internal
NULL

#' @rdname tool
QS = function(u){
  if (u == 0) ker = 1
  else ker = 25 * ( sin(6 * pi * u / 5) / (6 * pi * u / 5) - cos(6 * pi * u / 5) ) / (12 * pi^2 * u^2)
  return(ker)
}
