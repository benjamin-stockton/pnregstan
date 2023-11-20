#' Convert Angles to Unit Vectors
#'
#' @param theta a vector of angles
#'
#' @return U a unit vector or matrix of unit vectors in each row
#' @export
#'
#' @examples
#' theta <- seq(0, 2*pi, length.out = 10)
#' (U <- angle_to_unit_vec(theta))
angle_to_unit_vec <- function(theta) {
  if (length(theta) == 1 && is.numeric(theta) && theta <= 2*pi) {
    U <- c(cos(theta), sin(theta))
    return(U)
  }
  else if (length(theta) > 1 && is.numeric(theta)) {
    U <- cbind(cos(theta), sin(theta))
  } 
  return(U)
}