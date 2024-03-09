#' @title Sample from von Mises - Fisher distribution.
#' 
#' @description \code{rvMF} returns random draws from von Mises - Fisher distribution.
#'
#' @details The parameter theta is such that \eqn{dim(theta)} is the sphere dimension, \eqn{|theta|} the intensity parameter and \eqn{\frac{theta}{|theta|}} the mean directional parameter.
#' 
#' @param size as the number of draws needed.
#' @param theta as the distribution parameter.
#' @return A matrix whose each row is a random draw from the distribution.
#' @examples
#' # Draw 1000 vectors from vM-F with parameter 1, (1,0)
#' rvMF(1000,c(1,0))
#' 
#' # Draw 10 vectors from vM-F with parameter sqrt(14), (2,1,3)
#' rvMF(10,c(2,1,3))
#' 
#' # Draw from the vMF distribution with mean direction proportional 
#' # to c(1, -1) and concentration parameter 3
#' rvMF(10, 3 * c(1, -1) / sqrt(2))
#'   
#' @keywords distribution 
#' @keywords directional statistics 
#' @keywords coordinates
#' @keywords simulations
#' @references  
#' Wood, A. T. (1994). Simulation of the von Mises Fisher distribution. \emph{Communications in statistics-simulation and computation}, 23(1), 157-164. \doi{10.1080/03610919408813161}.
#' @references 
#' Hornik, K., & Grun, B. (2014). \pkg{movMF}: An \R package for fitting mixtures of von Mises-Fisher distributions. \emph{Journal of Statistical Software}, 58(10), 1-31. \doi{10.18637/jss.v058.i10}.
#' @export

rvMF <- function(size, theta){
  return(cpprvMF(size,theta))
}