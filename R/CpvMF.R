#' @title  Normalization constant of von Mises - Fisher distribution.
#' 
#' @description \code{CpvMF} returns the normalization constant of von Mises - Fisher density.
#'
#' @details  The probability density function of the von Mises - Fisher distribution is defined by :
#' \deqn{f(z|theta) = C_p(|theta|)\exp{(z theta)}}
#' \eqn{|theta|} is the intensity parameter and \eqn{\frac{theta}{|theta|}} the mean directional parameter. The normalization constant \eqn{C_p()} depends 
#' on the Bessel function of the first kind. See more details \href{https://en.wikipedia.org/wiki/Von_Mises-Fisher_distribution}{here}.
#' 
#' @param p as sphere dimension.
#' @param k as the intensity parameter.
#' @return the normalization constant.
#' @examples
#' 
#' CpvMF(2,3.1)
#' 
#' @keywords distribution 
#' @keywords directional statistics 
#' @keywords coordinates
#' @keywords simulations
#' @seealso 
#' \code{\link{rvMF}} and \code{\link{dvMF}}
#' 
#' @references  
#' Wood, A. T. (1994). Simulation of the von Mises Fisher distribution. \emph{Communications in statistics-simulation and computation}, 23(1), 157-164. \url{https://www.tandfonline.com/doi/abs/10.1080/03610919408813161}.
#' @references 
#' Hornik, K., & Grun, B. (2014). \pkg{movMF}: An \R package for fitting mixtures of von Mises-Fisher distributions. \emph{Journal of Statistical Software}, 58(10), 1-31. \url{https://epub.wu.ac.at/4893/}.
#' @importFrom Rcpp sourceCpp
#' @export


CpvMF <- function(p, k){
  return(cppCpvMF(p,k))
}