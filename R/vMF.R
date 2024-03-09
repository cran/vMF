#' @title Sample from von Mises - Fisher distribution
#'
#' @description vMF samples from von Mises-Fisher distribution and performs some operations. 
#' Unlike the \href{https://cran.r-project.org/package=movMF}{movMF} package, vMF does not consider mixtures of von Mises-Fisher distribution.
#' vFM particularly focuses on sampling from the distribution and performs it very quickly. This is useful to carry out fast simulations in directional statistics.
#' vMF also computes the density and normalization constant of the von Mises-Fisher distribution. 
#' 
#' 
#' @author 
#' Aristide Houndetoungan <\email{ariel92and@@gmail.com}>
#'   
#' @references  
#' Wood, A. T. (1994). Simulation of the von Mises Fisher distribution. \emph{Communications in statistics-simulation and computation}, 23(1), 157-164. \doi{10.1080/03610919408813161}.
#' @references 
#' Hornik, K., & Grun, B. (2014). \pkg{movMF}: An \R package for fitting mixtures of von Mises-Fisher distributions. \emph{Journal of Statistical Software}, 58(10), 1-31. \doi{10.18637/jss.v058.i10}.
#' @useDynLib vMF, .registration = TRUE
"_PACKAGE"
NULL