## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----comp1, echo = TRUE, eval = TRUE------------------------------------------
library(vMF)
library(movMF)

n      <- 5 # Number of draws
set.seed(123)
xvMF   <- rvMF(n,c(1,0,0))
set.seed(123)
xmovMF <- rmovMF(n,c(1,0,0))
all.equal(c(xvMF), c(xmovMF))
xvMF
xmovMF

## ----comp2, echo = TRUE, eval = TRUE------------------------------------------
n      <- 30
set.seed(123)
ddpcr::quiet(runif(n))
xvMF   <- rvMF(n,c(1,0,0))
set.seed(123)
xmovMF <- rmovMF(n,c(1,0,0))
all.equal(c(xvMF), c(xmovMF))
xvMF[1:5,]
xmovMF[1:5,]

## ----ex1, echo = TRUE, eval = TRUE--------------------------------------------
library(rbenchmark) 

fcompare <- function(n) {
  benchmark("vMF" = rvMF(n,c(1,0,0)), "movMF" = rmovMF(1,c(1,0,0)))
}

fcompare(1)
fcompare(10)
fcompare(100)

## ----ex11a, echo = FALSE, eval = TRUE-----------------------------------------
#load data to save time during the building
load("out.Rdata")

## ----ex11b, echo = TRUE, eval = FALSE-----------------------------------------
#  out  <- unlist(lapply(1:200, function(x) fcompare(x)$elapsed[1]/fcompare(x)$elapsed[2]))

## ----ex11c, echo = TRUE, eval = TRUE, fig.height = 4, fig.align = "center"----
library(ggplot2)
ggplot(data = data.frame(n = 1:200, time = out), aes(x = n, y = time)) +
  geom_point(col = "blue") + geom_hline(yintercept = 1, col = 2)

## ----ex2a, echo = TRUE, eval = FALSE------------------------------------------
#  set.seed(123)
#  P                  <- 4
#  initial            <- rmovMF(1, rep(0, P))
#  # Fonction based on vMF to simulate theta
#  SamplevMF          <- function(n) {
#    output           <- matrix(0, n + 1, P)
#    output[1, ]      <- initial
#    for (i in 1:n) {
#      output[i + 1,] <- rvMF(1, output[i,])
#    }
#    return(output)
#  }
#  
#  # Fonction based on movMF to simulate theta
#  SamplemovMF        <-function(n){
#    output           <- matrix(0, n + 1, P)
#    output[1, ]      <- initial
#    for (i in 1:n) {
#      output[i + 1,] <- rmovMF(1, output[i,])
#    }
#    return(output)
#  }
#  benchmark("vMF" = SamplevMF(1000), "movMF" = SamplemovMF(1000))

## ----ex2b, echo = FALSE, eval = TRUE------------------------------------------
print(outbench)

