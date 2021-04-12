##' This is the function that computes the likelihood for LLM-SSM with daily effects with a time varying drift (second component of alpha)
##' Daily effects are weekdays (5), saturdays (1) and Sundays (1)
##'
##' For daily data
##' @param vec A vector with values of hyperparameters
##' @param YY A vector with the values of the time series. Make sure to define externally
##' @param DOW A vector with the day of the week, same length as YY
##'
##' @author Sonia Mazzi
##' @export
lik.llm.we.vard <- function(vec, y, DOW){
  sigma.eps <- 1
  sigma.gnu <- vec[1]
  sigma.eta <- vec[2]
  delta <- vec[3]
  sigma.chi <- vec[4]
  K_chi <- 1 + 1 + 5^2
  ll <- length(y)
  aux <- c(1, 1, 0, 0, 0, delta, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
  TT <- matrix(aux, 4, 4, byrow = T)
  lTT <- nrow(TT)
  A <- cbind(diag(-1, lTT), rep(0, lTT))
  ncA0 <- ncol(A)
  QQ <- matrix(0, ncA0, ncA0)
  DD <- rep(0, ll)
  elem <- c((sigma.gnu^2), 0, 0, 0,
            0, (sigma.eta^2), 0, 0,
            0, 0, sigma.chi^2 * K_chi/(K_chi-1), - sigma.chi^2 * 1/K_chi,
            0, 0, - sigma.chi^2 * 1/K_chi, sigma.chi^2 * K_chi/(K_chi-1) )
  HHt <- matrix(elem, lTT, lTT, byrow = T)
  GGt <- sigma.eps^2
  P0 <- matrix(0, lTT, lTT)
  P <- P0
  #
  Z_sat <- c(1, 0, 1, 0)
  Z_sun <- c(1, 0, 0, 1)
  Z_other <- c(1, 0, 0, 0)
  Z_choices <- matrix(c(Z_sat, Z_sun, Z_other), byrow = TRUE, 3, 4)
  for (i in 1:ll){
    Z_row <- ifelse(DOW[i] == "Saturday", 1 ,
                ifelse(DOW[i] == "Sunday", 2, 3
                       )
                )
    Z <- matrix(Z_choices[Z_row,], 1, 4)
    aux2 <- matrix(c(rep(0, (ncA0-1)), y[i]), 1, ncA0)
    E <- aux2 - Z %*% A
    D <- Z %*% P %*% t(Z) + GGt
    DD[i] <- D
    Dinv <- 1/D
    K <- TT %*% P %*% t(Z) %*% Dinv
    A <- TT %*% A + K %*% E
    L <- TT - K %*% Z
    P <- L %*% P %*% t(TT) + HHt
    QQ <- QQ + t(E) %*% Dinv %*% E
    }
  #
  SS <- QQ[(1:lTT), (1:lTT)]
  qqq <- as.matrix(QQ[((lTT + 1) : ncA0), ((lTT + 1) : ncA0)])
  ss <- QQ[(1 : lTT), ((lTT + 1) : ncA0)]
  Sinv <- solve(SS)
  gamma.est <- Sinv %*% ss
  sigma2.est <- (qqq - t(ss) %*% Sinv %*% ss) / ll
  sigma2.est <- as.numeric(sigma2.est)
  sigmatilde2 <- (ll / (ll - ncA0 + 1)) * sigma2.est
  loglik <- (-0.5) * ((ll - ncA0 + 1) * (1 + log(sigmatilde2)) + sum(log(abs(DD))))
  -loglik
  }
