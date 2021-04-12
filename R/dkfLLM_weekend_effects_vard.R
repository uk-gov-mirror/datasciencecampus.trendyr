##' Runs the Diffuse Kalman Filter (DFK).
##' locally linear trend - time varying drift - Saturday and Sunday effects
##'
##' The state space form is
##' y(x)  =  Z(x)*alpha(x) + G*u(x)
##' alpha(x+1) = TT*alpha(x) + H*u(x)
##' the DKF is initialized with A0 and P0 (Q0=0).
##'
##' sigma.eps^2 is concentrated out
##'
##' @param vec A vector with values of hyperparameters
##' @param YY A vector with the values of the time series
##' @param DOW A vector with the day of the week, same length as YY
##'
##' @author Sonia Mazzi
##' @export
dkf.llm.we.vard <- function(vec, YY, DOW){
  sigma.eps <- 1
  sigma.gnu <- vec[1]
  sigma.eta <- vec[2]
  delta <- vec[3]
  sigma.chi <- vec[4]
  K_chi <- 1 + 1 + 5^2
  rho <- 1
  y <- YY
  ll <- length(y)
  aux <- c(1, 1, 0, 0, 0, delta, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
  TT <- matrix(aux, 4, 4, byrow = T)
  lTT <- nrow(TT)
  A <- cbind(diag(-1, lTT), rep(0, lTT))
  ncA0 <- ncol(A)
  QQ <- matrix(0, ncA0, ncA0)
  DD <- rep(0, ll)
  EE <- matrix(0, 1, (ncA0*ll))
  AA <- matrix(rep(0, ncA0 * lTT * (ll + 1)), lTT, (ncA0 * (ll + 1)))
  AA[,(1:ncA0)] <- A
  PP <- matrix(rep(0, (lTT * lTT * (ll + 1))), lTT, (lTT * (ll + 1)))
  KK <- matrix(rep(0, (lTT * ll)), lTT, ll)
  elem <- c((sigma.gnu^2), 0, 0, 0,
            0, (sigma.eta^2), 0, 0,
            0, 0, sigma.chi^2 * K_chi/(K_chi-1), - sigma.chi^2 * 1/K_chi,
            0, 0, - sigma.chi^2 * 1/K_chi, sigma.chi^2 * K_chi/(K_chi-1) )
  HHt <- matrix(elem, lTT, lTT, byrow = T)
  P <- matrix(0 ,lTT, lTT)
  PP[ ,(1 : lTT)] <- P
  ee <- rep(0, ll)
  mse.ee <- rep(0, ll)
  GGt <- sigma.eps^2
  ZZ <- matrix(0, ll, 4)
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
    ZZ[i,] <- Z
    aux2 <- matrix(c(rep(0, (ncA0 - 1)), y[i]), 1, ncA0)
    E <- aux2 - Z %*% A
    EE[,((i-1) * ncA0 + 1) : (i * ncA0)] <- E
    EEgam <- matrix(E[, (1 : lTT)], 1, lTT)
    D <- Z %*% P %*% t(Z) + GGt
    DD[i] <- D
    Dinv <- 1/D
    K <- TT %*% P %*% t(Z) %*% Dinv
    KK[,i] <- K
    A <- TT %*% A + K %*% E
    AA[, ((i * ncA0 + 1) : ((i + 1) * ncA0))] <- A
    L <- TT - K %*% Z
    P <- L %*% P %*% t(TT) + HHt
    PP[, ((i * lTT + 1) : ((i + 1) * lTT))] <- P
    if (i > 9){
      SS <- QQ[(1:lTT), (1:lTT)]
      SSinv <- solve(SS)
      ss <- QQ[(1:lTT), ((lTT+1) : ncA0)]
      gamma <- SSinv %*% ss
      aux <- matrix(c(-gamma, 1), ncA0, 1)
      ee[i] <- E %*% aux
      mse.ee[i] <- EEgam %*% SSinv %*% t(EEgam)
      }
    QQ <- QQ + t(E) %*% Dinv %*% E
    }
    Alast <- A
    Plast <- P
    SS <- QQ[(1:lTT), (1:lTT)]
    qqq <- as.matrix(QQ[((lTT+1) : ncA0), ((lTT+1) : ncA0)])
    ss <- QQ[(1:lTT), ((lTT+1) : ncA0)]
    Sinv <- solve(SS)
    gamma.est <- Sinv %*% ss
    sigma2.est <- (qqq - t(ss) %*% Sinv %*% ss)/ll
    sigma2.est <- as.numeric(sigma2.est)
    mse.gamma.est <- sigma2.est*Sinv
    mse.ee <- sigma2.est * (DD + mse.ee)
    ee.std <- ee/sqrt(mse.ee)
    sigmatilde2 <- (ll / (ll - ncA0 + 1)) * sigma2.est
    loglik <- (-0.5) * ((ll - ncA0 + 1) * (1 + log(sigmatilde2)) + sum(log(abs(DD))))
    list(gamma.est = gamma.est, mse.gamma = mse.gamma.est, loglik = loglik, ncA0 = ncA0, ll = ll, delta = delta,
         D = DD, E = EE, Z = ZZ, A = AA, P = PP, K = KK, TT = TT, y = y, ee = ee, mse.ee = mse.ee,
         ee.std = ee.std,
         Alast = Alast,
         Plast = Plast,
         sigma2 = sigma2.est,
         sigma.eps2 = sigma2.est * sigma.eps^2,
         sigma.gnu2 = sigma2.est * sigma.gnu^2,
         sigma.eta2 = sigma2.est * sigma.eta^2,
         sigma.chi2 = sigma2.est * sigma.chi^2,
         sigma.eps = sqrt(sigma2.est * sigma.eps^2),
         sigma.gnu = sqrt(sigma2.est * sigma.gnu^2),
         sigma.eta = sqrt(sigma2.est * sigma.eta^2),
         sigma.chi = sqrt(sigma2.est * sigma.chi^2))
    }
