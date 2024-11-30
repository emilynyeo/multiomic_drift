function (n = 96, p = 179, G = 179) 
{
  mes <- floor(runif(n) + 3)
  time <- NULL
  id <- NULL
  nb2 <- c(1:n)
  for (i in 1:n) {
    time <- c(time, seq(1, mes[i], by = 1))
    id <- c(id, rep(nb2[i], length(seq(1, mes[i], by = 1))))
  }
  bruit <- floor(0 * p)
  bruit <- bruit + (p - bruit)%%G
  nices <- NULL
  for (i in 1:G) {
    nices <- c(nices, rep(i, (p - bruit)/G))
  }
  comportements <- matrix(0, length(time), G)
  comportements[, 1] <- 2.44 + 0.04 * (time - ((time - 6)^2)/(time/3))
  comportements[, 2] <- 0.5 * time - 0.1 * (time - 5)^2
  comportements[, 3] <- 0.25 * time - 0.05 * (time - 6)^2
  comportements[, 4] <- cos((time - 1)/3)
  comportements[, 5] <- 0.1 * time + sin(0.6 * time + 1.3)
  comportements[, 6] <- -0.1 * time^2
  X <- matrix(0, length(time), p)
  for (i in 1:(p - bruit)) {
    X[, i] <- comportements[, nices[i]] + rnorm(length(time), 
                                                0, 0.2)
  }
  for (j in 1:n) {
    w <- which(id == j)
    X[w, 1:(p - bruit)] <- X[w, 1:(p - bruit)] + rnorm(1, 
                                                       0, 0.1)
  }
  for (i in (p - bruit):p) {
    X[, i] <- rnorm(length(time), 0, 3)
  }
  f <- 1.3 * X[, 1]^2 + 2 * sqrt(abs(X[, which(nices == 2)[1]]))
  sigma <- cbind(c(0.5, 0.6), c(0.6, 3))
  Btilde <- matrix(0, length(unique(id)), 2)
  for (i in 1:length(unique(id))) {
    Btilde[i, ] <- rmvnorm(1, mean = rep(0, 2), sigma = sigma)
  }
  Z <- as.matrix(cbind(rep(1, length(f)), 2 * runif(length(f))))
  effets <- NULL
  for (i in 1:length(unique(id))) {
    w <- which(id == unique(id)[i])
    effets <- c(effets, Z[w, , drop = FALSE] %*% Btilde[i, 
    ])
  }
  gam <- 0.8
  BM <- NULL
  m <- length(unique(id))
  for (i in 1:m) {
    w <- which(id == unique(id)[i])
    W <- rep(0, length(w))
    t <- time[w]
    for (j in 2:length(w)) {
      W[j] <- W[j - 1] + sqrt(gam * (t[j] - t[j - 1])) * 
        rnorm(1, 0, 1)
    }
    BM <- c(BM, W)
  }
  sigma2 <- 0.5
  Y <- f + effets + rnorm(length(f), 0, sigma2) + BM
  return(list(Y = Y, X = X, Z = Z, id = id, time = time))
}