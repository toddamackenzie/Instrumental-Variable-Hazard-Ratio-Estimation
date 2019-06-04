
library(survival)
library(rootSolve)

########################################################
# This function delivers estimator of constant HR based on the Instrumental variable W
IVHR <- function(Time, Status, X, W) {
  n <- length(Time)
  ord <- order(-Time)
  Time <- Time[ord]
  Status <- Status[ord]
  X <- X[ord]
  W <- W[ord]
  Est.Equat <- function(beta) {
    HR <- exp(beta *X)
    S.0 <- cumsum(HR)
    S.W1 <- cumsum(W*HR)
    sum(Status * (W - S.W1/S.0))
  }
  out.solution <- multiroot(Est.Equat, start=0)
  beta.hat <- ifelse(out.solution$estim.precis < 0.00001, out.solution$root, NA)
  HR <- exp(beta.hat *X)
  S.0 <- cumsum(HR)
  S.W1 <- cumsum(W*HR)
  S.X1 <- cumsum(X*HR)
  S.W1X1 <- cumsum(W*X*HR)
  S.W2 <- cumsum(W*W*HR)
  Var.E.E <- sum(Status * (W - S.W1/S.0)^2)
  Var.E.E.2 <- sum(Status * (S.W2/S.0 - (S.W1/S.0)^2))
  Deriv <- -sum(Status * (S.W1X1/S.0 - (S.W1*S.X1)/S.0^2))
  SE <- sqrt(Var.E.E/Deriv^2)
  SE.2 <- sqrt(Var.E.E.2/Deriv^2)
  list(Est.log.HR=beta.hat, SE=SE, SE.2=SE.2)
}

#######################################################
# This function delivers estimator of log linear HR based on the Instrumental variable W
IVHR.TD.Lin <- function(Time, Status, X, W) {
  n <- length(Time)
  ord <- order(Time) # ascending
  Time <- Time[ord]
  Status <- Status[ord]
  X <- X[ord]
  W <- W[ord]
  Est.Equat <- function(beta) {
    IV.Score <- c(0,0)
    for (i in 1:n) {
      if (Status[i]) {
        HR <- exp((beta[1] + beta[2]* Time[i])*X)
        AtRisk <- i:n
        S.0 <- sum(HR[AtRisk])
        S.W1 <- sum((W*HR)[AtRisk])
        IV.Score[1] <- IV.Score[1] + (W[i] - S.W1/S.0)
        IV.Score[2] <- IV.Score[2] + (W[i] - S.W1/S.0) * Time[i]
      }
    }
    IV.Score
  }
  out.solution <- multiroot(Est.Equat, start=c(0,0))
  if (is.na(out.solution$estim.precis)) {
    beta.hat <- c(NA, NA)
    Cov <- matrix(nrow=2,ncol=2, rep(NA,4))
  }
  if (!is.na(out.solution$estim.precis)) {
    if (out.solution$estim.precis > 0.00001) {
      beta.hat <- c(NA, NA)
      Cov <- matrix(nrow=2,ncol=2, rep(NA,4))
    }
    if (out.solution$estim.precis <= 0.00001) {
      beta.hat <- out.solution$root
      Dev.XP <- Deriv <- matrix(nrow=2,ncol=2,c(0,0,0,0))
      for (i in 1:n) {
        if (Status[i]) {
          HR <- exp((beta.hat[1] + beta.hat[2]*Time[i]) * X)
          AtRisk <- i:n
          S.0 <- sum(HR[AtRisk])
          S.W1 <- sum((W*HR)[AtRisk])
          S.X1 <- sum((X*HR)[AtRisk])
          Dev <- W[i] - S.W1/S.0
          Dev.XP <- Dev.XP + Dev^2 * matrix(nrow=2,ncol=2, c(1, Time[i], Time[i], Time[i]^2))
          S.W1X1 <- sum((W*X*HR)[AtRisk])
          Deriv <- Deriv - (S.W1X1/S.0 - (S.W1*S.X1)/S.0^2) * matrix(nrow=2,ncol=2, c(1, Time[i], Time[i], Time[i]^2))
        }
      }
      Inv.Deriv <- solve(Deriv)
      Cov <- Inv.Deriv %*% Dev.XP %*% Inv.Deriv
    }
  }    
  list(Est.log.HR=beta.hat, Cov=Cov)
}

#######################################################
# This function delivers estimator of log linear HR based on maximum partial likelihood
Cox.TD.Lin <- function(Time, Status, X) {
  n <- length(Time)
  ord <- order(Time) # ascending
  Time <- Time[ord]
  Status <- Status[ord]
  X <- X[ord]
  Est.Equat <- function(beta) {
    Score <- c(0,0)
    for (i in 1:n) {
      if (Status[i]) {
        HR <- exp((beta[1] + beta[2]* Time[i])*X)
        AtRisk <- i:n
        S.0 <- sum(HR[AtRisk])
        S.X1 <- sum((X*HR)[AtRisk])
        Score[1] <- Score[1] + (X[i] - S.X1/S.0)
        Score[2] <- Score[2] + (X[i] - S.X1/S.0) * Time[i]
      }
    }
    Score
  }
  out.solution <- multiroot(Est.Equat, start=c(0,0))
  if (is.na(out.solution$estim.precis)) beta.hat <- c(NA, NA)
  if (!is.na(out.solution$estim.precis)) {
    if (out.solution$estim.precis <= 0.00001) beta.hat <- out.solution$root
    if (out.solution$estim.precis > 0.00001) beta.hat <- c(NA, NA)
    }
  Dev.XP <- Deriv <- matrix(nrow=2,ncol=2,c(0,0,0,0))
  for (i in 1:n) {
    if (Status[i]) {
      HR <- exp((beta.hat[1] + beta.hat[2]*Time[i]) * X)
      AtRisk <- i:n
      S.0 <- sum(HR[AtRisk])
      S.X1 <- sum((X*HR)[AtRisk])
      Dev <- X[i] - S.X1/S.0
      Dev.XP <- Dev.XP + Dev^2 * matrix(nrow=2,ncol=2, c(1, Time[i], Time[i], Time[i]^2))
      S.X2 <- sum((X^2*HR)[AtRisk])
      Deriv <- Deriv - (S.X2/S.0 - (S.X1*S.X1)/S.0^2) * matrix(nrow=2,ncol=2, c(1, Time[i], Time[i], Time[i]^2))
    }
  }
  Inv.Deriv <- solve(Deriv)
  Cov <- Inv.Deriv %*% Dev.XP %*% Inv.Deriv
  list(Est.log.HR=beta.hat, Cov=Cov)
}

