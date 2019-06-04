remove(list=ls())
library(survival)
library(rootSolve)
library(MASS)

IVHR.multi <- function(Time, Status, X, IV, Z=NULL) {
  # X is the exposure
  # IV is the instrumental variable
  # Z is an optional matrix of covariates
  n <- length(Time)
  ord <- order(-Time)
  Time <- Time[ord]
  Status <- Status[ord]
  X <- as.matrix(as.matrix(X)[ord,])
  IV <- as.matrix(as.matrix(IV)[ord,])
  if (dim(X)[2] != dim(IV)[2]) stop("X and IV must have equal number of columns.")
  if (dim(X)[1] != dim(IV)[1]) stop("X and IV must have equal number of rows.")
  if (is.null(Z)) Z <- matrix(nrow=n, ncol=0)
  Z <- as.matrix(Z)[ord,]
  XX <- cbind(X, Z)
  W <- cbind(IV, Z)
  S.W1 <- matrix(nrow=n, ncol=dim(W)[2])
  Est.Equat <- function(beta) {
    HR <- exp(XX %*% beta)
    S.0 <- cumsum(HR)
    for (i.col in 1:dim(W)[2])
      S.W1[,i.col] <- cumsum(W[,i.col]*HR)
    colSums(Status * (W - S.W1/S.0))
    }
#  to.min <- function(beta) {Est.Equat(beta)^2}
#  out.solution <- nlminb(start=0, to.min)
#  beta.hat <- ifelse(out.solution$convergence == 0, out.solution$par, NA)
  out.solution <- multiroot(Est.Equat, start= rep(0, dim(W)[2]))
  no.root.found <- is.na(out.solution$estim.precis) | out.solution$estim.precis >= 0.00001
  if(no.root.found) {
    beta.hat <- rep(NA, dim(W)[2])
    Sandwich <- matrix(NA, nrow=dim(W)[2], ncol=dim(W)[2])
    } 
  if(!no.root.found) {
    beta.hat <- out.solution$root
    # Variance
    HR <- exp(XX %*% beta.hat)
    S.0 <- cumsum(HR)
    W1 <- S.X1 <- matrix(nrow=n, ncol=dim(W)[2])
    for (i.col in 1:dim(W)[2]) {
      S.W1[,i.col] <- cumsum(W[,i.col]*HR)
      S.X1[,i.col] <- cumsum(XX[,i.col]*HR)
      }
    S.W1X1 <- array(dim=c(n, dim(W)[2], dim(W)[2]))
    Var.E.E <-  Deriv <- matrix(nrow=dim(W)[2], ncol=dim(W)[2])
   for (i.col in 1:dim(W)[2]) {
      for (j.col in 1:dim(W)[2]) {
        S.W1X1[,i.col, j.col] <- cumsum(W[,i.col]*XX[,j.col]*HR)
        Var.E.E[i.col, j.col] <- sum(Status * (W[, i.col] - S.W1[, i.col]/S.0) * (W[, j.col] - S.W1[, j.col]/S.0))
        Deriv[i.col, j.col] <- -sum(Status * (S.W1X1[, i.col, j.col]/S.0 - (S.W1[, i.col]*S.X1[, j.col])/S.0^2))
        }
      }
    Inv.Deriv <- ginv(Deriv)
    Sandwich <- t(Inv.Deriv) %*% Var.E.E %*% Inv.Deriv
    }
  list(Est.log.HR=beta.hat, SE=diag(Sandwich)^0.5, Sandwich=Sandwich)
  }


# An illustration using simulated data
# Data is simulated assuming time-to-event depends on covariates, including omitted covariate, according to multivariable Cox
# Unlike the marginal Cox  model or additive hazards effect of omitted covariate it is designed for 

log.odds.IV.on.X <- log(5)
log.odds.Omit.on.X <- log(3)
log.HR.X.on.Time <- 0.5 # that which we wish to estimate
log.HR.Omitted.on.Time <- log(0.5)

# Choosing a big sample size here, but no reason to think small sample sizes will not work too
n <- 100000
Censoring.rate <- 0.75
IV <- rnorm(n)
Omitted <- rnorm(n)
Cov.1 <- rnorm(n) # Confounder, affects X and Time,
Cov.2 <- rnorm(n) # Affects X
Cov.3 <- rnorm(n) # Affects Time only
X <- 1/(1+exp(-(log.odds.IV.on.X*IV + log.odds.Omit.on.X*Omitted + 0.6*Cov.1 + 0.3 * Cov.2 + 0.0*Cov.3))) < runif(n)
Y <- -log(runif(n)) / exp( log.HR.X.on.Time*X + log.HR.Omitted.on.Time*Omitted - 0.5*Cov.1 - 0.0*Cov.2 + 0.5*Cov.3) 
C <- runif(n)
c.rescale <- sort(Y/C)[n*(1-Censoring.rate)]
C <- C * c.rescale
T <- pmin(Y, C)
St  <- ifelse(Y <= C, 1, 0)

Cov <- cbind(Cov.1, Cov.2, Cov.3)

# If we know the omitted covariate
round(o <- coxph(Surv(T, St) ~ X + Omitted + Cov)$coef,2)

# Omitted covariate causes bias
round(o <- coxph(Surv(T, St) ~ X + Cov)$coef, 2)

# Instrumental variable aproach based on the method described in the paper in Health Services Outcomes and Research Methodology
o.IV <- IVHR.multi(T, St, X, IV, Cov)
round(cbind(o.IV$Est.log.HR, o.IV$SE), 2)


