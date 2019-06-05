
library(survival)
library(MASS)
library(rootSolve)

########################################################
# This is the function which delivers estimates based on the Instrumental variable W
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
  #  to.min <- function(beta) {Est.Equat(beta)^2}
  #  out.solution <- nlminb(start=0, to.min)
  #  beta.hat <- ifelse(out.solution$convergence == 0, out.solution$par, NA)
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



###########################################################

Inverse.SDF.0 <- function(u) {-log(u)} # inverse survival function for unit exponential - no treatment
Inverse.SDF.1.loglin <- function(u, Coefs) { # inverse survival function for treated based on 3 piece stepwise time dependent HR
  h <- -log(u)
  if (Coefs[2] == 0) time <- h/Coefs[1]
  if (Coefs[2] != 0) time <- log(1+Coefs[2]*h/exp(Coefs[1]))/Coefs[2]
  time
}



# If the second coefficient is negative, hazard may be so close to zero that there is a mass at infinity
n <- 5000
W <- rnorm(n)
X <- 1/(1 + exp(-W)) < runif(n)
Y0 <- -log(runif(n))
Y1 <- Inverse.SDF.1.loglin(runif(n), Coefs=c(-1,+0.5))
Y <- ifelse(X, Y1, Y0)
Cox.TD.Lin(Y, rep(1,n), X)$Est.log.HR
IVHR.TD.Lin(Y, rep(1,n), X, W)$Est.log.HR



##########################################################
# Set rho at 0.75, sample size at 2000, and cens freq at 50% 
N <- 1000
Cens.Freq <- 0.5
rho <- 0.50
##########################################################

Intercept.Grid <- seq(from=-1,to=+1, length=11)
Slope.Grid <- seq(from=-0.5, to=0.5, length=11)

Results.Table <- NULL
Out.Var.Names <- c("Intercept", "Slope", "rho", "log.odds.ratio.X.vs.W", "log.odds.ratio.X.vs.U", "cox.intercept", "cox.slope", "iv.intercept", "iv.slope")
File_Name <- paste("C:/Users/Todd MacKenzie/Dropbox/ivhr time dependent/Sim Results/LogLinearTD Out-", Sys.Date(), ".txt", sep="")
while(1) {
  # Parameter setting
	log.odds.ratio.X.vs.U <- log(10)*(-1+2*runif(1))
	log.odds.ratio.X.vs.W <- log(5) + log(5)*runif(1)  # odds at least 5 up to 25
	Coefs <- c( sample(Intercept.Grid, 1), sample(Slope.Grid, 1) )
  Settings <- c(Coefs, rho, log.odds.ratio.X.vs.W, log.odds.ratio.X.vs.U)
	# Instrument (continuous), omitted covariate (continuous) and exposure (binary)
	W <- rnorm(N)
	U <- rnorm(N) # omitted covariate or latent factor affecting survival
	X <- runif(N) < 1/(1+exp(-log.odds.ratio.X.vs.U * U - log.odds.ratio.X.vs.W * W))
	# Generation of the Potential Outcomes using Gaussian Copula
	Norm.0 <- rho * U + sqrt(1-rho^2) * rnorm(N)
	Norm.1 <- rho * U + sqrt(1-rho^2) * rnorm(N)
	Survival.0 <- Inverse.SDF.0(pnorm(Norm.0))
	Survival.1 <- Inverse.SDF.1.loglin(pnorm(Norm.1), Coefs) #NaN results if it is a sub-distribution due to haz -> 0 too fast
  Survival.1 <- ifelse(!is.na(Survival.1), Survival.1, 2 * max(Survival.0)) 
  # Observed Potential Outcome
	Survival <- ifelse(X, Survival.1, Survival.0)
  # Censoring times; censoring fixed at 50%
	Censoring.Time <- runif(N)
	Censoring.Time <- Censoring.Time * quantile(Survival/Censoring.Time, Cens.Freq)
	# Now the observed time-to-event and status
	Time.to.Event <- pmin(Survival, Censoring.Time)
	Status <- ifelse(Survival < Censoring.Time, 1, 0)
	# Estimation by Cox
	o.cox <- Cox.TD.Lin(Time.to.Event, Status, X)
  o.iv <- IVHR.TD.Lin(Time.to.Event, Status, X, W)
  Estimates <- c(o.cox$Est.log.HR, o.iv$Est.log.HR)
  NewRow <- c(Settings, Estimates)
  Results.Table <- rbind(Results.Table, NewRow)
  if (dim(Results.Table)[1] %% 10000 == 0) write.table(Results.Table, File_Name, row.names=FALSE, col.names=Out.Var.Names)
  if (dim(Results.Table)[1] %% 100 == 0) cat(dim(Results.Table)[1], "\t")
}

dimnames(Results.Table)[[2]] <- c("Intercept", "Slope", "rho", "log.odds.ratio.X.vs.W", "log.odds.ratio.X.vs.U", "cox.intercept", "cox.slope", "iv.intercept", "iv.slope")

HR <- Results.Table[, c("Intercept", "Slope")]
Cox.Est <- Results.Table[, c("cox.intercept", "cox.slope")]
IV.Est <- Results.Table[, c("iv.intercept", "iv.slope")]

#png("C:/Users/Todd MacKenzie/Dropbox/ivhr time dependent/figures/sim_loglin.png")

par(mfrow=c(2,3))
rng.intercept <- (4/3) * range(Intercept.Grid)
rng.slope <- (4/3) * range(Slope.Grid)
FU.Name <- c("Intercept", "Slope")
for (i in 1:2) {
  if (i == 1) {
    Grid <- Intercept.Grid
    rng <- rng.intercept
  }
  if (i == 2) {
    Grid <- Slope.Grid
    rng <- rng.slope
  }
  keep <- Results.Table[, "log.odds.ratio.X.vs.U"] > log(5)  
  plot(Grid, tapply(IV.Est[keep,i], HR[keep,i], median, na.rm=TRUE), pch=16, xlab="True Coefficient", ylab="Estimate", main=paste(FU.Name[i],"\n+ Confounding"), xlim=rng, ylim=rng)
  points(Grid, tapply(Cox.Est[keep,i], HR[keep,i], median, na.rm=TRUE), col=2, pch=16)
  abline(a=0,b=1, lty=3)
  keep <- Results.Table[, "log.odds.ratio.X.vs.U"] > -log(2) & Results.Table[, "log.odds.ratio.X.vs.U"] < log(2)  
  plot(Grid*exp(0.05), tapply(IV.Est[keep,i], HR[keep,i], median, na.rm=TRUE)*exp(0.05), pch=16, xlab="True Coefficient", ylab="Estimate", main="Little or No Confounding", xlim=rng, ylim=rng)
  points(Grid*exp(-0.05), tapply(Cox.Est[keep,i], HR[keep,i], median, na.rm=TRUE)*exp(-0.05), col=2, pch=16)
  abline(a=0,b=1, lty=3)
  keep <- Results.Table[, "log.odds.ratio.X.vs.U"] < -log(5)  
  plot(Grid, tapply(IV.Est[keep,i], HR[keep,i], median, na.rm=TRUE), pch=16, xlab="True Coefficient", ylab="Estimate", main="- Confounding", xlim=rng, ylim=rng)
  points(Grid, tapply(Cox.Est[keep,i], HR[keep,i], median, na.rm=TRUE), col=2, pch=16)
  abline(a=0,b=1, lty=3)
}

dev.off()


