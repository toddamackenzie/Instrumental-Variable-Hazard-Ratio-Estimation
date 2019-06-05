
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
  if (out.solution$estim.precis < 0.00001) beta.hat <- out.solution$root
  if (out.solution$estim.precis >= 0.00001) beta.hat <- rep(NA,2) 
  Dev.XP <- Deriv <- matrix(nrow=2,ncol=2,c(0,0,0,0))
  for (i in 1:n) {
    if (Status[i]) {
      HR <- exp((beta.hat[1] + beta.hat[2]*Time[i]) * X)
      AtRisk <- i:n
      S.0 <- sum(HR[AtRisk])
      S.W1 <- sum((W*HR)[AtRisk])
      Dev <- W[i] - S.W1/S.0
      Dev.XP <- Dev.XP + Dev^2 * matrix(nrow=2,ncol=2, c(1, Time[i], Time[i], Time[i]^2))
      S.W1X1 <- sum((W*X*HR)[AtRisk])
      Deriv <- Deriv - (S.W1X1/S.0 - (S.W1*S.X1)/S.0^2) * matrix(nrow=2,ncol=2, c(1, Time[i], Time[i], Time[i]^2))
    }
  }
  Inv.Deriv <- solve(Deriv)
  Cov <- Inv.Deriv %*% Dev.XP %*% Inv.Deriv
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
  if (out.solution$estim.precis < 0.00001) beta.hat <- out.solution$root
  else beta.hat <- rep(NA,2)
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
Inverse.SDF.1.3step <- function(u, HR, TimeCuts) { # inverse survival function for treated based on 3 piece stepwise time dependent HR
  h <- -log(u)
  Lambda.Cuts <- cumsum(c(HR[1]*TimeCuts[1], HR[2]*(TimeCuts[2]-TimeCuts[1])))
  t1 <- h/HR[1]
  t2 <- TimeCuts[1] + (h - Lambda.Cuts[1])/HR[2]
  t3 <- TimeCuts[2] + (h - Lambda.Cuts[2])/HR[3]
  ifelse(h < Lambda.Cuts[1], t1, ifelse(h < Lambda.Cuts[2], t2, t3))
  }
Inverse.SDF.1.loglin <- function(u, Coefs) { # inverse survival function for treated based on 3 piece stepwise time dependent HR
  h <- -log(u)
  log(1+Coefs[2]*h/exp(Coefs[1]))/Coefs[2]
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

TimeCuts <- -log(1 - c(0.10, 0.25))
n.periods <- 1 + length(TimeCuts)
HR.Grid <- c(1/3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.5, 1.7, 2, 2.5, 3)

##########################################################
# Set sample size at 2000, and cens freq at 50% for 
N <- 1000
Cens.Freq <- 0.5
rho <- 0.50
log.odds.ratio.X.vs.U.grid <- log(c(1/5, 1, 5))
##########################################################

Results.Table <- NULL
Out.Var.Names <-c("log.HR.1", "log.HR.2", "log.HR.3", "rho", "log.odds.ratio.X.vs.W", "log.odds.ratio.X.vs.U", 
  "Cox.Est.1", "Cox.Est.2", "Cox.Est.3", "IV.Est.1", "IV.Est.2", "IV.Est.3") 
File_Name <- paste("C:/Users/Todd MacKenzie/Dropbox/ivhr time dependent/Sim Results/Stepwise TD Out-", Sys.Date(), ".txt", sep="")
while(1) {
	# Parameter setting
	log.odds.ratio.X.vs.U <- log(10)*(-1+2*runif(1)) 
	log.odds.ratio.X.vs.W <- log(5) + log(5)*runif(1)  # odds at least 5 up to 25, stength of the instrument
	HR <- sample(HR.Grid, 3)
	Settings <- c(log(HR), rho, log.odds.ratio.X.vs.W, log.odds.ratio.X.vs.U)
	# Instrument (continuous), omitted covariate (continuous) and exposure (binary)
	W <- rnorm(N)
	U <- rnorm(N) # omitted covariate or latent factor affecting survival
	X <- runif(N) < 1/(1+exp(-log.odds.ratio.X.vs.U * U - log.odds.ratio.X.vs.W * W))
	# Generation of the Potential Outcomes using Gaussian Copula
	Norm.0 <- rho * U + sqrt(1-rho^2) * rnorm(N)
	Norm.1 <- rho * U + sqrt(1-rho^2) * rnorm(N)
	Survival.0 <- Inverse.SDF.0(pnorm(Norm.0))
	Survival.1 <- Inverse.SDF.1.3step(pnorm(Norm.1), HR, TimeCuts)
	# Observed Potential Outcome
	Survival <- ifelse(X, Survival.1, Survival.0)
	# Censoring times; censoring fixed at 50%
	Censoring.Time <- runif(N)
	Censoring.Time <- Censoring.Time * quantile(Survival/Censoring.Time, Cens.Freq)
	# Now the observed time-to-event and status
	Time.to.Event <- pmin(Survival, Censoring.Time)
	Status <- ifelse(Survival < Censoring.Time, 1, 0)
	# Estimation by Cox
	pl1 <- coxph(Surv(Time.to.Event, Status & Time.to.Event < TimeCuts[1]) ~ X)
	pl2 <- coxph(Surv(Time.to.Event, Status & Time.to.Event < TimeCuts[2]) ~ X, subset=Survival > TimeCuts[1])
	pl3 <- coxph(Surv(Time.to.Event, Status) ~ X, subset=Time.to.Event > TimeCuts[2])
	iv1 <- IVHR(Time.to.Event, Status & Time.to.Event < TimeCuts[1], X, W)
	kp <- Time.to.Event >= TimeCuts[1]
	iv2 <- IVHR(Time.to.Event[kp], (Status & Time.to.Event < TimeCuts[2])[kp], X[kp], W[kp])
	kp <- Time.to.Event >= TimeCuts[2]
	iv3 <- IVHR(Time.to.Event[kp], Status[kp], X[kp], W[kp])
	Estimates <- c(pl1$coef, pl2$coef, pl3$coef, iv1$Est.log.HR, iv2$Est.log.HR, iv3$Est.log.HR)
	NewRow <- c(Settings, Estimates)
	Results.Table <- rbind(Results.Table, NewRow)
    if (dim(Results.Table)[1] %% 10000 == 0) write.table(Results.Table, File_Name, row.names=FALSE, col.names=Out.Var.Names)
	if (dim(Results.Table)[1] %% 100 == 0) cat(dim(Results.Table)[1], "\t")
}

dimnames(Results.Table)[[2]] <- Out.Var.Names 


# 
Results.Table <- read.delim("C:/Users/Todd MacKenzie/Dropbox/ivhr time dependent/Sim Results/Stepwise TD Out-2018-12-17.txt", sep=" ")

HR <- Results.Table[, c("log.HR.1", "log.HR.2", "log.HR.3")]
Cox.Est <- Results.Table[, c("Cox.Est.1", "Cox.Est.2", "Cox.Est.3")]
IV.Est <- Results.Table[, c("IV.Est.1", "IV.Est.2", "IV.Est.3")]

#png("C:/Users/Todd MacKenzie/Dropbox/ivhr time dependent/figures/sim_3_step.png")

par(mfrow=c(3,3))
FU.Name <- c("Early", "Middle", "Late")
for (i in 1:3) {
  keep <- Results.Table[, "log.odds.ratio.X.vs.U"] >= log(5)  
  plot(HR.Grid, exp(tapply(IV.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), log="xy", pch=16, xlab="True HR", ylab="Estimated HR", main=paste(FU.Name[i], ":", "+ Confounding"))
  points(HR.Grid, exp(tapply(Cox.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=2, pch=16)
  abline(a=0,b=1, lty=3)
  keep <- Results.Table[, "log.odds.ratio.X.vs.U"] > -log(2) & Results.Table[, "log.odds.ratio.X.vs.U"] < log(2)  
  plot(HR.Grid*exp(0.05), exp(tapply(IV.Est[keep,i], HR[keep,i], median, na.rm=TRUE))*exp(0.05), log="xy", pch=16,xlab="True HR", ylab="Estimated HR", main=paste(FU.Name[i], ":", "No Confounding"))
  points(HR.Grid*exp(-0.05), exp(tapply(Cox.Est[keep,i], HR[keep,i], median, na.rm=TRUE))*exp(-0.05), col=2, pch=16)
  abline(a=0,b=1, lty=3)
  keep <- Results.Table[, "log.odds.ratio.X.vs.U"] <= -log(5)  
  plot(HR.Grid, exp(tapply(IV.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), log="xy", pch=16, xlab="True HR", ylab="Estimated HR", main=paste(FU.Name[i], ":", "- Confounding"))
  points(HR.Grid, exp(tapply(Cox.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=2, pch=16)
  abline(a=0,b=1, lty=3)
}
dev.off()







