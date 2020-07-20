
library(survival)
library(MASS)
library(rootSolve)
library(copula)
library(coxme)

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

Wang.Estimator = function(Y, Delta, D, Z)
{
        n<- length(D)           
        delta.D = ifelse( sd(D-Z)== 0, 1, mean(D*Z)/mean(Z) - mean(D*(1-Z))/mean(1-Z))
        
        part1 = (2*Z - 1) / delta.D
        part2 = sapply(1:n, function(i) sum(D * as.numeric(Y >= Y[i]) * part1))
        part3 = sapply(1:n, function(i) sum(as.numeric(Y >= Y[i]) * part1))
        
        
        KM = survfit(Surv(Y, Delta) ~ D, type="kaplan-meier")
        surv.diff = diff(summary(KM)$table[,"*rmean"])
        
        # crude.hr = coxph(Surv(Y, Delta) ~ D + X)$coef[1]
        
        if(surv.diff < 0){
           e.psi = sum(Delta * D * part1 * (1-part2/(3*part3 - 2*part2))) /
                sum(3 * Delta * (1-D) * part1 * part2/(3*part3 - 2*part2))

            if(!is.finite(e.psi)){
                index = (3*part3 - 2*part2) != 0
                e.psi = sum((Delta * D * part1 * (1-part2/(3*part3 - 2*part2)))[index]) /
                    sum((3 * Delta * (1-D) * part1 * part2/(3*part3 - 2*part2))[index])
            }
            k = 2
            while(is.finite(e.psi) & (e.psi < 1e-6 | e.psi > 1e6) & k<100){
                k = k + 1
                # g(D) = (k+1)-k*D
                e.psi = sum(Delta * D * part1 * (1-part2/((k+1)*part3 - k*part2))) / 
                    sum((k+1) * Delta * (1-D) * part1 * part2/((k+1)*part3 - k*part2))
            }
        }else{
            # g(D) = 2D+1
            e.psi = sum(3 * Delta * D * part1 * (1-3*part2/(2*part2 + part3))) / 
                sum(Delta * (1-D) * part1 * 3 * part2 / (2*part2 + part3))
            if(!is.finite(e.psi)){
                index = (2*part2 + part3) != 0
                e.psi = sum((3*Delta * D*part1* (1-3*part2/(2*part2 + part3)))[index]) / 
                    sum((Delta * (1-D) * part1 * 3 * part2 / (2*part2 + part3))[index])
            }
            k = 2
            while(is.finite(e.psi) & (e.psi < 1e-6 | e.psi > 1e6) & k<100){
                k = k + 1
                # g(D) = k*D + 1
                e.psi = sum((k+1) * Delta * D * part1 * (1-(k+1)*part2/(k*part2 + part3))) / 
                    sum( Delta * (1-D) * part1 * (k+1) * part2/(k*part2 + part3))
            }
        }
        
        return(log(e.psi))       
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


TimeCuts <- -log(1 - c(0.25, 0.50))
n.periods <- 1 + length(TimeCuts)
HR.Grid <- exp(seq(from=log(1/3), to=log(3), length=25))

##########################################################
# Set sample size at 2000, and cens freq at 50% for 
N <- 400
Cens.Freq <- 0.50
log.odds.ratio.X.vs.Om.grid <- log(c(1/5, 1, 5))

##########################################################

Results.Table <- NULL
Out.Var.Names <-c("log.HR.1", "log.HR.2", "log.HR.3", "rho", "copula.type", "log.odds.ratio.X.vs.W", "log.odds.ratio.X.vs.Om", 
  "Cox.Est.1", "Cox.Est.2", "Cox.Est.3", "IV.Est.1", "IV.Est.2", "IV.Est.3", "Wang.Est.1", "Wang.Est.2", "Wang.Est.3", "RIF.Est.1", "RIF.Est.2", "RIF.Est.3") 
File_Name <- paste("Sim Results/Stepwise TD Out-", Sys.Date(), ".txt", sep="")
while(1) {
	# Parameter setting
	rho <- ifelse(runif(1) < 0.5, 0.99, 0.50)
	copula.type <- ceiling(3*runif(1))
	log.odds.ratio.X.vs.Om <- log(10)*(-1+2*runif(1)) 
	log.odds.ratio.X.vs.W <- log(2) + log(25)*runif(1)  # odds at least 2 up to 50, stength of the instrument
	HR <- sample(HR.Grid, 3)
	Settings <- c(log(HR), rho, copula.type, log.odds.ratio.X.vs.W, log.odds.ratio.X.vs.Om)
	# Generation of the Potential Outcomes using Gaussian Copula
	if (copula.type == 1) {
  	UY <- rCopula(N, normalCopula(rho, dim=2))
	}
	# Generation of Clayton copula 
	if (copula.type == 2) {
	  Clayton.param <- 2 * rho / (1 - rho) # Kendall corr
	  UY <- rCopula(N, claytonCopula(Clayton.param, dim=2))
	}
	if (copula.type == 3) {
	  Gumbel.param <- 1/(1 - rho) # Kendall corr
	  UY <- rCopula(N, gumbelCopula(Gumbel.param, dim=2))
	}

	Survival.0 <- Inverse.SDF.0(UY[,1])
	Survival.1 <- Inverse.SDF.1.3step(UY[,2], HR, TimeCuts)
	# Omitted covariate or latent factor affecting survival
	Om <- qnorm(UY[,1])
	# Instrument (continuous), omitted covariate (continuous) and exposure (binary)
	W <- rnorm(N)
	X <- runif(N) < 1/(1+exp(-log.odds.ratio.X.vs.Om * Om - log.odds.ratio.X.vs.W * W))
	# Observed Potential Outcome
	Survival <- ifelse(X, Survival.1, Survival.0)
	# Censoring times; censoring fixed at 50%
	Censoring.Time <- runif(N)
	Censoring.Time <- Censoring.Time * quantile(Survival/Censoring.Time, 1-Cens.Freq)
	# Now the observed time-to-event and status
	Time.to.Event <- pmin(Survival, Censoring.Time)
	Status <- ifelse(Survival < Censoring.Time, 1, 0)
	# Estimation by Cox
	pl1 <- coxph(Surv(Time.to.Event, Status & Time.to.Event < TimeCuts[1]) ~ X)
	pl2 <- coxph(Surv(Time.to.Event, Status & Time.to.Event < TimeCuts[2]) ~ X, subset=Survival > TimeCuts[1])
	pl3 <- coxph(Surv(Time.to.Event, Status) ~ X, subset=Time.to.Event > TimeCuts[2])
  Cox.Est <- c(pl1$coef, pl2$coef, pl3$coef)
	# Estimation by our approach	
	iv1 <- IVHR(Time.to.Event, Status & Time.to.Event < TimeCuts[1], X, W)
	kp <- Time.to.Event >= TimeCuts[1]
	iv2 <- IVHR(Time.to.Event[kp], (Status & Time.to.Event < TimeCuts[2])[kp], X[kp], W[kp])
	kp <- Time.to.Event >= TimeCuts[2]
	iv3 <- IVHR(Time.to.Event[kp], Status[kp], X[kp], W[kp])
	IV.Est <- c(iv1$Est.log.HR, iv2$Est.log.HR, iv3$Est.log.HR)
	# Linbo Wang approach
	WW <- W > 0 # Dichotomize W since Wang is for binary instrument
	wang1 <- Wang.Estimator(Time.to.Event, Status & Time.to.Event < TimeCuts[1], X, WW)
	kp <- Time.to.Event >= TimeCuts[1]
	wang2 <- Wang.Estimator(Time.to.Event[kp], (Status & Time.to.Event < TimeCuts[2])[kp], X[kp], WW[kp])
	kp <- Time.to.Event >= TimeCuts[2]
	wang3 <- Wang.Estimator(Time.to.Event[kp], Status[kp], X[kp], WW[kp])
	Wang.IV.Est <- c(wang1, wang2, wang3)
	# Estimation by residual inclusion & frailty other approach; do 1/3 of time
  RIF.Estimates <- rep(NA, 3)
	if (runif(1)<0.25) {	
  	Res <- lm(X ~ W)$res
  	ID <- 1:N
  	RIF1 <- coxme(Surv(Time.to.Event, Status & Time.to.Event < TimeCuts[1]) ~ X + Res  + (1|ID))$coef[1]
  	RIF2 <- coxme(Surv(Time.to.Event, Status & Time.to.Event < TimeCuts[2]) ~ X + Res  + (1|ID), subset=Survival > TimeCuts[1])$coef[1]
  	RIF3 <- coxme(Surv(Time.to.Event, Status) ~ X + Res  + (1|ID), subset=Time.to.Event > TimeCuts[2])$coef[1]
  	RIF.Estimates <- c(RIF1, RIF2, RIF3)
  }

	NewRow <- c(Settings, Cox.Est, IV.Est, Wang.IV.Est, RIF.Estimates)
	Results.Table <- rbind(Results.Table, NewRow)
    if (dim(Results.Table)[1] %% 10000 == 0) write.table(Results.Table, File_Name, row.names=FALSE, col.names=Out.Var.Names)
	if (dim(Results.Table)[1] %% 100 == 0) cat(dim(Results.Table)[1], "\t")
}

dimnames(Results.Table)[[2]] <- Out.Var.Names 

#Results.Table <- read.delim("Sim Results/Stepwise TD Out-2020-02-11.txt", sep=" ")

HR <- Results.Table[, c("log.HR.1", "log.HR.2", "log.HR.3")]
Cox.Est <- Results.Table[, c("Cox.Est.1", "Cox.Est.2", "Cox.Est.3")]
IV.Est <- Results.Table[, c("IV.Est.1", "IV.Est.2", "IV.Est.3")]
RIF.Est <- Results.Table[, c("RIF.Est.1", "RIF.Est.2", "RIF.Est.3")]
Wang.Est <- Results.Table[, c("Wang.Est.1", "Wang.Est.2", "Wang.Est.3")]

# Modelling the absolute bias with respect to simulation variables
abs.diff1 <- abs(Results.Table[,"IV.Est.1"] - Results.Table[, "log.HR.1"])
abs.diff2 <- abs(Results.Table[,"IV.Est.2"] - Results.Table[, "log.HR.2"])
abs.diff3 <- abs(Results.Table[,"IV.Est.3"] - Results.Table[, "log.HR.3"])


summary(lm(abs.diff1 ~ log.odds.ratio.X.vs.W + log.odds.ratio.X.vs.Om + rho + log.HR.1, data=data.frame(Results.Table)))$coef
summary(lm(abs.diff2 ~ log.odds.ratio.X.vs.W + log.odds.ratio.X.vs.Om + rho + log.HR.2, data=data.frame(Results.Table)))$coef
summary(lm(abs.diff3 ~ log.odds.ratio.X.vs.W + log.odds.ratio.X.vs.Om + rho + log.HR.3, data=data.frame(Results.Table)))$coef

Figure <- function(Subset=NULL, fig.file.name=NULL) {
  if (is.null(Subset)) Subset <- rep(TRUE, dim(HR)[1])
  if (!is.null(fig.file.name))
    png(paste("Figures/", fig.file.name, " ", Sys.Date(), ".png", sep=""), units="in", width=5, height=5, res=300)
  par(mfrow=c(3,3))
  FU.Name <- c("Early", "Middle", "Late")
  y.lim <- c(0.1, 10)
  for (i in 1:3) {
    keep <- Subset & Results.Table[, "log.odds.ratio.X.vs.Om"] >= log(5)  
    par(mar=c(ifelse(i==3,4,3),4,3,1))
    plot(HR.Grid*0.97, exp(tapply(IV.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), log="xy", pch=16, xlab="True HR", ylab="Estimated HR", main=paste(FU.Name[i], ":", "+ Confounding", sep=""), ylim=y.lim)
    points(HR.Grid*0.99, exp(tapply(Cox.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=2)
    points(HR.Grid*1.01, exp(tapply(Wang.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=3)
    points(HR.Grid*1.03, exp(tapply(RIF.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=4)
    abline(a=0,b=1, lty=3, lwd=0.7)
    keep <- Subset & Results.Table[, "log.odds.ratio.X.vs.Om"] > -log(2) & Results.Table[, "log.odds.ratio.X.vs.Om"] < log(2)  
    par(mar=c(ifelse(i==3,4,3),3,3,1))
    plot(HR.Grid*0.97, exp(tapply(IV.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), log="xy", pch=16,xlab="True HR", ylab="Estimated HR", main=paste(FU.Name[i], ":", "No Confounding", sep=""), ylim=y.lim)
    points(HR.Grid*0.99, exp(tapply(Cox.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=2)
    points(HR.Grid*1.01, exp(tapply(Wang.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=3)
    points(HR.Grid*1.03, exp(tapply(RIF.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=4)
    abline(a=0,b=1, lty=3, lwd=0.7)
    keep <- Subset & Results.Table[, "log.odds.ratio.X.vs.Om"] <= -log(5)  
    par(mar=c(ifelse(i==3,4,3),3,3,1))
    plot(HR.Grid*0.97, exp(tapply(IV.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), log="xy", pch=16, xlab="True HR", ylab="Estimated HR", main=paste(FU.Name[i], ":", "- Confounding", sep=""), ylim=y.lim)
    points(HR.Grid*0.99, exp(tapply(Cox.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=2)
    points(HR.Grid*1.01, exp(tapply(Wang.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=3)
    points(HR.Grid*1.03, exp(tapply(RIF.Est[keep,i], HR[keep,i], median, na.rm=TRUE)), col=4)
    abline(a=0,b=1, lty=3, lwd=0.7)
  }
  if (!is.null(fig.file.name)) dev.off()
}  

# By strength of instrument
Figure(Results.Table[, "log.odds.ratio.X.vs.W"]<log(5), fig.file.name = "IV O.R. lt 5")
Figure(Results.Table[, "log.odds.ratio.X.vs.W"]>=log(5), fig.file.name = "IV O.R. ge 5")



Figure(fig.file.name = "All")

Cop.Rho <- paste(c("Normal","Clayton", "Gumbel")[Results.Table[,"copula.type"]], Results.Table[,"rho"])
Figure(Cop.Rho=="Normal 0.5", fig.file.name = "Normal  r=.50")
Figure(Cop.Rho=="Normal 0.99", fig.file.name = "Normal  r=.99")
Figure(Cop.Rho=="Clayton 0.5", fig.file.name = "Clayton  r=.50")
Figure(Cop.Rho=="Clayton 0.99", fig.file.name = "Clayton  r=.99")
Figure(Cop.Rho=="Gumbel 0.5", fig.file.name = "Gumbel  r=.50")
Figure(Cop.Rho=="Gumbel 0.99", fig.file.name = "Gumbel  r=.99")







