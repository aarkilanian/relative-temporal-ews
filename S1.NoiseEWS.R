##### Initial setup #####

# Inital parameters of model

N.init <- 100   # Initial population size
r <- 0.5        # Population growth rate
K <- 100        # Carrying capacity
h <- 1          # Half-saturation constant

c.fold <- 0.045   # Forcing level for fold model
c.trans <- 0.0025 # Forcing level for transcritical model

timestep <- 0.3 # Time step used in Euler approximation

tmax <- 500    # Max length of time series
tmin <- 0      # Min length of time series
reps <- 100    # Number of repetitions to run, script can handle up to 100

##### Simulate four populations, one with noise and one without for each bifurcation type #####

# Make containers to hold time series data

fold.noise <- vector("list", reps)
fold.nonoise <- vector("list", reps)
trans.noise <- vector("list", reps)
trans.nonoise <- vector("list", reps)

for (j in 1:(reps)) {
  
  # Setup empty vectors with appropriate size
  N.fold.noise <- numeric(tmin + tmax / timestep + 1)
  N.fold.nonoise <- numeric(tmin + tmax / timestep + 1)
  N.trans.noise <- numeric(tmin + tmax / timestep + 1)
  N.trans.nonoise <- numeric(tmin + tmax / timestep + 1)
  
  # Logicals monitoring if population is extinct
  fold.noise.dead <- FALSE
  fold.nonoise.dead <- FALSE
  trans.noise.dead <- FALSE
  trans.nonoise.dead <- FALSE
  
  curTime <- 0
  
  for (i in 1:length(N.fold.noise)) {
    
    if (i == 1) {
      
      N.fold.noise[i] <- N.init
      N.fold.nonoise[i] <- N.init
      N.trans.noise[i] <- N.init
      N.trans.nonoise[i] <- N.init
      
      # increment current timestep
      curTime <- 0.3
      
      next()
    }
    
    # Simulating fold without noise
    if(fold.nonoise.dead){} # If extinct don't simulate
    else {
      N.past <- N.fold.nonoise[i-1]
      
      growth <- N.past * r * (1 - (N.past / K)) * timestep 
      harvest <- (c.fold * curTime * N.past^2 / (N.past^2 + h^2)) * timestep
      
      N.cur <- N.past + growth - harvest
      
      # If extinct, keep it extinct and prevent going into negative pop size
      if (N.cur < 0) { 
        fold.nonoise.dead <- TRUE
        N.cur <- 0
      }
      
      # Store N
      N.fold.nonoise[i] <- N.cur
    }
    
    # Simulating fold with noise
    if(fold.noise.dead){} # If extinct don't simulate
    else {
      N.past <- N.fold.noise[i-1]
      
      growth <- N.past * r * (1 - (N.past / K)) * timestep 
      noise  <- rnorm(1,0,1.5) * timestep
      harvest <- (c.fold * curTime * N.past^2 / (N.past^2 + h^2)) * timestep
      
      N.cur <- N.past + growth + noise - harvest
      
      # If extinct, keep it extinct and prevent going into negative pop size
      if (N.cur < 0) { 
        fold.noise.dead <- TRUE
        N.cur <- 0
      }
      
      # Store N
      N.fold.noise[i] <- N.cur
    }
    
    # Simulating transcritical without noise
    if(trans.nonoise.dead){} # If extinct don't simulate
    else {
      N.past <- N.trans.nonoise[i-1]
      
      growth <- N.past * r * (1 - (N.past / K)) * timestep 
      harvest <- c.trans * curTime * N.past * timestep
      
      N.cur <- N.past + growth - harvest
      
      # If extinct, keep it extinct and prevent going into negative pop size
      if (N.cur < 0) { 
        trans.nonoise.dead <- TRUE
        N.cur <- 0
      }
      
      # Store N
      N.trans.nonoise[i] <- N.cur
    }
    
    # Simulating transcritical without noise
    if(trans.noise.dead){} # If extinct don't simulate
    else {
      N.past <- N.trans.noise[i-1]
      
      growth <- N.past * r * (1 - (N.past / K)) * timestep 
      noise  <- rnorm(1,0,1.5) * timestep
      harvest <- c.trans * curTime * N.past * timestep
      
      N.cur <- N.past + growth + noise - harvest
      
      # If extinct, keep it extinct and prevent going into negative pop size
      if (N.cur < 0) { 
        trans.noise.dead <- TRUE
        N.cur <- 0
      }
      
      # Store N
      N.trans.noise[i] <- N.cur
    }
    
    # increment current timestep
    curTime <- curTime + timestep
    
  }
  
  # Store replicates
  fold.noise[[j]] <- N.fold.noise
  fold.nonoise[[j]] <- N.fold.nonoise
  trans.noise[[j]] <- N.trans.noise
  trans.nonoise[[j]] <- N.trans.nonoise
  
}

##### EWS Analysis setup #####

require(ggplot2)
require(dplyr)
require(reshape2)
require(mgcv)

is.there <- function(x=0, L.CI, U.CI){ 
  pos <- ifelse(x<U.CI, 1, -1)
  negs<-ifelse(x>L.CI, -1, 1) 
  return(pos+negs)
}

gam_smoothing<-function(years, timeseries,knots){
  
  if(length(which(timeseries<=0))==0){
    
    gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian(link="log"))}else{
      
      gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian)}
  
  time.series.fit<-predict(gam1, newdata=data.frame(years=years), type="response")
  
  
  
  X0<-predict(gam1, newdata=data.frame(years=years), type= "lpmatrix")
  
  X1<-predict(gam1, newdata=data.frame(years=years), type= "lpmatrix")
  
  Xi<-(X1-X0)
  
  df <- Xi%*%coef(gam1)              ## ith smooth derivative 
  
  df.sd <- rowSums(Xi%*%gam1$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  #plot(years,df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)))##plot 'em
  
  #lines(years,df+2*df.sd,lty=2);lines(years,df-2*df.sd,lty=2)
  
  splines<-data.frame(years=years,deriv=df,U.CI=df+2*df.sd, L.CI=df-2*df.sd)	
  
  splines$sign<-is.there(0, splines$L.CI, splines$U.CI)/2
  
  splines$fit<-time.series.fit
  
  return(splines)}

EWS_analysis <- function (timeseries, winsize = 50, detrending = c("no", "gaussian", 
                                                                   "loess", "linear", "first-diff"), bandwidth = NULL, span = NULL, 
                          degree = NULL, logtransform = FALSE, interpolate = FALSE, 
                          AR_n = FALSE, powerspectrum = FALSE) {
  timeseries <- data.matrix(timeseries)
  if (dim(timeseries)[2] == 1) {
    Y = timeseries
    timeindex = 1:dim(timeseries)[1]
  }
  else if (dim(timeseries)[2] == 2) {
    Y <- timeseries[, 2]
    timeindex <- timeseries[, 1]
  }
  else {
    warning("not right format of timeseries input")
  }
  detrending <- match.arg(detrending)
  if (detrending == "gaussian") {
    if (is.null(bandwidth)) {
      bw <- round(bw.nrd0(timeindex))
    }
    else {
      bw <- round(length(Y) * bandwidth/100)
    }
    smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = bw, 
                    range.x = range(timeindex), x.points = timeindex)
    nsmY <- Y - smYY$y
    smY <- smYY$y
  }
  else if (detrending == "linear") {
    nsmY <- resid(lm(Y ~ timeindex))
    smY <- fitted(lm(Y ~ timeindex))
  }
  else if (detrending == "no"){
    smY <- Y
    nsmY <- Y
  }
  mw <- round(length(Y) * winsize/100)
  omw <- length(nsmY) - mw + 1
  nMR <- matrix(data = NA, nrow = mw, ncol = omw)
  
  for (i in 1:omw) {
    Ytw <- nsmY[i:(i + mw - 1)]
    nMR[, i] <- Ytw
  }
  nARR <- numeric()
  nSD <- numeric()
  nSD <- apply(nMR, 2, sd, na.rm = TRUE)
  for (i in 1:ncol(nMR)) {
    nYR <- ar.ols(nMR[, i], aic = FALSE, order.max = 1, dmean = FALSE, 
                  intercept = FALSE)
    nARR[i] <- nYR$ar
  }
  
  timevec <- seq(1, length(nARR))
  KtAR <- cor.test(timevec, nARR, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  KtSD <- cor.test(timevec, nSD, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  
  out <- c(KTauSD = round(KtSD$estimate, digits = 3), 
           KTauAr1 = round(KtAR$estimate,digits = 3))
  return(out)
}

##### Fold model subsetting #####

fold.noise.sub <- vector("list", reps)
fold.nonoise.sub <- vector("list", reps)

for (i in 1:reps) {
  
  # Find bounds for subsetting
  lower.noise <- 1
  
  smoothy <- gam_smoothing(1:1667, fold.noise[[i]], -1)$fit
  change <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
  upper.noise <- which(change < -0.01)[[1]]
  
  lower.nonoise <- 1
  
  smoothy <- gam_smoothing(1:1667, fold.nonoise[[i]], -1)$fit
  change <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
  upper.nonoise <- which(change < -0.01)[[1]]
  
  # Subset time series to not include abundance after bifurcation
  fold.noise.sub[[i]] <- fold.noise[[i]][lower.noise:upper.noise]
  fold.nonoise.sub[[i]] <- fold.nonoise[[i]][lower.nonoise:upper.nonoise]
  
}

##### Transcritical model subsetting #####

trans.noise.sub <- vector("list", reps)
trans.nonoise.sub <- vector("list", reps)

for (i in 1:reps) {
  
  # Find bounds for subsetting
  lower.noise <- 1
  
  smoothy <- gam_smoothing(1:1667, trans.noise[[i]], -1)$fit
  change <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
  upper.noise <- which(change < -0.05)[[1]]
  
  lower.nonoise <- 1
  
  smoothy <- gam_smoothing(1:1667, trans.nonoise[[i]], -1)$fit
  change <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
  upper.nonoise <- which(change < -0.05)[[1]]
  
  # Subset time series to not include abundance after bifurcation
  trans.noise.sub[[i]] <- trans.noise[[i]][lower.noise:upper.noise]
  trans.nonoise.sub[[i]] <- trans.nonoise[[i]][lower.nonoise:upper.nonoise]
  
}

##### EWS analysis #####

sd.fold.noise <- numeric(100)
sd.fold.nonoise <- numeric(100)
sd.trans.noise <- numeric(100)
sd.trans.nonoise <- numeric(100)

ar.fold.noise <- numeric(100)
ar.fold.nonoise <- numeric(100)
ar.trans.noise <- numeric(100)
ar.trans.nonoise <- numeric(100)

for (i in 1:reps) {
  
  # fold with noise
  res <- EWS_analysis(fold.noise.sub[[i]], detrending = "gaussian")
  sd.fold.noise[i] <- res[1]
  ar.fold.noise[i] <- res[2]
  
  # fold without noise
  res <- EWS_analysis(fold.nonoise.sub[[i]], detrending = "gaussian")
  sd.fold.nonoise[i] <- res[1]
  ar.fold.nonoise[i] <- res[2]
  
  # trans with noise
  res <- EWS_analysis(trans.noise.sub[[i]], detrending = "gaussian")
  sd.trans.noise[i] <- res[1]
  ar.trans.noise[i] <- res[2]
  
  # trans without noise
  res <- EWS_analysis(trans.nonoise.sub[[i]], detrending = "gaussian")
  sd.trans.nonoise[i] <- res[1]
  ar.trans.nonoise[i] <- res[2]
  
}

##### Plot results in line graphs #####

par(mfrow = c(2,1))

# Standard deviation plots
plot(sd.fold.noise, col = 'blue', typ = 'l', ylim = c(-1,1),
     main = "Standard deviation based early warning signals", ylab = "Kendall's tau (fold model)", xlab = "Replicates")
lines(sd.fold.nonoise, col = "black")
abline(h = 0, col = 'red', lty = 3)

plot(sd.trans.noise, col = 'blue', typ = 'l', ylim = c(-1,1),
     ylab = "Kendall's tau (trans model)", xlab = "Replicates")
lines(sd.trans.nonoise, col = "black")
abline(h = 0, col = 'red', lty = 3)

# Autocorrelation plots
plot(ar.fold.noise, col = 'blue', typ = 'l', ylim = c(-1,1),
     main = "Autocorrelation based early warning signals", ylab = "Kendall's tau (fold model)", xlab = "Replicates")
lines(ar.fold.nonoise, col = "black")
abline(h = 0, col = 'red', lty = 3)

plot(ar.trans.noise, col = 'blue', typ = 'l', ylim = c(-1,1),
     ylab = "Kendall's tau (trans model)", xlab = "Replicates")
lines(ar.trans.nonoise, col = "black")
abline(h = 0, col = 'red', lty = 3)

##### Plot results in boxplots #####

res.frame <- data.frame(sd.fold.noise, sd.fold.nonoise, ar.fold.noise, ar.fold.nonoise,
                        sd.trans.noise, sd.trans.nonoise, ar.trans.noise, ar.trans.nonoise)


