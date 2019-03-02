# 2.
# This script sets up simulation data for analysis of the effect of length of 
# sampling on EWS across all sampling resolutions for strongest and weakest
# forcing levels
#
# Code developed by Alex Arkilanian and Gaurav Baruah
#--------------------------------------------------------------------------

# Clean
rm(list = ls())

##### Packages & Functions #####

load("AllPop.RData")

reps <- 100

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Fold models #####

## Select and subset data #####
none.2 <- vector("list", reps)
weak.2 <- vector("list", reps)
medium.2 <- vector("list", reps)
strong.2 <- vector("list", reps)

none.1 <- vector("list", reps)
weak.1 <- vector("list", reps)
medium.1 <- vector("list", reps)
strong.1 <- vector("list", reps)

none.50 <- vector("list", reps)
weak.50 <- vector("list", reps)
medium.50 <- vector("list", reps)
strong.50 <- vector("list", reps)

none.25 <- vector("list", reps)
weak.25 <- vector("list", reps)
medium.25 <- vector("list", reps)
strong.25 <- vector("list", reps)

for(i in 1:reps){

N.fold     <- Ns.fold[[i]]
N.fold.2   <- Ns.fold.2[[i]]
N.fold.1   <- Ns.fold.1[[i]]
N.fold.50  <- Ns.fold.50[[i]]
N.fold.25  <- Ns.fold.25[[i]]

# Subset full resolution data by forcing (strong and weak) - new
N.pop.none   <- N.fold[1:1667]
N.pop.weak   <- N.fold[1668:3334]
N.pop.medium <- N.fold[3335:5001]
N.pop.strong <- N.fold[5002:6668]

# Lower bounds
lower.none      <- 1
lower.weak      <- 1668
lower.medium    <- 3335
lower.strong    <- 5002

lower.2.none   <- 1
lower.2.weak   <- ceiling(1668 * (0.025/(1.3 * 2))) + 1
lower.2.medium <- ceiling(3335 * (0.025/(1.3 * 2))) + 2
lower.2.strong <- ceiling(5002 * (0.025/(1.3 * 2))) + 3

lower.1.none   <- 1
lower.1.weak   <- ceiling(1668 * (0.025/(1.3 * 1))) + 1
lower.1.medium <- ceiling(3335 * (0.025/(1.3 * 1))) + 2
lower.1.strong <- ceiling(5002 * (0.025/(1.3 * 1))) + 3

lower.50.none   <- 1
lower.50.weak   <- ceiling(1668 * (0.025/(1.3 * 0.5))) + 1
lower.50.medium <- ceiling(3335 * (0.025/(1.3 * 0.5))) + 2
lower.50.strong <- ceiling(5002 * (0.025/(1.3 * 0.5))) + 3

lower.25.none   <- 1
lower.25.weak   <- ceiling(1668 * (0.025/(1.3 * 0.25))) + 1
lower.25.medium <- ceiling(3335 * (0.025/(1.3 * 0.25))) + 2
lower.25.strong <- ceiling(5002 * (0.025/(1.3 * 0.25))) + 3

# Calculate upper bounds - new
upper.none <- 1667

smoothy  <- gam_smoothing(1:1667, N.pop.weak, -1)$fit
change   <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
upper.weak <- which(change < -0.01)[[1]]

smoothy  <- gam_smoothing(1:1667, N.pop.medium, -1)$fit
change   <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
upper.medium <- which(change < -0.01)[[1]]

smoothy  <- gam_smoothing(1:1667, N.pop.strong, -1)$fit
change   <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
upper.strong <- which(change < -0.01)[[1]]

# Upper bounds - Conversion: convert to time points and then to corresponding res

upper.2.none    <- floor(upper.none * 0.025 / (1.3 * 2))
upper.2.weak    <- floor(upper.weak * 0.025 / (1.3 * 2)) + 1
upper.2.medium  <- floor(upper.medium * 0.025 / (1.3 * 2)) + 2
upper.2.strong  <- floor(upper.strong * 0.025 / (1.3 * 2)) + 3

upper.1.none    <- floor(upper.none * 0.025 / (1.3 * 1))
upper.1.weak    <- floor(upper.weak * 0.025 / (1.3 * 1)) + 1
upper.1.medium  <- floor(upper.medium * 0.025 / (1.3 * 1)) + 2
upper.1.strong  <- floor(upper.strong * 0.025 / (1.3 * 1)) + 3

upper.50.none   <- floor(upper.none * 0.025 / (1.3 * 0.5))
upper.50.weak   <- floor(upper.weak * 0.025 / (1.3 * 0.5)) + 1
upper.50.medium <- floor(upper.medium * 0.025 / (1.3 * 0.5)) + 2
upper.50.strong <- floor(upper.strong * 0.025 / (1.3 * 0.5)) + 3

upper.25.none   <- floor(upper.none * 0.025 / (1.3 * 0.25))
upper.25.weak   <- floor(upper.weak * 0.025 / (1.3 * 0.25)) + 1
upper.25.medium <- floor(upper.medium * 0.025 / (1.3 * 0.25)) + 2
upper.25.strong <- floor(upper.strong * 0.025 / (1.3 * 0.25)) + 3

# Subset data
none.2[[i]]    <- N.fold.2[lower.2.none:(lower.2.none + upper.2.none)]
weak.2[[i]]    <- N.fold.2[lower.2.weak:(lower.2.weak + upper.2.weak)]
medium.2[[i]]  <- N.fold.2[lower.2.medium:(lower.2.medium + upper.2.medium)]
strong.2[[i]]  <- N.fold.2[lower.2.strong:(lower.2.strong + upper.2.strong)]

none.1[[i]]    <- N.fold.1[lower.1.none:(lower.1.none + upper.1.none)]
weak.1[[i]]    <- N.fold.1[lower.1.weak:(lower.1.weak + upper.1.weak)]
medium.1[[i]]  <- N.fold.1[lower.1.medium:(lower.1.medium + upper.1.medium)]
strong.1[[i]]  <- N.fold.1[lower.1.strong:(lower.1.strong + upper.1.strong)]

none.50[[i]]    <- N.fold.50[lower.50.none:(lower.50.none + upper.50.none)]
weak.50[[i]]   <- N.fold.50[lower.50.weak:(lower.50.weak + upper.50.weak)]
medium.50[[i]]  <- N.fold.50[lower.50.medium:(lower.50.medium + upper.50.medium)]
strong.50[[i]] <- N.fold.50[lower.50.strong:(lower.50.strong + upper.50.strong)]

none.25[[i]]    <- N.fold.25[lower.25.none:(lower.25.none + upper.25.none)]
weak.25[[i]]   <- N.fold.25[lower.25.weak:(lower.25.weak + upper.25.weak)]
medium.25[[i]]  <- N.fold.25[lower.25.medium:(lower.25.medium + upper.25.medium)]
strong.25[[i]] <- N.fold.25[lower.25.strong:(lower.25.strong + upper.25.strong)]

}

## 1/4 G Sampling #####

# Strong forcing
m <- length(strong.25[[i]]) - 7

res.strong.25.sd <- matrix(ncol = m, nrow = reps)
res.strong.25.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(strong.25[[i]][1:(length(strong.25[[i]])-j)], detrending = "gaussian")
    res.strong.25.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.strong.25.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Medium forcing
m <- length(medium.25[[1]]) - 7

res.medium.25.sd <- matrix(ncol = m, nrow = reps)
res.medium.25.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(medium.25[[i]][1:(length(medium.25[[i]])-j)], detrending = "gaussian")
    res.medium.25.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.medium.25.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Weak forcing
m <- length(weak.25[[i]]) - 7

res.weak.25.sd <- matrix(ncol = m, nrow = reps)
res.weak.25.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(weak.25[[i]][1:(length(weak.25[[i]])-j)], detrending = "gaussian")
    res.weak.25.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.weak.25.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

## 1/2 G Sampling #####

# Strong forcing
m <- length(strong.50[[i]]) - 7

res.strong.50.sd <- matrix(ncol = m, nrow = reps)
res.strong.50.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(strong.50[[i]][1:(length(strong.50[[i]])-j)], detrending = "gaussian")
    res.strong.50.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.strong.50.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Medium forcing
m <- length(medium.50[[i]]) - 7

res.medium.50.sd <- matrix(ncol = m, nrow = reps)
res.medium.50.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(medium.50[[i]][1:(length(medium.50[[i]])-j)], detrending = "gaussian")
    res.medium.50.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.medium.50.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Weak forcing
m <- length(weak.50[[i]]) - 7

res.weak.50.sd <- matrix(ncol = m, nrow = reps)
res.weak.50.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(weak.50[[i]][1:(length(weak.50[[i]])-j)], detrending = "gaussian")
    res.weak.50.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.weak.50.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

## 1 G Sampling #####

# Strong forcing
m <- length(strong.1[[i]]) - 7

res.strong.1.sd <- matrix(ncol = m, nrow = reps)
res.strong.1.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(strong.1[[i]][1:(length(strong.1[[i]])-j)], detrending = "gaussian")
    res.strong.1.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.strong.1.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Medium forcing
m <- length(medium.1[[i]]) - 7

res.medium.1.sd <- matrix(ncol = m, nrow = reps)
res.medium.1.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(medium.1[[i]][1:(length(medium.1[[i]])-j)], detrending = "gaussian")
    res.medium.1.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.medium.1.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Weak forcing
m <- length(weak.1[[i]]) - 7

res.weak.1.sd <- matrix(ncol = m, nrow = reps)
res.weak.1.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(weak.1[[i]][1:(length(weak.1[[i]])-j)], detrending = "gaussian")
    res.weak.1.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.weak.1.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

## 2 G Sampling #####

# Strong forcing
m <- length(strong.2[[i]]) - 7

res.strong.2.sd <- matrix(ncol = m, nrow = reps)
res.strong.2.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(strong.2[[i]][1:(length(strong.2[[i]])-j)], detrending = "gaussian")
    res.strong.2.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.strong.2.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Medium forcing
m <- length(medium.2[[i]]) - 7

res.medium.2.sd <- matrix(ncol = m, nrow = reps)
res.medium.2.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(medium.1[[i]][1:(length(medium.1[[i]])-j)], detrending = "gaussian")
    res.medium.2.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.medium.2.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Weak forcing
m <- length(weak.2[[i]]) - 7

res.weak.2.sd <- matrix(ncol = m, nrow = reps)
res.weak.2.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(weak.2[[i]][1:(length(weak.2[[i]])-j)], detrending = "gaussian")
    res.weak.2.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.weak.2.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

## Save Results #####

save(res.strong.25.au, res.strong.50.au, res.strong.1.au, res.strong.2.au,
     res.strong.25.sd, res.strong.50.sd, res.strong.1.sd, res.strong.2.sd,
     res.medium.25.au, res.medium.50.au, res.medium.1.au, res.medium.2.au,
     res.medium.25.sd, res.medium.50.sd, res.medium.1.sd, res.medium.2.sd,
     res.weak.25.au, res.weak.50.au, res.weak.1.au, res.weak.2.au,
     res.weak.25.sd, res.weak.50.sd, res.weak.1.sd, res.weak.2.sd,
     file = "analyzedFold.RData")

##### Trans models #####

## Select and subset data #####
none.2 <- vector("list", reps)
weak.2 <- vector("list", reps)
medium.2 <- vector("list", reps)
strong.2 <- vector("list", reps)

none.1 <- vector("list", reps)
weak.1 <- vector("list", reps)
medium.1 <- vector("list", reps)
strong.1 <- vector("list", reps)

none.50 <- vector("list", reps)
weak.50 <- vector("list", reps)
medium.50 <- vector("list", reps)
strong.50 <- vector("list", reps)

none.25 <- vector("list", reps)
weak.25 <- vector("list", reps)
medium.25 <- vector("list", reps)
strong.25 <- vector("list", reps)

for(i in 1:reps){
  
  N.trans     <- Ns.trans[[i]]
  N.trans.2   <- Ns.trans.2[[i]]
  N.trans.1   <- Ns.trans.1[[i]]
  N.trans.50  <- Ns.trans.50[[i]]
  N.trans.25  <- Ns.trans.25[[i]]
  
  # Subset full resolution data by forcing (strong and weak) - new
  N.pop.none   <- N.trans[1:1667]
  N.pop.weak   <- N.trans[1668:3334]
  N.pop.medium <- N.trans[3335:5001]
  N.pop.strong <- N.trans[5002:6668]
  
  # Lower bounds
  lower.none      <- 1
  lower.weak      <- 1668
  lower.medium    <- 3335
  lower.strong    <- 5002
  
  lower.2.none   <- 1
  lower.2.weak   <- ceiling(1668 * (0.025/(1.4 * 2))) + 1
  lower.2.medium <- ceiling(3335 * (0.025/(1.4 * 2))) + 2
  lower.2.strong <- ceiling(5002 * (0.025/(1.4 * 2))) + 3
  
  lower.1.none   <- 1
  lower.1.weak   <- ceiling(1668 * (0.025/(1.4 * 1))) + 1
  lower.1.medium <- ceiling(3335 * (0.025/(1.4 * 1))) + 2
  lower.1.strong <- ceiling(5002 * (0.025/(1.4 * 1))) + 3
  
  lower.50.none   <- 1
  lower.50.weak   <- ceiling(1668 * (0.025/(1.4 * 0.5))) + 1
  lower.50.medium <- ceiling(3335 * (0.025/(1.4 * 0.5))) + 2
  lower.50.strong <- ceiling(5002 * (0.025/(1.4 * 0.5))) + 3
  
  lower.25.none   <- 1
  lower.25.weak   <- ceiling(1668 * (0.025/(1.4 * 0.25))) + 1
  lower.25.medium <- ceiling(3335 * (0.025/(1.4 * 0.25))) + 2
  lower.25.strong <- ceiling(5002 * (0.025/(1.4 * 0.25))) + 3
  
  # Calculate upper bounds - new
  upper.none <- 1667
  
  smoothy  <- gam_smoothing(1:1667, N.pop.weak, -1)$fit
  change   <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
  upper.weak <- which(change < -0.05)[[1]]
  
  smoothy  <- gam_smoothing(1:1667, N.pop.medium, -1)$fit
  change   <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
  upper.medium <- which(change < -0.05)[[1]]
  
  smoothy  <- gam_smoothing(1:1667, N.pop.strong, -1)$fit
  change   <- (smoothy[2:1667] - smoothy[1:1666]) / smoothy[2:1667]
  upper.strong <- which(change < -0.05)[[1]]
  
  # Upper bounds - Conversion: convert to time points and then to corresponding res
  
  upper.2.none    <- floor(upper.none * 0.025 / (1.4 * 2))
  upper.2.weak    <- floor(upper.weak * 0.025 / (1.4 * 2)) + 1
  upper.2.medium  <- floor(upper.medium * 0.025 / (1.4 * 2)) + 2
  upper.2.strong  <- floor(upper.strong * 0.025 / (1.4 * 2)) + 3
  
  upper.1.none    <- floor(upper.none * 0.025 / (1.4 * 1))
  upper.1.weak    <- floor(upper.weak * 0.025 / (1.4 * 1)) + 1
  upper.1.medium  <- floor(upper.medium * 0.025 / (1.4 * 1)) + 2
  upper.1.strong  <- floor(upper.strong * 0.025 / (1.4 * 1)) + 3
  
  upper.50.none   <- floor(upper.none * 0.025 / (1.4 * 0.5))
  upper.50.weak   <- floor(upper.weak * 0.025 / (1.4 * 0.5)) + 1
  upper.50.medium <- floor(upper.medium * 0.025 / (1.4 * 0.5)) + 2
  upper.50.strong <- floor(upper.strong * 0.025 / (1.4 * 0.5)) + 3
  
  upper.25.none   <- floor(upper.none * 0.025 / (1.4 * 0.25))
  upper.25.weak   <- floor(upper.weak * 0.025 / (1.4 * 0.25)) + 1
  upper.25.medium <- floor(upper.medium * 0.025 / (1.4 * 0.25)) + 2
  upper.25.strong <- floor(upper.strong * 0.025 / (1.4 * 0.25)) + 3
  
  # Subset data
  none.2[[i]]    <- N.trans.2[lower.2.none:(lower.2.none + upper.2.none)]
  weak.2[[i]]    <- N.trans.2[lower.2.weak:(lower.2.weak + upper.2.weak)]
  medium.2[[i]]  <- N.trans.2[lower.2.medium:(lower.2.medium + upper.2.medium)]
  strong.2[[i]]  <- N.trans.2[lower.2.strong:(lower.2.strong + upper.2.strong)]
  
  none.1[[i]]    <- N.trans.1[lower.1.none:(lower.1.none + upper.1.none)]
  weak.1[[i]]    <- N.trans.1[lower.1.weak:(lower.1.weak + upper.1.weak)]
  medium.1[[i]]  <- N.trans.1[lower.1.medium:(lower.1.medium + upper.1.medium)]
  strong.1[[i]]  <- N.trans.1[lower.1.strong:(lower.1.strong + upper.1.strong)]
  
  none.50[[i]]    <- N.trans.50[lower.50.none:(lower.50.none + upper.50.none)]
  weak.50[[i]]   <- N.trans.50[lower.50.weak:(lower.50.weak + upper.50.weak)]
  medium.50[[i]]  <- N.trans.50[lower.50.medium:(lower.50.medium + upper.50.medium)]
  strong.50[[i]] <- N.trans.50[lower.50.strong:(lower.50.strong + upper.50.strong)]
  
  none.25[[i]]    <- N.trans.25[lower.25.none:(lower.25.none + upper.25.none)]
  weak.25[[i]]   <- N.trans.25[lower.25.weak:(lower.25.weak + upper.25.weak)]
  medium.25[[i]]  <- N.trans.25[lower.25.medium:(lower.25.medium + upper.25.medium)]
  strong.25[[i]] <- N.trans.25[lower.25.strong:(lower.25.strong + upper.25.strong)]
  
}


## 1/4 G Sampling #####

# Strong forcing
m <- length(strong.25[[i]]) - 7

res.strong.25.sd <- matrix(ncol = m, nrow = reps)
res.strong.25.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(strong.25[[i]][1:(length(strong.25[[i]])-j)], detrending = "gaussian")
    res.strong.25.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.strong.25.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Medium forcing
m <- length(medium.25[[1]]) - 7

res.medium.25.sd <- matrix(ncol = m, nrow = reps)
res.medium.25.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(medium.25[[i]][1:(length(medium.25[[i]])-j)], detrending = "gaussian")
    res.medium.25.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.medium.25.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Weak forcing
m <- length(weak.25[[i]]) - 7

res.weak.25.sd <- matrix(ncol = m, nrow = reps)
res.weak.25.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(weak.25[[i]][1:(length(weak.25[[i]])-j)], detrending = "gaussian")
    res.weak.25.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.weak.25.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

## 1/2 G Sampling #####

# Strong forcing
m <- length(strong.50[[i]]) - 7

res.strong.50.sd <- matrix(ncol = m, nrow = reps)
res.strong.50.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(strong.50[[i]][1:(length(strong.50[[i]])-j)], detrending = "gaussian")
    res.strong.50.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.strong.50.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Medium forcing
m <- length(medium.50[[i]]) - 5

res.medium.50.sd <- matrix(ncol = m, nrow = reps)
res.medium.50.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(medium.50[[i]][1:(length(medium.50[[i]])-j)], detrending = "gaussian")
    res.medium.50.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.medium.50.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Weak forcing
m <- length(weak.50[[i]]) - 5

res.weak.50.sd <- matrix(ncol = m, nrow = reps)
res.weak.50.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(weak.50[[i]][1:(length(weak.50[[i]])-j)], detrending = "gaussian")
    res.weak.50.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.weak.50.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

## 1 G Sampling #####

# Strong forcing
m <- length(strong.1[[i]]) - 7

res.strong.1.sd <- matrix(ncol = m, nrow = reps)
res.strong.1.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(strong.1[[i]][1:(length(strong.1[[i]])-j)], detrending = "gaussian")
    res.strong.1.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.strong.1.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Medium forcing
m <- length(medium.1[[i]]) - 7

res.medium.1.sd <- matrix(ncol = m, nrow = reps)
res.medium.1.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(medium.1[[i]][1:(length(medium.1[[i]])-j)], detrending = "gaussian")
    res.medium.1.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.medium.1.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Weak forcing
m <- length(weak.1[[i]]) - 7

res.weak.1.sd <- matrix(ncol = m, nrow = reps)
res.weak.1.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(weak.1[[i]][1:(length(weak.1[[i]])-j)], detrending = "gaussian")
    res.weak.1.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.weak.1.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}


## 2 G Sampling #####

# Strong forcing
m <- length(strong.2[[i]]) - 7

res.strong.2.sd <- matrix(ncol = m, nrow = reps)
res.strong.2.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(strong.2[[i]][1:(length(strong.2[[i]])-j)], detrending = "gaussian")
    res.strong.2.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.strong.2.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Medium forcing
m <- length(medium.2[[i]]) - 7

res.medium.2.sd <- matrix(ncol = m, nrow = reps)
res.medium.2.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(medium.1[[i]][1:(length(medium.1[[i]])-j)], detrending = "gaussian")
    res.medium.2.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.medium.2.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

# Weak forcing
m <- length(weak.2[[i]]) - 7

res.weak.2.sd <- matrix(ncol = m, nrow = reps)
res.weak.2.au <- matrix(ncol = m, nrow = reps)

for(i in 1:reps){
  for(j in 0:(m-1)){
    
    # Perform analysis and store value
    res <- EWS_analysis(weak.2[[i]][1:(length(weak.2[[i]])-j)], detrending = "gaussian")
    res.weak.2.sd[i,(j+1)] <- res[1]  # Standard Deviation
    res.weak.2.au[i,(j+1)] <- res[2]  # Autocorrelation
    
  }
}

## Save Results #####

save(res.strong.25.au, res.strong.50.au, res.strong.1.au, res.strong.2.au,
     res.strong.25.sd, res.strong.50.sd, res.strong.1.sd, res.strong.2.sd,
     res.medium.25.au, res.medium.50.au, res.medium.1.au, res.medium.2.au,
     res.medium.25.sd, res.medium.50.sd, res.medium.1.sd, res.medium.2.sd,
     res.weak.25.au, res.weak.50.au, res.weak.1.au, res.weak.2.au,
     res.weak.25.sd, res.weak.50.sd, res.weak.1.sd, res.weak.2.sd,
     file = "analyzedTrans.RData")
