# 3.
# This script sets up microcosm data for analysis of the effect of length of 
# sampling on EWS across all sampling resolutions for strongest and weakest
# forcing levels
#
# Code developed by Alex Arkilanian
#-----------------------------------------------------------------------

## Starting Stuff #####

# Clean
rm(list = ls())
 
# Packages
require(ggplot2)
require(earlywarnings)
require(reshape2)
require(mgcv)

# Functions

interp<-function(days, obs){
  if(length(na.omit(obs))<length(obs)){
    return(approx(days, obs, n = length(obs), xout=days, method = "linear")$y)}else{
      return(obs)}
}

# Take a timeseries and return kendall's tau of ews
# Timeseries is abundance vector
# Trait is trait vector (body size in this case)
# max.length is the max length that vectors should have
genericEWSTrait <- function (days, timeseries, trait, max.length = length(timeseries), winsize = 50, detrending = "gaussian", bandwidth = NULL, span = NULL, 
                        degree = NULL, logtransform = FALSE, interpolate = FALSE, 
                        AR_n = FALSE, powerspectrum = FALSE) 
{
  AR <- integer()
  SD <- integer()
  AR.trait <- integer()
  SD.trait <- integer()
  SD.AR.trait <- integer()
  trait.only <- integer()
  length <- integer()
  
  timeseries <- interp(days, timeseries)
  trait <- interp(days, trait)
  timeindex <- 1:max.length
  
  for(j in 1:(max.length-5)){
    
  # Perform appropriate detrending based on what is given to function
  bw <- round(bw.nrd0(timeindex))
  
  smoothed <- ksmooth(timeindex, timeseries, kernel = "normal", bandwidth = bw, 
                  range.x = range(timeindex), x.points = timeindex)
  diff <- timeseries - smoothed$y
  smoothSeries <- smoothed$y
  
  # Rolling window count
  win.length <- round(length(timeseries) * winsize/100) # Moving window size
  win.number <- length(smoothSeries) - win.length + 1   # Numbe of rolling windows
  low <- 6
  high <- win.number
  
  res <- matrix(data = NA, nrow = win.length, ncol = win.number)
  for (i in 1:win.number) { 
    tmp <- diff[i:(i + win.length - 1)]
    res[, i] <- tmp 
    }
  
  # Rolling window trait
  win.length.t <- round(length(timeseries) * winsize/100)
  win.number.t <- length(smoothSeries) - win.length.t + 1
  low.t <- 6
  high.t <- win.number.t
  
  res.t <- matrix(data = NA, nrow = win.length.t, ncol = win.number.t)
  for (i in 1:win.number.t) { 
    tmp <- trait[i:(i + win.length.t - 1)]
    res.t[, i] <- tmp }
  
  # Empty containers to hold ews data
  nSD <- numeric()
  nARR <- numeric()
  
  # Calculate ews over res dataset
  nSD <- apply(res, 2, sd, na.rm = TRUE)
  
  for (i in 1:ncol(res)) {
    nARR[i] <- ar.ols(res[, i], aic = FALSE, order.max = 1, dmean = FALSE,
                  intercept = FALSE)$ar
  }
  
  trait.rolled <- apply(res.t, 2, mean, na.rm = TRUE)
  
  # Z-transformations
  trait.z  <- ((trait.rolled-mean(trait.rolled))/(sd(trait.rolled, na.rm = TRUE))) # body size
  ar.z     <- (nARR-mean(nARR))/sd(nARR, na.rm = TRUE)              # autocorrelation
  sd.z     <- (nSD-mean(nSD))/sd(nSD, na.rm = TRUE)              # standard deviation
  
    # Correlation test statistic (kendall's tau) and confidence intervals - abundance based data
    time   <- 1:win.number;
    
    tau.ar <- cor.test(ar.z, time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate 
    
    tau.sd <- cor.test(sd.z, time, method = "kendall", conf.level = 0.95, alternative ="two.sided")$estimate
    
    # Correlation test statistic (kendall's tau) and confidence intervals - including trait data
    tau.trait.ar    <- cor.test((trait.z+ar.z), time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate
    
    tau.trait.sd    <- cor.test((trait.z+sd.z), time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate
    
    tau.trait.sd.ar <- cor.test((trait.z+sd.z+ar.z), time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate
    
    # Correlation test statistic (kendall's tau) and confidence intervals - only trait data
    tau.trait.only <- cor.test(trait.z, time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate
    
    # Store values
    AR[j] <- tau.ar
    SD[j] <- tau.sd
    AR.trait[j] <- tau.trait.ar
    SD.trait[j] <- tau.trait.sd
    SD.AR.trait[j] <- tau.trait.sd.ar
    trait.only[j] <- tau.trait.only
    
    length[j] <- length(timeindex)
  
    # Remove one data point
    timeseries <- timeseries[-1]
    trait <- trait[-1]
    timeindex <- 1:(max.length - j)
  }
  out <- data.frame(length, AR, SD, AR.trait, SD.trait, SD.AR.trait, trait.only)
  return(out)
}

genericEWSTrait.half <- function (days, timeseries, trait, max.length = length(timeseries), winsize = 50, detrending = "gaussian", bandwidth = NULL, span = NULL, 
                                  degree = NULL, logtransform = FALSE, interpolate = FALSE, 
                                  AR_n = FALSE, powerspectrum = FALSE) 
{
  
  AR <- integer()
  SD <- integer()
  AR.trait <- integer()
  SD.trait <- integer()
  SD.AR.trait <- integer()
  trait.only <- integer()
  length <- integer()
  
  # Interpolation
  timeseries <- interp(days, timeseries)
  trait <- interp(days, trait)
  
  # Subsample, take every second element
  timeseries <- timeseries[seq(1, length(timeseries), 2)] 
  trait <- trait[seq(1, length(trait), 2)] 
  
  # Adjust max.length

  timeindex <- 1:max.length
  
  for(j in 1:(max.length-5)){
    
    # Perform appropriate detrending based on what is given to function
    bw <- round(bw.nrd0(timeindex))
    
    smoothed <- ksmooth(timeindex, timeseries, kernel = "normal", bandwidth = bw, 
                        range.x = range(timeindex), x.points = timeindex)
    diff <- timeseries - smoothed$y
    smoothSeries <- smoothed$y
    
    # Rolling window count
    win.length <- round(length(timeseries) * winsize/100) # Moving window size
    win.number <- length(smoothSeries) - win.length + 1   # Numbe of rolling windows
    low <- 6
    high <- win.number
    
    res <- matrix(data = NA, nrow = win.length, ncol = win.number)
    for (i in 1:win.number) { 
      tmp <- diff[i:(i + win.length - 1)]
      res[, i] <- tmp 
    }
    
    # Rolling window trait
    win.length.t <- round(length(timeseries) * winsize/100)
    win.number.t <- length(smoothSeries) - win.length.t + 1
    low.t <- 6
    high.t <- win.number.t
    
    res.t <- matrix(data = NA, nrow = win.length.t, ncol = win.number.t)
    for (i in 1:win.number.t) { 
      tmp <- trait[i:(i + win.length.t - 1)]
      res.t[, i] <- tmp }
    
    # Empty containers to hold ews data
    nSD <- numeric()
    nARR <- numeric()
    
    # Calculate ews over res dataset
    nSD <- apply(res, 2, sd, na.rm = TRUE)
    
    for (i in 1:ncol(res)) {
      nARR[i] <- ar.ols(res[, i], aic = FALSE, order.max = 1, dmean = FALSE,
                        intercept = FALSE)$ar
    }
    
    trait.rolled <- apply(res.t, 2, mean, na.rm = TRUE)
    
    # Z-transformations
    trait.z  <- ((trait.rolled-mean(trait.rolled))/(sd(trait.rolled, na.rm = TRUE))) # body size
    ar.z     <- (nARR-mean(nARR))/sd(nARR, na.rm = TRUE)              # autocorrelation
    sd.z     <- (nSD-mean(nSD))/sd(nSD, na.rm = TRUE)              # standard deviation
    
    # Correlation test statistic (kendall's tau) and confidence intervals - abundance based data
    time   <- 1:win.number;
    
    tau.ar <- cor.test(ar.z, time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate 
    
    tau.sd <- cor.test(sd.z, time, method = "kendall", conf.level = 0.95, alternative ="two.sided")$estimate
    
    # Correlation test statistic (kendall's tau) and confidence intervals - including trait data
    tau.trait.ar    <- cor.test((trait.z+ar.z), time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate
    
    tau.trait.sd    <- cor.test((trait.z+sd.z), time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate
    
    tau.trait.sd.ar <- cor.test((trait.z+sd.z+ar.z), time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate
    
    # Correlation test statistic (kendall's tau) and confidence intervals - only trait data
    tau.trait.only <- cor.test(trait.z, time, method = "kendall", conf.level = 0.95, alternative = "two.sided")$estimate
    
    # Store values
    AR[j] <- tau.ar
    SD[j] <- tau.sd
    AR.trait[j] <- tau.trait.ar
    SD.trait[j] <- tau.trait.sd
    SD.AR.trait[j] <- tau.trait.sd.ar
    trait.only[j] <- tau.trait.only
    
    length[j] <- length(timeindex)
    
    # Remove one data point
    timeseries <- timeseries[-1]
    trait <- trait[-1]
    timeindex <- 1:(max.length - j)
  }
  out <- data.frame(length, AR, SD, AR.trait, SD.trait, SD.AR.trait, trait.only)
  return(out)
  
}

# Data
exp.data <- read.csv("microcosm_data.csv")
tipping.points <- read.csv("Tipping points.csv")[,-1]

## Data manipulations #####

# Reverse body size direction
exp.data$mean.size <- -1 * exp.data$mean.size

# Subset data
con.data     <- subset(exp.data, Treatment %in% "Constant")
con.count    <- dcast(data = con.data, Day ~ Tube, value.var = "Count")
con.size <- dcast(data = con.data, Day ~ Tube, value.var = "mean.size")

slow.data     <- subset(exp.data, Treatment %in% "Slow")
slow.count    <- dcast(data = slow.data, Day ~ Tube, value.var = "Count")
slow.size    <- dcast(data = slow.data, Day ~ Tube, value.var = "mean.size")

medium.data     <- subset(exp.data, Treatment %in% "Medium")
medium.count    <- dcast(data = medium.data, Day ~ Tube, value.var = "Count")
medium.size    <- dcast(data = medium.data, Day ~ Tube, value.var = "mean.size")

fast.data     <- subset(exp.data, Treatment %in% "Fast")
fast.count    <- dcast(data = fast.data, Day ~ Tube, value.var = "Count")
fast.size    <- dcast(data = fast.data, Day ~ Tube, value.var = "mean.size")

# Remove constant which are deteriorating
con.count <- con.count[,-11]
con.size <- con.size[,-11]

# Remove others which are not deteriorating
slow.count <- slow.count[1:25,c(1,2,3,4,5,8)]
slow.size <- slow.size[1:25,c(1,2,3,5,6,8)]

medium.count <- medium.count[1:23,c(1,2,5,6,10)]
medium.size <- medium.size[1:23,c(1,2,5,6,10)]

fast.count <- fast.count[1:19,c(1,2,3,5,7,8,9,10)]
fast.size <- fast.size[1:19,c(1,2,3,5,7,8,9,10)]

## EWS & length of sampling #####

# Constant
con.mess <- mapply(genericEWSTrait, con.count[1], con.count[2:11], con.size[2:11])

length <- numeric()
ar <- numeric()
sd <- numeric()
ar.trait <- numeric()
sd.trait <- numeric()
ar.sd.trait <- numeric()
trait.only <- numeric()
j <- 0

for (i in seq(1, 290, 29)){
  j <- j + 1
  mess <- unname(unlist(con.mess[,j]))
  length[i:(i+28)] <- mess[1:29]
  ar[i:(i+28)] <- mess[30:58]
  sd[i:(i+28)] <- mess[59:87]
  ar.trait[i:(i+28)] <- mess[87:115]
  sd.trait[i:(i+28)] <- mess[116:144]
  ar.sd.trait[i:(i+28)] <- mess[145:173]
  trait.only[i:(i+28)] <- mess[174:202]
}

constant <- data.frame(strength = "none", length, length.gen = length * 0.5, 
                       metric = append(rep("ar", length(ar)), rep("sd", length(sd))), 
                       no.trait = append(ar, sd), trait = append(ar.trait, sd.trait), 
                       trait.only = trait.only)

# Slow
slow.mess <- mapply(genericEWSTrait, slow.count[1], slow.count[2:6], slow.size[2:6])

length <- numeric()
ar <- numeric()
sd <- numeric()
ar.trait <- numeric()
sd.trait <- numeric()
ar.sd.trait <- numeric()
trait.only <- numeric()
j <- 0

for (i in seq(1, 100, 20)){
  j <- j + 1
  mess <- unname(unlist(slow.mess[,j]))
  length[i:(i+19)] <- mess[1:20]
  ar[i:(i+19)] <- mess[21:40]
  sd[i:(i+19)] <- mess[41:60]
  ar.trait[i:(i+19)] <- mess[61:80]
  sd.trait[i:(i+19)] <- mess[81:100]
  ar.sd.trait[i:(i+19)] <- mess[101:120]
  trait.only[i:(i+19)] <- mess[121:140]
}

slow <- data.frame(strength = "slow", length, length.gen = length * 0.5, 
                       metric = append(rep("ar", length(ar)), rep("sd", length(sd))), 
                       no.trait = append(ar, sd), trait = append(ar.trait, sd.trait), 
                       trait.only = trait.only)
# Medium
med.mess <- mapply(genericEWSTrait, medium.count[1], medium.count[2:5], medium.size[2:5])

length <- numeric()
ar <- numeric()
sd <- numeric()
ar.trait <- numeric()
sd.trait <- numeric()
ar.sd.trait <- numeric()
trait.only <- numeric()
j <- 0

for (i in seq(1, 72, 18)){
  j <- j + 1
  mess <- unname(unlist(med.mess[,j]))
  length[i:(i+17)] <- mess[1:18]
  ar[i:(i+17)] <- mess[19:36]
  sd[i:(i+17)] <- mess[37:54]
  ar.trait[i:(i+17)] <- mess[55:72]
  sd.trait[i:(i+17)] <- mess[73:90]
  ar.sd.trait[i:(i+17)] <- mess[91:108]
  trait.only[i:(i+17)] <- mess[109:126]
}

medium <- data.frame(strength = "medium", length, length.gen = length * 0.5, 
                       metric = append(rep("ar", length(ar)), rep("sd", length(sd))), 
                       no.trait = append(ar, sd), trait = append(ar.trait, sd.trait),
                       trait.only <- trait.only)
# Fast
fast.mess <- mapply(genericEWSTrait, fast.count[1], fast.count[2:8], fast.size[2:8])

length <- numeric()
ar <- numeric()
sd <- numeric()
ar.trait <- numeric()
sd.trait <- numeric()
ar.sd.trait <- numeric()
trait.only <- numeric()
j <- 0

for (i in seq(1, 98, 14)){
  j <- j + 1
  mess <- unname(unlist(fast.mess[,j]))
  length[i:(i+13)] <- mess[1:14]
  ar[i:(i+13)] <- mess[15:28]
  sd[i:(i+13)] <- mess[29:42]
  ar.trait[i:(i+13)] <- mess[43:56]
  sd.trait[i:(i+13)] <- mess[57:70]
  ar.sd.trait[i:(i+13)] <- mess[71:84]
  trait.only[i:(i+13)] <- mess[85:98]
}

fast <- data.frame(strength = "fast", length, length.gen = length * 0.5, 
                       metric = append(rep("ar", length(ar)), rep("sd", length(sd))), 
                       no.trait = append(ar, sd), trait = append(ar.trait, sd.trait),
                       trait.only = trait.only)

# Save data

save(constant, slow, medium, fast, file = "TraitData.RData")

## EWS & resolution of sampling #####

# Constant
con.mess <- mapply(genericEWSTrait.half, con.count[1], con.count[2:11], con.size[2:11])

length <- numeric()
ar <- numeric()
sd <- numeric()
ar.trait <- numeric()
sd.trait <- numeric()
ar.sd.trait <- numeric()
trait.only <- numeric()
j <- 0

for (i in seq(1, 120, 12)){
  j <- j + 1
  mess <- unname(unlist(con.mess[,j]))
  length[i:(i+11)] <- mess[1:12]
  ar[i:(i+11)] <- mess[13:24]
  sd[i:(i+11)] <- mess[25:36]
  ar.trait[i:(i+11)] <- mess[37:48]
  sd.trait[i:(i+11)] <- mess[49:60]
  ar.sd.trait[i:(i+11)] <- mess[61:72]
  trait.only[i:(i+11)] <- mess[73:84]
}

con.half <- data.frame(strength = "none", length, length.gen = length, 
                       metric = append(rep("ar", length(ar)), rep("sd", length(sd))), 
                       no.trait = append(ar, sd), trait = append(ar.trait, sd.trait),
                       trait.only = trait.only)
# Slow
slow.mess <- mapply(genericEWSTrait.half, slow.count[1], slow.count[2:6], slow.size[2:6])

length <- numeric()
ar <- numeric()
sd <- numeric()
ar.trait <- numeric()
sd.trait <- numeric()
ar.sd.trait <- numeric()
trait.only <- numeric()
j <- 0

for (i in seq(1, 40, 8)){
  j <- j + 1
  mess <- unname(unlist(slow.mess[,j]))
  length[i:(i+7)] <- mess[1:8]
  ar[i:(i+7)] <- mess[9:16]
  sd[i:(i+7)] <- mess[17:24]
  ar.trait[i:(i+7)] <- mess[25:32]
  sd.trait[i:(i+7)] <- mess[33:40]
  ar.sd.trait[i:(i+7)] <- mess[41:48]
  trait.only[i:(i+7)] <- mess[49:56]
}

slow.half <- data.frame(strength = "slow", length, length.gen = length, 
                       metric = append(rep("ar", length(ar)), rep("sd", length(sd))), 
                       no.trait = append(ar, sd), trait = append(ar.trait, sd.trait),
                       trait.only = trait.only)
# Medium
med.mess <- mapply(genericEWSTrait.half, medium.count[1], medium.count[2:5], medium.size[2:5])

length <- numeric()
ar <- numeric()
sd <- numeric()
ar.trait <- numeric()
sd.trait <- numeric()
ar.sd.trait <- numeric()
trait.only <- numeric()
j <- 0

for (i in seq(1, 28, 7)){
  j <- j + 1
  mess <- unname(unlist(med.mess[,j]))
  length[i:(i+6)] <- mess[1:7]
  ar[i:(i+6)] <- mess[8:14]
  sd[i:(i+6)] <- mess[15:21]
  ar.trait[i:(i+6)] <- mess[22:28]
  sd.trait[i:(i+6)] <- mess[29:35]
  ar.sd.trait[i:(i+6)] <- mess[36:42]
  trait.only[i:(i+6)] <- mess[43:49]
}

medium.half <- data.frame(strength = "medium", length, length.gen = length, 
                       metric = append(rep("ar", length(ar)), rep("sd", length(sd))), 
                       no.trait = append(ar, sd), trait = append(ar.trait, sd.trait),
                       trait.only = trait.only)
# Fast
fast.mess <- mapply(genericEWSTrait.half, fast.count[1], fast.count[2:8], fast.size[2:8])

length <- numeric()
ar <- numeric()
sd <- numeric()
ar.trait <- numeric()
sd.trait <- numeric()
ar.sd.trait <- numeric()
trait.only <- numeric()
j <- 0

for (i in seq(1, 35, 5)){
  j <- j + 1
  mess <- unname(unlist(fast.mess[,j]))
  length[i:(i+4)] <- mess[1:5]
  ar[i:(i+4)] <- mess[6:10]
  sd[i:(i+4)] <- mess[11:15]
  ar.trait[i:(i+4)] <- mess[16:20]
  sd.trait[i:(i+4)] <- mess[21:25]
  ar.sd.trait[i:(i+4)] <- mess[26:30]
  trait.only[i:(i+4)] <- mess[31:35]
}

fast.half <- data.frame(strength = "fast", length, length.gen = length, 
                       metric = append(rep("ar", length(ar)), rep("sd", length(sd))), 
                       no.trait = append(ar, sd), trait = append(ar.trait, sd.trait),
                       trait.only = trait.only)
# Save

save(con.half, slow.half, medium.half, fast.half, file = "TraitHalfData.RData")

