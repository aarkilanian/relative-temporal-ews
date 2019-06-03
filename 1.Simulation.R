# 1.
# This script runs population simulations and generates the abundance data
# on which EWS analyses are performed.
#
# Code developed by Alex Arkilanian
#----------------------------------------------------------------------------

## Starting stuff #####

# Clean
rm(list = ls())

# Packages
require(ggplot2)

## Inital parameters of model #####

N.init <- 100   # Initial population size
r <- 0.5       # Population growth rate
K <- 100        # Carrying capacity
h <- 1          # Half-saturation constant

c.list.fold <- c(0, 0.015, 0.025, 0.035)            # Forcing parameter for fold
c.list.trans <- c(0, 0.00085, 0.0010, 0.0015)          # Forcing parameter for trans

# Replace these with fractions of generation time
dt.list <- c(0.3, (1.4 * 0.25), (1.4 * 0.5), (1.4 * 1), (1.4 * 2)) # List of resolutions

tmax <- 1000   # Max length of time series
tmin <- 0   # Min length of time series

reps <- 100    # Number of repetitions to run, script can handle up to 100

## Initial dataframe and logical vectors #####

# Number of data points within each time step for each forcing
p.list <- c((-tmin + tmax)/dt.list[1] + 1, # real
            floor((-tmin + tmax)/dt.list[2] + 1),
            floor((-tmin + tmax)/dt.list[3] + 1),
            floor((-tmin + tmax)/dt.list[4] + 1),
            floor((-tmin + tmax)/dt.list[5] + 1))

# Forcing in order for fold models
c.fold <- c(rep(c.list.fold, each = p.list[1]))

# Forcing in order for trans models
c.trans <- c(rep(c.list.trans, each = p.list[1]))

# Time values
time <- c(rep(seq(tmin, tmax, dt.list[1]), length(c.list.fold)))

# Logical vectors marking when to reinitialize population size
# Determined as 
reInit <- time == tmin
reInit.2 <- rep(c(TRUE, rep(FALSE, (p.list[5] - 1))), 9)
reInit.1 <- rep(c(TRUE, rep(FALSE, (p.list[4] - 1))), 9)
reInit.50 <- rep(c(TRUE, rep(FALSE, (p.list[3] - 1))), 9)
reInit.25 <- rep(c(TRUE, rep(FALSE, (p.list[2] - 1))), 9)

# When to sample from high-res to low-res populations
# Determined as when time point is multiple of resolution step
TF.2 <- (time %% dt.list[5]) == 0
TF.1 <- (time %% dt.list[4]) == 0
TF.50 <- (time %% dt.list[3]) == 0
TF.25 <- (time %% dt.list[2]) == 0

## Running model at dt = dt.list[1] and subsetting to low resolution #####

# Loop initializations
dt <- dt.list[1]

Ns.fold.2 <- vector("list", reps)
Ns.trans.2 <- vector("list", reps)

Ns.fold.1 <- vector("list", reps)
Ns.trans.1 <- vector("list", reps)

Ns.fold.50 <- vector("list", reps)
Ns.trans.50 <- vector("list", reps)

Ns.fold.25 <- vector("list", reps)
Ns.trans.25 <- vector("list", reps)

Ns.fold <- vector("list", reps)
Ns.trans <- vector("list", reps)

# Main loop

for (j in 1:(reps)) {
  
  a <- b <- c <- d <- 1
  
  fold.dead <- FALSE   
  trans.dead <- FALSE
  
  ### Initialize containers
  # Real-time containers
  N.fold <- numeric(p.list[1]*4)
  N.trans <- numeric(p.list[1]*4)
  
  # Quarter generation time containers
  N.fold.25 <- numeric(p.list[2]*4)
  N.trans.25 <- numeric(p.list[2]*4)
  
  # Half generation time containers
  N.fold.50 <- numeric(p.list[3]*4)
  N.trans.50 <- numeric(p.list[3]*4)
  
  # Full generation time containers
  N.fold.1 <- numeric(p.list[4]*4)
  N.trans.1 <- numeric(p.list[4]*4)
  
  # Double generation time containers
  N.fold.2 <- numeric(p.list[5]*4)
  N.trans.2 <- numeric(p.list[5]*4)
  
  for (i in 1:length(time)) {
    
    # Populations initilization and reinitialization
    if (reInit[i]) {
      N.fold[i] <- N.init
      N.trans[i] <- N.init
      fold.dead <- FALSE
      trans.dead <- FALSE
    } else {
      # Retrieve values here rather than retrieving multiple times in loop
      curTime <- time[i]        
      curHarv.fold <- c.fold[i]
      curHarv.trans <- c.trans[i]
      
      ## Fold bifurcation
      if(fold.dead){} # If pop extinct don't do anything
      else {
        N.past <- N.fold[i-1] 
        growth <- N.past * r * (1 - (N.past / K)) * dt 
        noise  <- rnorm(1,0,1.5) * dt
        harvest <- (curHarv.fold * curTime * N.past^2 / (N.past^2 + h^2)) * dt
        
        N.cur <- N.past + growth + noise - harvest
        
        # If dead, keep it dead and prevent going into negative pop size
        if (N.cur < 0) { 
          fold.dead <- TRUE
          N.cur <- 0
        }
        
        # Store N
        N.fold[i]      <- N.cur
      }
      
      ## Transcritical bifurcation
      if(trans.dead){} # If pop extinct don't do anything
      else {
        N.past <- N.trans[i-1]
        growth <- N.past * r * (1 - (N.past / K)) * dt
        noise <- rnorm(1,0,1.5) * dt
        harvest <- curHarv.trans * curTime * N.past * dt
        
        N.cur  <- N.past + growth + noise - harvest
        
        # If dead, keep it dead and prevent going into negative pop size 
        if (N.cur < 0) {
          trans.dead <- TRUE
          N.cur <- 0
        }
        
        # Store N
        N.trans[i]      <- N.cur
      }
    }
    
    # Population subsetting
    if(TF.2[i]) {
      # Store N
      N.fold.2[a]  <- N.fold[i]
      N.trans.2[a] <- N.trans[i]
      a <- a + 1
    }
    if(TF.1[i])    {
      # Store N
      N.fold.1[b]  <- N.fold[i]
      N.trans.1[b] <- N.trans[i]
      b <- b + 1
    }
    if(TF.50[i]) {
      # Store N
      N.fold.50[c]  <- N.fold[i]
      N.trans.50[c] <- N.trans[i]
      c <- c + 1
    }
    if(TF.25[i])    {
      # Store N
      N.fold.25[d]  <- N.fold[i]
      N.trans.25[d] <- N.trans[i]
      d <- d + 1
    }
    
  } 
  
  # N results
  Ns.fold.2[[j]]    <- N.fold.2
  Ns.fold.1[[j]]    <- N.fold.1
  Ns.fold.50[[j]]   <- N.fold.50
  Ns.fold.25[[j]]   <- N.fold.25
  
  Ns.trans.2[[j]]   <- N.trans.2
  Ns.trans.1[[j]]   <- N.trans.1
  Ns.trans.50[[j]]  <- N.trans.50
  Ns.trans.25[[j]]  <- N.trans.25
  
  Ns.fold[[j]]      <- N.fold
  Ns.trans[[j]]     <- N.trans
  
}

# Clean up

save(  Ns.fold.2, Ns.fold.1, Ns.fold.50, Ns.fold.25,
       Ns.trans.2, Ns.trans.1, Ns.trans.50, Ns.trans.25,
       Ns.fold, Ns.trans,
       file = "AllPop.RData")
