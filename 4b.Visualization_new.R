##### Load functions and packages #####

rm(list = ls())
require(ggplot2)
require(gridExtra)
require(reshape2)
require(plyr)
require(tidyr)

namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}

se <- function(x){
  sqrt(stats::var(x)/length(x))
}

false.neg <- function(x){
  (length(which(x <= 0))) / length(x)
}

##### Process fold simulation data #####

load("analyzedFold.RData")

matrix.list <- namedList(res.weak.2.sd, res.weak.2.au, res.weak.1.sd, res.weak.1.au, 
                         res.weak.50.sd, res.weak.50.au, res.weak.25.sd, res.weak.25.au,
                         res.medium.2.sd, res.medium.2.au, res.medium.1.sd, res.medium.1.au, 
                         res.medium.50.sd, res.medium.50.au, res.medium.25.sd, res.medium.25.au,
                         res.strong.2.sd, res.strong.2.au, res.strong.1.sd, res.strong.1.au, 
                         res.strong.50.sd, res.strong.50.au, res.strong.25.sd, res.strong.25.au)

biggy.fold <- data.frame(strength = factor(),
                         resolution = factor(),
                         length = numeric(),
                         metric = factor(),
                         tau = numeric())

levels(biggy.fold$strength) <- c("weak","medium","strong")
levels(biggy.fold$resolution) <- c("2", "1", "25", "50")
levels(biggy.fold$metric) <- c("sd","ar(1)")

biggy.fold[1,] <- c("weak","2",1,"sd", 1.0)

j <- 1
for(i in matrix.list){
  taus <- numeric()
  for(h in 1:length(i[1,])){taus <- append(taus, i[,h])}
  biggy.fold <- rbind(biggy.fold,
                      data.frame(
                        strength = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[2], times = length(i)),
                        resolution = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[3], times = length(i)),
                        length = rep((length(i[1,]) + 3):4, each = 100),
                        metric = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[4], times = length(i)),
                        tau = taus
                      )
  )
  j <- j + 1
}
biggy.fold <- biggy.fold[-1,]

biggy.fold$tau <- as.numeric(biggy.fold$tau)
biggy.fold$length <- as.numeric(biggy.fold$length)

biggy.fold$length.gen.tmp <- ifelse(biggy.fold$resolution %in% "2", 2,
                                    ifelse(biggy.fold$resolution %in% "1", 1,
                                           ifelse(biggy.fold$resolution %in% "50", 0.5, 0.25)))

biggy.fold$length.gen <- biggy.fold$length * biggy.fold$length.gen.tmp

# 
# biggy.summary.fold <- aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, mean)
# 
# biggy.summary.fold$se <- unlist(aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, se)[5])
# 
# biggy.summary.fold$top <- unlist(aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, quantile, probs = 0.75)[5])
# 
# biggy.summary.fold$bot <- unlist(aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, quantile, probs = 0.25)[5])
# 
# biggy.summary.fold$false.neg <- unlist(aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, false.neg)[5])
# 

for(i in which(biggy.fold$strength %in% "none")){
  biggy.fold$false.neg[i] <- (1 - biggy.fold$false.neg[i])
}


##### Process trans simulation data #####

load("analyzedTrans.RData")

matrix.list <- namedList(res.weak.2.sd, res.weak.2.au, res.weak.1.sd, res.weak.1.au, 
                         res.weak.50.sd, res.weak.50.au, res.weak.25.sd, res.weak.25.au,
                         res.medium.2.sd, res.medium.2.au, res.medium.1.sd, res.medium.1.au, 
                         res.medium.50.sd, res.medium.50.au, res.medium.25.sd, res.medium.25.au,
                         res.strong.2.sd, res.strong.2.au, res.strong.1.sd, res.strong.1.au, 
                         res.strong.50.sd, res.strong.50.au, res.strong.25.sd, res.strong.25.au)

biggy.trans <- data.frame(strength = factor(),
                          resolution = factor(),
                          length = numeric(),
                          metric = factor(),
                          tau = numeric())

levels(biggy.trans$strength) <- c("weak","medium","strong")
levels(biggy.trans$resolution) <- c("2", "1", "25", "50")
levels(biggy.trans$metric) <- c("sd","ar(1)")

biggy.trans[1,] <- c("weak","2",1,"sd", 1.0)

j <- 1
for(i in matrix.list){
  taus <- numeric()
  for(h in 1:length(i[1,])){taus <- append(taus, i[,h])}
  biggy.trans <- rbind(biggy.trans,
                       data.frame(
                         strength = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[2], times = length(i)),
                         resolution = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[3], times = length(i)),
                         length = rep((length(i[1,]) + 3):4, each = 100),
                         metric = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[4], times = length(i)),
                         tau = taus
                       )
  )
  j <- j + 1
}
biggy.trans <- biggy.trans[-1,]

biggy.trans$tau <- as.numeric(biggy.trans$tau)
biggy.trans$length <- as.numeric(biggy.trans$length)

biggy.trans$length.gen.tmp <- ifelse(biggy.trans$resolution %in% "2", 2,
                                     ifelse(biggy.trans$resolution %in% "1", 1,
                                            ifelse(biggy.trans$resolution %in% "50", 0.5, 0.25)))

biggy.trans$length.gen <- biggy.trans$length * biggy.trans$length.gen.tmp

# biggy.summary.trans <- aggregate(data = biggy.trans, tau ~ strength * resolution * length.gen * metric, mean)
# 
# biggy.summary.trans$se <- unlist(aggregate(data = biggy.trans, tau ~ strength * resolution * length.gen * metric, se)[5])
# 
# biggy.summary.trans$top <- unlist(aggregate(data = biggy.trans, tau ~ strength * resolution * length.gen * metric, quantile, probs = 0.75)[5])
# 
# biggy.summary.trans$bot <- unlist(aggregate(data = biggy.trans, 
#                                             tau ~ strength * resolution * length.gen * metric, quantile, probs = 0.25)[5])
# 
# biggy.summary.trans$false.neg <- unlist(aggregate(data = biggy.trans, tau ~ strength * resolution * length.gen * metric, false.neg)[5])


##### Stitch datasets together #####

biggy.fold$method <- "fold"
biggy.trans$method <- "trans"
merged <- rbind(biggy.fold, biggy.trans)
merged$method <- as.factor(merged$method)

merged <- subset(merged, strength %in% c("weak", "medium", "strong"))
merged$metric <- revalue(merged$metric, c("sd"="sd", "au"="ar(1)"))
merged$strength <- revalue(merged$strength, c("weak"="slow", "medium"="moderate", "strong"="fast"))
merged$length.gen.bin <- cut(merged$length.gen, breaks = c(0,5,10,15,20,25,30,35,40,45),
                             labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44"))
merged$length.gen.bin <- with(merged, factor(length.gen.bin, levels = rev(levels(length.gen.bin))))

merged.2 <- subset(merged, !(merged$strength %in% "none"))

merged.2$resolution <- revalue(merged$resolution, c("2"="2", "1"="1", "50"="0.5", "25"="0.25"))
merged.2$resolution <- factor(merged.2$resolution, levels = c("2", "1", "0.5", "0.25"))


##### Figure 1: resolution of sampling in simulations #####

# Simulation resolution analysis

merged <- subset(merged, !(merged$strength %in% "none"))



merged$resolution <- revalue(merged$resolution, c("2"="2", "1"="1", "50"="0.5", "25"="0.25"))
merged$resolution <- factor(merged$resolution, levels = c("2", "1", "0.5", "0.25"))
merged$resolution <- factor(merged$resolution, levels(merged$resolution)[c(4, 3, 2, 1)])


#slow forcing and fold
slow.forcing.fold.res2<-subset(merged, strength=="slow" & 
                                 resolution == 2 & length ==12 & method =="fold")
slow.forcing.fold.res1<-subset(merged, strength=="slow" & 
                                 resolution == 1 & length ==26 & method =="fold")

slow.forcing.fold.res0.5<-subset(merged, strength=="slow" & 
                                   resolution == 0.5 & length ==55 & method =="fold")

slow.forcing.fold.res0.25<-subset(merged, strength=="slow" & 
                                    resolution == 0.25 & length ==113 & method =="fold")


#fast forcing and fold
fast.forcing.fold.res2<-subset(merged, strength=="fast" & 
                                 resolution == 2 & length == 6  & method =="fold")
fast.forcing.fold.res1<-subset(merged, strength=="fast" & 
                                 resolution == 1 & length ==12 & method =="fold")

fast.forcing.fold.res0.5<-subset(merged, strength=="fast" & 
                                   resolution == 0.5 & length ==24 & method =="fold")

fast.forcing.fold.res0.25<-subset(merged, strength=="fast" & 
                                    resolution == 0.25 & length ==49 & method =="fold")

# moderate forcing and fold
mod.forcing.fold.res2<-subset(merged, strength=="moderate" & 
                                resolution == 2 & length == 8 & method =="fold")
mod.forcing.fold.res1<-subset(merged, strength=="moderate" & 
                                resolution == 1 & length == 18 & method =="fold")

mod.forcing.fold.res0.5<-subset(merged, strength=="moderate" & 
                                  resolution == 0.5 & length ==37 & method =="fold")

mod.forcing.fold.res0.25<-subset(merged, strength=="moderate" & 
                                   resolution == 0.25 & length ==75 & method =="fold")
fold.data<-rbind(mod.forcing.fold.res2,mod.forcing.fold.res1,mod.forcing.fold.res0.5,mod.forcing.fold.res0.25,
                 fast.forcing.fold.res2,fast.forcing.fold.res1,fast.forcing.fold.res0.5,fast.forcing.fold.res0.25,
                 slow.forcing.fold.res2,slow.forcing.fold.res1,slow.forcing.fold.res0.5,slow.forcing.fold.res0.25)

# transcritical data resolution

#slow forcing and trans
slow.forcing.trans.res2<-subset(merged, strength=="slow" & 
                                  resolution == 2 & length ==9 & method =="trans")
slow.forcing.trans.res1<-subset(merged, strength=="slow" & 
                                  resolution == 1 & length ==20 & method =="trans")

slow.forcing.trans.res0.5<-subset(merged, strength=="slow" & 
                                    resolution == 0.5 & length ==44 & method =="trans")

slow.forcing.trans.res0.25<-subset(merged, strength=="slow" & 
                                     resolution == 0.25 & length ==86 & method =="trans")


#fast forcing and trans
fast.forcing.trans.res2<-subset(merged, strength=="fast" & 
                                  resolution == 2 & length ==4 & method =="trans")
fast.forcing.trans.res1<-subset(merged, strength=="fast" & 
                                  resolution == 1 & length ==9 & method =="trans")

fast.forcing.trans.res0.5<-subset(merged, strength=="fast" & 
                                    resolution == 0.5 & length ==18 & method =="trans")

fast.forcing.trans.res0.25<-subset(merged, strength=="fast" & 
                                     resolution == 0.25 & length == 38 & method =="trans")

# moderate forcing and trans
mod.forcing.trans.res2<-subset(merged, strength=="moderate" & 
                                 resolution == 2 & length == 5 & method =="trans")
mod.forcing.trans.res1<-subset(merged, strength=="moderate" & 
                                 resolution == 1 & length == 12 & method =="trans")

mod.forcing.trans.res0.5<-subset(merged, strength=="moderate" & 
                                   resolution == 0.5 & length ==28 & method =="trans")

mod.forcing.trans.res0.25<-subset(merged, strength=="moderate" & 
                                    resolution == 0.25 & length ==53 & method =="trans")

trans.data<-rbind(mod.forcing.trans.res2,mod.forcing.trans.res1,mod.forcing.trans.res0.5,mod.forcing.trans.res0.25,
                  fast.forcing.trans.res2,fast.forcing.trans.res1,fast.forcing.trans.res0.5,fast.forcing.trans.res0.25,
                  slow.forcing.trans.res2,slow.forcing.trans.res1,slow.forcing.trans.res0.5,slow.forcing.trans.res0.25)



#combined trans + fold data

data.resolution<-rbind(trans.data,fold.data)

library(wesanderson)
best_color_paletter<- c(wes_palettes$Darjeeling1, wes_palettes$Rushmore1)

sim.res <- ggplot(data = data.resolution, aes(x = resolution, y = tau, fill = method:metric)) +
  geom_hline(yintercept = 0, col = 'gray', size = 0.5) +
  geom_point(aes(col = method:metric), position = position_jitterdodge(dodge.width = 0.7), alpha = 0.4, size = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.2) +
  theme_light() +
  #theme(strip.background = element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=8)) +
  #theme(legend.position = "none", text = element_text(size=8)) +
  scale_fill_manual(values=best_color_paletter)+
  scale_y_continuous(name = expression(paste("Kendall's ", tau)), limits = c(-1,1)) +
  scale_x_discrete(name = "Resolution of timeseries (/generations)") +
  facet_wrap(~ metric * strength)

sim.res