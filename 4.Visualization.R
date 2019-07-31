# 4.
# This script takes the analysis results from simulated and experimental
# populations, summarizes them, and visualizes them into final publication
# plots.
#
# Code developed by Alex Arkilanian
#--------------------------------------------------------------------------

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

levels(biggy.fold$strength) <- c("none","weak","medium","strong")
levels(biggy.fold$resolution) <- c("2", "1", "25", "50")
levels(biggy.fold$metric) <- c("sd","ar(1)")

biggy.fold[1,] <- c("none","2",1,"sd", 1.0)

j <- 1
for(i in matrix.list){
  taus <- numeric()
  for(h in 1:length(i[1,])){taus <- append(taus, i[,h])}
  biggy.fold <- rbind(biggy.fold,
                 data.frame(
                   strength = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[2], times = length(i)),
                   resolution = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[3], times = length(i)),
                   length = rep((length(i[1,]) + 5):6, each = 100),
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


biggy.summary.fold <- aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, mean)

biggy.summary.fold$se <- unlist(aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, se)[5])

biggy.summary.fold$top <- unlist(aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, quantile, probs = 0.75)[5])

biggy.summary.fold$bot <- unlist(aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, quantile, probs = 0.25)[5])

biggy.summary.fold$false.neg <- unlist(aggregate(data = biggy.fold, tau ~ strength * resolution * length.gen * metric, false.neg)[5])


for(i in which(biggy.summary.fold$strength %in% "none")){
  biggy.summary.fold$false.neg[i] <- (1 - biggy.summary.fold$false.neg[i])
}

rm(list=setdiff(ls(), c("biggy.fold", "namedList", "se", "false.neg", "biggy.summary.fold")))

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

levels(biggy.trans$strength) <- c("none","weak","medium","strong")
levels(biggy.trans$resolution) <- c("2", "1", "25", "50")
levels(biggy.trans$metric) <- c("sd","ar(1)")

biggy.trans[1,] <- c("none","2",1,"sd", 1.0)

j <- 1
for(i in matrix.list){
  taus <- numeric()
  for(h in 1:length(i[1,])){taus <- append(taus, i[,h])}
  biggy.trans <- rbind(biggy.trans,
                      data.frame(
                        strength = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[2], times = length(i)),
                        resolution = rep(unlist(strsplit(names(matrix.list)[j], '.', fixed = T))[3], times = length(i)),
                        length = rep((length(i[1,]) + 5):6, each = 100),
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

biggy.summary.trans <- aggregate(data = biggy.trans, tau ~ strength * resolution * length.gen * metric, mean)

biggy.summary.trans$se <- unlist(aggregate(data = biggy.trans, tau ~ strength * resolution * length.gen * metric, se)[5])

biggy.summary.trans$top <- unlist(aggregate(data = biggy.trans, tau ~ strength * resolution * length.gen * metric, quantile, probs = 0.75)[5])

biggy.summary.trans$bot <- unlist(aggregate(data = biggy.trans, tau ~ strength * resolution * length.gen * metric, quantile, probs = 0.25)[5])

biggy.summary.trans$false.neg <- unlist(aggregate(data = biggy.trans, tau ~ strength * resolution * length.gen * metric, false.neg)[5])


for(i in which(biggy.summary.trans$strength %in% "none")){
  biggy.summary.trans$false.neg[i] <- (1 - biggy.summary.trans$false.neg[i])
}


rm(list=setdiff(ls(), c("biggy.fold", "namedList", "biggy.trans", "se", "false.neg", "biggy.summary.fold", "biggy.summary.trans")))

##### Stitch datasets together #####

biggy.summary.fold$method <- "fold"
biggy.summary.trans$method <- "trans"
merged <- rbind(biggy.summary.fold, biggy.summary.trans)
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

##### Figure 2: Length of sampling in simulations #####

sim.false.binned <- ggplot(data = subset(merged, resolution %in% "25"), aes(x = as.numeric(length.gen.bin), y = false.neg*100, col = metric:method)) +
  theme_light() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major.y = element_blank()) +
  geom_smooth(data = subset(merged, resolution %in% "25" & method %in% "fold"), method = "lm", size = 0.3, linetype = "solid", se = FALSE) +
  geom_smooth(data = subset(merged, resolution %in% "25" & method %in% "trans"), method = "lm", size = 0.3, linetype = "dashed", se = FALSE) +
  geom_point(data = subset(merged, resolution %in% "25" & method %in% "fold"), position = position_dodge(width = 0.8), alpha = 0.3, shape = 3, size = 0.5) +
  geom_point(data = subset(merged, resolution %in% "25" & method %in% "trans"), position = position_dodge(width = 0.8), alpha = 0.3, size = 0.5) +
  scale_color_manual(values = c("#e00000", "#e00000", "#dbaf00", "#dbaf00")) +
  facet_wrap(~ strength, nrow = 1, scales = "free_x") +
  scale_y_continuous("Rate of false negatives") +
  scale_x_continuous("Length of sampling timeseries (generations)",
                     breaks = c(1,2,3,4,5,6,7,8,9),
                     labels = c("40-44", "35-39", "30-34", "25-29", "20-24", "15-19", "10-14", "5-9", "0-4"))

tau.binned <-ggplot(data = subset(merged, resolution %in% "25"), aes(x = as.numeric(length.gen.bin), y = tau, col = metric:method)) +
  theme_light() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "white"), text = element_text(size=8)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major.y = element_blank()) +
  geom_hline(yintercept = 0, col = "grey", size = 0.5) +
  geom_smooth(data = subset(merged, resolution %in% "25" & method %in% "fold"), method = "lm", size = 0.3, linetype = "solid", se = FALSE) +
  geom_smooth(data = subset(merged, resolution %in% "25" & method %in% "trans"), method = "lm", size = 0.3, linetype = "dashed", se = FALSE) +
  geom_point(data = subset(merged, resolution %in% "25" & method %in% "fold"), position = position_dodge(width = 0.8), alpha = 0.3, shape = 3, size = 0.5) +
  geom_point(data = subset(merged, resolution %in% "25" & method %in% "trans"), position = position_dodge(width = 0.8), alpha = 0.3, shape = 1, size = 0.5) +
  scale_color_manual(values = c("#e00000", "#e00000", "#dbaf00", "#dbaf00")) +
  facet_wrap(~ strength, nrow = 1, scales = "free_x") +
  scale_x_continuous("Length of sampling timeseries (generations)",
                     breaks = c(1,2,3,4,5,6,7,8,9),
                     labels = c("40-44", "35-39", "30-34", "25-29", "20-24", "15-19", "10-14", "5-9", "0-4")) +
  scale_y_continuous(expression(paste("Kendall's ", tau)))

grid.arrange(tau.binned, sim.false.binned)

##### Figure 1: resolution of sampling in simulations #####

# Simulation resolution analysis

merged <- subset(merged, !(merged$strength %in% "none"))

merged$resolution <- revalue(merged$resolution, c("2"="2", "1"="1", "50"="0.5", "25"="0.25"))
merged$resolution <- factor(merged$resolution, levels = c("2", "1", "0.5", "0.25"))
merged$resolution <- factor(merged$resolution, levels(merged$resolution)[c(4, 3, 2, 1)])

sim.res <- ggplot(data = merged, aes(x = resolution, y = tau, fill = method:metric)) +
  geom_hline(yintercept = 0, col = 'gray', size = 0.5) +
  geom_point(aes(col = method:metric), position = position_jitterdodge(dodge.width = 0.7), alpha = 0.4, size = 0.2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.2) +
  theme_light() +
  theme(strip.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=8)) +
  theme(text = element_text(size=8)) +
  scale_fill_manual(values = c("#FFDC54", "#ff9a8e", "#B28E00", "#e00000")) +
  scale_color_manual(values = c("#FFDC54", "#ff9a8e", "#B28E00", "#e00000")) +
  scale_y_continuous(name = expression(paste("Kendall's ", tau)), limits = c(-0.5,0.5)) +
  scale_x_discrete(name = "Resolution of timeseries (/generations)") +
  facet_wrap(~ metric * strength)

##### Process microcosm data #####

load("TraitData.RData")
load("TraitHalfData.RData")

# Make large dataframe

full.res <- rbind(constant, slow, medium, fast)
full.res$resolution <- 50
half.res <- rbind(con.half, slow.half, medium.half, fast.half)
half.res$resolution <- 1
biggy <- rbind(full.res, half.res)
biggy$resolution <- as.factor(biggy$resolution)

biggy.summary <- aggregate(data = biggy, trait ~ length.gen * strength * metric * resolution, mean)
biggy.summary$trait.top <- unlist(aggregate(data = biggy, trait ~ length.gen * strength * metric * resolution, quantile, probs = 0.75)[5])
biggy.summary$trait.bot <- unlist(aggregate(data = biggy, trait ~ length.gen * strength * metric * resolution, quantile, probs = 0.25)[5])

biggy.summary$no.trait <- unlist(aggregate(data = biggy, no.trait ~ length.gen * strength * metric * resolution, mean)[5])
biggy.summary$no.trait.top <- unlist(aggregate(data = biggy, no.trait ~ length.gen * strength * metric * resolution, quantile, probs = 0.75)[5])
biggy.summary$no.trait.bot <- unlist(aggregate(data = biggy, no.trait ~ length.gen * strength * metric * resolution, quantile, probs = 0.25)[5])

biggy.summary$false.neg.no <- unlist(aggregate(data = biggy, no.trait ~ length.gen * strength * metric * resolution, false.neg)[5])
for(i in which(biggy.summary$strength %in% "none")){
  biggy.summary$false.neg.no[i] <- (1 - biggy.summary$false.neg.no[i])
}
biggy.summary$false.neg.trait <- unlist(aggregate(data = biggy, trait ~ length.gen * strength * metric * resolution, false.neg)[5])
for(i in which(biggy.summary$strength %in% "none")){
  biggy.summary$false.neg.trait[i] <- (1 - biggy.summary$false.neg.trait[i])
}

levels(biggy.summary$metric) <- c("ar(1)", "sd")
biggy.summary$metric <- factor(biggy.summary$metric, levels = c("sd","ar(1)"))

# Lowest resolution lengths vector for plotting
lowres.lengths <- unique(subset(biggy, biggy$resolution %in% 1)$length.gen)

# Melted dataframe for boxplots
biggy.melted <- melt(biggy.summary[,-c(6,7,9,10,11,12)], id.vars = c("length.gen", "strength", "resolution", "metric"), value.name = "tau")

# Dataframe of mean tau values for low resolution lengths without forcing
biggy.means <- aggregate(data = biggy, trait ~ length.gen * metric * resolution * strength, mean)
biggy.means$no.trait <- unlist(aggregate(data = biggy, no.trait ~ length.gen * metric * resolution * strength, mean)[5])
biggy.means <- subset(biggy.means, biggy.means$length.gen %in% lowres.lengths & biggy.means$strength %in% "none")

# Dataframe of mean tau values for low resolution lengths with forcing
biggy.means.force <- aggregate(data = subset(biggy, biggy$strength %in% c("slow", "medium", "fast")), trait ~ length.gen * metric * resolution, mean)
biggy.means.force$no.trait <- unlist(aggregate(data = subset(biggy, biggy$strength %in% c("slow", "medium", "fast")), no.trait ~ length.gen * metric * resolution, mean)[4])
biggy.means.force <- subset(biggy.means.force, biggy.means.force$length.gen %in% lowres.lengths)

# Data split by presence of forcing
biggy.melted$forcing <- ifelse(biggy.melted$strength %in% "none", "No forcing", "Forcing")

# Revalue strengths and remove constant forcing
biggy.summary$strength <- revalue(biggy.summary$strength, c("slow"="slow", "medium"="moderate", "fast"="fast"))
biggy.summary <- subset(biggy.summary, strength %in% c("slow", "moderate", "fast"))

biggy.summary.new <- gather(biggy.summary, trait, tau, c(trait, no.trait), factor_key = TRUE)
biggy.summary.new <- gather(biggy.summary.new, trait.neg, false.neg, c(false.neg.no, false.neg.trait), factor_key = TRUE)
biggy.summary.new <- subset(biggy.summary.new, !(biggy.summary.new$strength %in% "none"))


biggy.summary$strength <- revalue(biggy.summary$strength, c("slow"="slow", "medium"="moderate", "fast"="fast"))
biggy.summary <- subset(biggy.summary, strength %in% c("slow", "moderate", "fast"))

biggy.summary2 <- gather(biggy.summary, inclusion, tau, trait, no.trait)
biggy.summary2$inclusion <- as.factor(biggy.summary2$inclusion)

##### Figure 4: Length of sampling in microcosms #####

mic.tau.new <- ggplot(data = subset(biggy.summary.new, biggy.summary.new$resolution %in% 50), aes(x = length.gen, col = trait:metric, y = tau)) +
  theme_light() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "white"), text = element_text(size=8)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major.y = element_blank()) +
  geom_hline(yintercept = 0, col = "grey", size = 0.5) +
  geom_line(aes(linetype = trait), size = 0.3) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("#AB00FF", "#35A69C", "#AB00FF", "#35A69C")) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  scale_x_reverse() + 
  scale_y_continuous(expression(paste("Kendall's ", tau))) +
  facet_wrap(~ strength, scales = "free_x")

mic.false.new <- ggplot(data = subset(biggy.summary.new, biggy.summary.new$resolution %in% 50), aes(x = length.gen, col = trait.neg:metric, y = false.neg)) +
  theme_light() +
  theme(legend.position = "none", text = element_text(size=8)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), panel.grid.major.y = element_blank()) +
  geom_smooth(aes(linetype = trait.neg), method = "loess", se = FALSE, size = 0.3) +
  scale_color_manual(values = c("#AB00FF", "#35A69C", "#AB00FF", "#35A69C")) +
  scale_linetype_manual(values = c("dotted","solid")) +
  scale_x_reverse("Length of sampling timeseries (generations)") +
  scale_y_continuous("Rate of false negatives") +
  facet_wrap(~ strength, scales = "free_x")

grid.arrange(mic.tau.new, mic.false.new)

##### Figure 3: resolution of sampling in microcosms #####

biggy.summary2$resolution <- factor(biggy.summary2$resolution,
                                    levels(biggy.summary2$resolution)[c(2, 1)])

mic.res <- ggplot(data = biggy.summary2, aes(y = tau, x = resolution, fill = inclusion:metric)) +
  theme_light() +
  theme(strip.background = element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=8)) +
  #theme(strip.background = element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=8)) +
  geom_hline(yintercept = 0, col = 'gray', size = 0.5) +
  scale_fill_manual(values = c("#35A69C", "#AB00FF", "#1B544F", "#39004A")) +
  scale_colour_manual(values = c("#35A69C", "#AB00FF", "#1B544F", "#39004A")) +
  scale_y_continuous(limits = c(-1,1), name = expression(paste("Kendall's ", tau))) +
  scale_x_discrete(limits = rev(levels(biggy.summary2$resolution)), name = "Resolution of timeseries (/generations)", labels = c("0.5", "1")) +
  geom_point(aes(col = inclusion:metric), alpha = 0.5, position = position_jitterdodge(dodge.width = 0.7), size = 0.2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.2) +
  facet_wrap(~ metric*strength)

##### Print plots #####

### Figure 1 ###
# Fold on the left, trans on the right
# SD is in blue, AR(1) is in orange

pdf("plots/fig1.pdf", width = 6, height = 4)
sim.res
dev.off()

### Figure 2 ###
# Fold is solid line, trans is dashed line
# SD is in blue, AR(1) is in orange

pdf("plots/fig2.pdf", width = 6, height = 4)
grid.arrange(tau.binned, sim.false.binned)
dev.off()

# Figure 3
# SD is on left, AR(1) is on right
# Trait is darker, no trait is lighter

pdf("plots/fig3.pdf", width = 6, height = 4)
mic.res
dev.off()

# Figure 4
# With trait is in solid, without trait is in dashed
# SD is in yellow, AR(1) is in teal

pdf("plots/fig4.pdf", width = 6, height = 4)
grid.arrange(mic.tau.new, mic.false.new)
dev.off()

