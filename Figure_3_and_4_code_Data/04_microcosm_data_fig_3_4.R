## code written by Alex Arkilanian

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
