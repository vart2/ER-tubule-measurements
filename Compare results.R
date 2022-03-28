##Load packages. Use install.packages("pkg_name") in the console if first time.
library(tidyverse)
library(ggplot2)
library(qqplotr)

#Modify file paths to results files
file1 <- "/Users/vishal/Desktop/Project/Sample 5/2.5 blur, global, minferet/results>4.csv"
file2 <- "/Users/vishal/Desktop/Project/Sample 6/1.5 blur, global, minferet/results>4.csv"

name1 <- "Wild type"    #Name of first stack
name2 <- "Atlastin mutant"    #Name of second stack

scale1 <- 2     #Scale of images in nm/pixel for first results
scale2 <- 5     #Scale of images in nm/pixel for second results
thickness1 <- 35      #Slice thickness in nm for first results
thickness2 <- 35      #Slice thickness in nm for second results

# Load in results ----

#Read in results
results1 <- read.csv(file1, check.names = FALSE)[-1]
results2 <- read.csv(file2, check.names = FALSE)[-1]

#Convert results to nm
results1_nm <- results1
results1_nm[2:ncol(results1_nm)] <- results1_nm[2:ncol(results1_nm)]*scale1
results2_nm <- results2
results2_nm[2:ncol(results2_nm)] <- results2_nm[2:ncol(results2_nm)]*scale2

#Convert results to long format
results1_long <- results1 %>%      #convert results1 to long format
  pivot_longer(!slice, 
               names_to = "tubule", 
               values_to = "feret_px",
               values_drop_na = TRUE)
results1_long <- results1_long %>%      #Appends column for diameter in nm
  mutate(feret_nm=feret_px*scale1)
results2_long <- results2 %>%      #convert results2 to long format
  pivot_longer(!slice, 
               names_to = "tubule", 
               values_to = "feret_px",
               values_drop_na = TRUE)
results2_long <- results2_long %>%      #Appends column for diameter in nm
  mutate(feret_nm=feret_px*scale2)

#Combine results
feret1 <- data.frame(feret1_nm = results1_long$feret_nm)  #Get feret from results1
feret1$number <- row.names(feret1)
feret2 <- data.frame(feret2_nm = results2_long$feret_nm)  #Get feret from results2
feret2$number <-row.names(feret2)
res_combined <- merge(feret1, feret2, by = "number", all = TRUE)[-1]  #Merge results
names(res_combined)[names(res_combined)=="feret1_nm"] <- name1
names(res_combined)[names(res_combined)=="feret2_nm"] <- name2
res_combined_long <- res_combined %>%      #convert too long format
  pivot_longer(cols = everything(),
               names_to = "sample",
               values_to = "feret_nm",
               values_drop_na = TRUE)

# Create dataframe of tubule lengths ----
lengths1 <- c()
for (i in 2:ncol(results1)) {
  nonNAIndex <- which(!is.na(results1[i]))
  lengths1[i-1] <- max(nonNAIndex)-min(nonNAIndex)+1
}
res_lengths1 <- results1[1, 2:ncol(results1)]
res_lengths1 <- rbind(res_lengths1, lengths1)
res_lengths1 <- res_lengths1[-1, ]
res_lengths1 <- res_lengths1 %>%
  pivot_longer(cols = everything(),
               names_to = "tubule",
               values_to = "length_slices")
res_lengths1 <- res_lengths1 %>%
  mutate(length_nm=length_slices*thickness1)
res_lengths1 <- res_lengths1[-1]
res_lengths1 <- cbind(sample = name1, res_lengths1)
lengths2 <- c()
for (i in 2:ncol(results2)) {
  nonNAIndex <- which(!is.na(results2[i]))
  lengths2[i-1] <- max(nonNAIndex)-min(nonNAIndex)+1
}
res_lengths2 <- results2[1, 2:ncol(results2)]
res_lengths2 <- rbind(res_lengths2, lengths2)
res_lengths2 <- res_lengths2[-1, ]
res_lengths2 <- res_lengths2 %>%
  pivot_longer(cols = everything(),
               names_to = "tubule",
               values_to = "length_slices")
res_lengths2<- res_lengths2 %>%
  mutate(length_nm=length_slices*thickness2)
res_lengths2 <- res_lengths2[-1]
res_lengths2 <- cbind(sample = name2, res_lengths2)
res_lengths_combined <- rbind(res_lengths1, res_lengths2)  #Merge results

# Calculate stats ----
quantiles <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
o <- function(x) {
  subset(x, x < quantile(x, 0.05) | quantile(x, 0.95) < x)
}
mean.error <- function(x) sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))
mean_feret1 <- mean(results1_long$feret_nm)
mean_feret2 <- mean(results2_long$feret_nm)
means <- data.frame("sample"=c(name1, name2), "mean"=c(mean_feret1, mean_feret2))
sd_feret1 <- sd(results1_long$feret_nm)
sd_feret2 <- sd(results2_long$feret_nm)
sds <- data.frame("sample"=c(name1, name2), "sd"=c(sd_feret1, sd_feret2))
sem_feret1 <- mean.error(results1_long$feret_nm)
sem_feret2 <- mean.error(results2_long$feret_nm)
sems <- data.frame("sample"=c(name1, name2), "sem"=c(sem_feret1, sem_feret2))
percentiles_feret1 <- quantile(results1_long$feret_nm, c(0.05, 0.1, 0.5, 0.9, 0.95))
percentiles_feret2 <- quantile(results2_long$feret_nm, c(0.05, 0.1, 0.5, 0.9, 0.95))
percentiles <- data.frame("sample"=c(name1, name2),
                          "5th percentile"=c(percentiles_feret1[1], percentiles_feret2[1]),
                          "10th percentile"=c(percentiles_feret1[2], percentiles_feret2[2]),
                          "50th percentile"=c(percentiles_feret1[3], percentiles_feret2[3]),
                          "90th percentile"=c(percentiles_feret1[4], percentiles_feret2[4]),
                          "95th percentile"=c(percentiles_feret1[5], percentiles_feret2[5]),
                          check.names = FALSE)
stats <- merge(means, sds, by="sample")
stats <- merge(stats, sems, by="sample")
stats <- merge(stats, percentiles, by="sample")

# Calculate individual tubule statistics ----
tubule_stats1 <- results1_nm[1, 2:ncol(results1_nm)]
tubule_mean1 <- c()
for (i in 2:ncol(results1_nm)) {
  tubule_mean1[i-1] <- mean(results1_nm[-1, i], na.rm = TRUE)
}
tubule_sd1 <- c()
for (i in 2:ncol(results1_nm)) {
  tubule_sd1[i-1] <- sd(results1_nm[-1, i], na.rm = TRUE)
}
tubule_sem1 <- c()
for (i in 2:ncol(results1_nm)) {
  tubule_sem1[i-1] <- mean.error(results1_nm[-1, i])
}
"tubule_5%1" <- c()
for (i in 2:ncol(results1_nm)) {
  `tubule_5%1`[i-1] <- quantile(results1_nm[-1, i], 0.05, na.rm = TRUE)
}
"tubule_10%1" <- c()
for (i in 2:ncol(results1_nm)) {
  `tubule_10%1`[i-1] <- quantile(results1_nm[-1, i], 0.1, na.rm = TRUE)
}
"tubule_25%1" <- c()
for (i in 2:ncol(results1_nm)) {
  `tubule_25%1`[i-1] <- quantile(results1_nm[-1, i], 0.25, na.rm = TRUE)
}
"tubule_50%1" <- c()
for (i in 2:ncol(results1_nm)) {
  `tubule_50%1`[i-1] <- quantile(results1_nm[-1, i], 0.5, na.rm = TRUE)
}
"tubule_75%1" <- c()
for (i in 2:ncol(results1_nm)) {
  `tubule_75%1`[i-1] <- quantile(results1_nm[-1, i], 0.75, na.rm = TRUE)
}
"tubule_90%1" <- c()
for (i in 2:ncol(results1_nm)) {
  `tubule_90%1`[i-1] <- quantile(results1_nm[-1, i], 0.9, na.rm = TRUE)
}
"tubule_95%1" <- c()
for (i in 2:ncol(results1_nm)) {
  `tubule_95%1`[i-1] <- quantile(results1_nm[-1, i], 0.95, na.rm = TRUE)
}
tubule_stats1 <- rbind(tubule_stats1, tubule_mean1, tubule_sd1, tubule_sem1, `tubule_5%1`, `tubule_10%1`, `tubule_25%1`, `tubule_50%1`, `tubule_75%1`, `tubule_90%1`, `tubule_95%1`)
tubule_stats1 <- tubule_stats1[-1, ]
rownames(tubule_stats1) <- c("mean", "sd", "sem", "5%", "10%", "25%", "50%", "75%", "90%", "95%")
tubule_stats1 <- as.data.frame(t(tubule_stats1))

tubule_stats2 <- results2_nm[1, 2:ncol(results2_nm)]
tubule_mean2 <- c()
for (i in 2:ncol(results2_nm)) {
  tubule_mean2[i-1] <- mean(results2_nm[-1, i], na.rm = TRUE)
}
tubule_sd2 <- c()
for (i in 2:ncol(results2_nm)) {
  tubule_sd2[i-1] <- sd(results2_nm[-1, i], na.rm = TRUE)
}
tubule_sem2 <- c()
for (i in 2:ncol(results2_nm)) {
  tubule_sem2[i-1] <- mean.error(results2_nm[-1, i])
}
"tubule_5%2" <- c()
for (i in 2:ncol(results2_nm)) {
  `tubule_5%2`[i-1] <- quantile(results2_nm[-1, i], 0.05, na.rm = TRUE)
}
"tubule_10%2" <- c()
for (i in 2:ncol(results2_nm)) {
  `tubule_10%2`[i-1] <- quantile(results2_nm[-1, i], 0.1, na.rm = TRUE)
}
"tubule_25%2" <- c()
for (i in 2:ncol(results2_nm)) {
  `tubule_25%2`[i-1] <- quantile(results2_nm[-1, i], 0.25, na.rm = TRUE)
}
"tubule_50%2" <- c()
for (i in 2:ncol(results2_nm)) {
  `tubule_50%2`[i-1] <- quantile(results2_nm[-1, i], 0.5, na.rm = TRUE)
}
"tubule_75%2" <- c()
for (i in 2:ncol(results2_nm)) {
  `tubule_75%2`[i-1] <- quantile(results2_nm[-1, i], 0.75, na.rm = TRUE)
}
"tubule_90%2" <- c()
for (i in 2:ncol(results2_nm)) {
  `tubule_90%2`[i-1] <- quantile(results2_nm[-1, i], 0.9, na.rm = TRUE)
}
"tubule_95%2" <- c()
for (i in 2:ncol(results2_nm)) {
  `tubule_95%2`[i-1] <- quantile(results2_nm[-1, i], 0.95, na.rm = TRUE)
}
tubule_stats2 <- rbind(tubule_stats2, tubule_mean2, tubule_sd2, tubule_sem2, `tubule_5%2`, `tubule_10%2`, `tubule_25%2`, `tubule_50%2`, `tubule_75%2`, `tubule_90%2`, `tubule_95%2`)
tubule_stats2 <- tubule_stats2[-1, ]
rownames(tubule_stats2) <- c("mean", "sd", "sem", "5%", "10%", "25%", "50%", "75%", "90%", "95%")
tubule_stats2 <- as.data.frame(t(tubule_stats2))
tubule_stats1_comb <- cbind(sample = name1, tubule_stats1)
tubule_stats2_comb <- cbind(sample = name2, tubule_stats2)
tubule_stats_comb <- rbind(tubule_stats1_comb, tubule_stats2_comb)


# Graphs ----
#Frequency graph comparing diameters
ggplot(NULL, aes(x=feret_nm, y = (..count..)*100/sum(..count..))) +
  geom_freqpoly(data=results1_long, binwidth=scale1, aes(colour="green4")) +
  geom_freqpoly(data=results2_long, binwidth=scale2, aes(colour="dodgerblue4")) +
  geom_vline(xintercept = mean_feret1,
             linetype = "dashed",
             colour = "green4") +
  geom_vline(xintercept = mean_feret2,
             linetype = "dashed",
             colour = "dodgerblue4") +
  scale_colour_identity(name = NULL,
                        labels = c("green4"=name1, "dodgerblue4"=name2),
                        guide = "legend") +
  xlab("Diameter / nm") +
  ylab("Relative frequency / %") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

#Boxplot comparing diameters
ggplot(data=res_combined_long, aes(x=sample, y=feret_nm)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  #stat_summary(fun=mean, geom="point",colour="darkred", size=3,
  #             aes(y=feret_nm)) +
  #stat_summary(fun=mean, geom="point",colour="darkred", size=3,
  #             aes(y=feret_nm)) +
  #stat_summary(fun=mean, geom="text", vjust=-1,
  #             aes(x=name1, y=mean(results1_long$feret_nm)),
  #            label = paste("mean =", round(mean_feret1),"nm")) +
  #stat_summary(fun=mean, geom="text", vjust=-1,
  #             aes(x=name2, y=mean(results2_long$feret_nm)),
  #             label = paste("mean =", round(mean_feret2),"nm")) +
  xlab(NULL) +
  ylab("Diameter / nm") +
  scale_y_continuous(breaks = scales::breaks_pretty(10))

#Frequency graph comparing tubule lengths in nm
ggplot(NULL, aes(y = (..count..)*100/sum(..count..))) +
  geom_freqpoly(data=res_lengths1, binwidth=thickness1, aes(x=length_nm, colour="green4")) +
  geom_freqpoly(data=res_lengths2, binwidth=thickness2, aes(x=length_nm, colour="dodgerblue4")) +
  geom_vline(xintercept = mean(res_lengths1$length_nm),
             linetype = "dashed",
             colour = "green4") +
  geom_vline(xintercept = mean(res_lengths2$length_nm),
             linetype = "dashed",
             colour = "dodgerblue4") +
  scale_colour_identity(name = NULL,
                        labels = c("green4"=name1, "dodgerblue4"=name2),
                        guide = "legend") +
  xlab("Tubule length / nm") +
  ylab("Relative frequency / %") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

#Boxplot comparing tubule lengths in nm
ggplot(data=res_lengths_combined, aes(x=sample, y=length_nm)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  xlab(NULL) +
  ylab("Tubule length / nm") +
  scale_y_continuous(breaks = scales::breaks_pretty(10))

#Frequency graph comparing tubule lengths in slices
ggplot(NULL, aes(y = (..count..)*100/sum(..count..))) +
  geom_freqpoly(data=res_lengths1, binwidth=1, aes(x=length_slices, colour="green4")) +
  geom_freqpoly(data=res_lengths2, binwidth=1, aes(x=length_slices, colour="dodgerblue4")) +
  geom_vline(xintercept = mean(res_lengths1$length_slices),
             linetype = "dashed",
             colour = "green4") +
  geom_vline(xintercept = mean(res_lengths2$length_slices),
             linetype = "dashed",
             colour = "dodgerblue4") +
  scale_colour_identity(name = NULL,
                        labels = c("green4"=name1, "dodgerblue4"=name2),
                        guide = "legend") +
  xlab("Tubule length / slices") +
  ylab("Relative frequency / %") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

#Boxplot comparing tubule lengths in slices
ggplot(data=res_lengths_combined, aes(x=sample, y=length_slices)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  xlab(NULL) +
  ylab("Tubule length / slcies") +
  scale_y_continuous(breaks = scales::breaks_pretty(10))

#Frequency plot comparing mean tubule diameters
ggplot(NULL, aes(x=mean, y = (..count..)*100/sum(..count..))) +
  geom_freqpoly(data=tubule_stats1, binwidth = scale1, aes(colour="green4")) +
  geom_freqpoly(data=tubule_stats2, binwidth = scale2, aes(colour="dodgerblue4")) +
  geom_vline(xintercept = mean(tubule_stats1$mean),
             linetype = "dashed",
             colour = "green4") +
  geom_vline(xintercept = mean(tubule_stats2$mean),
             linetype = "dashed",
             colour = "dodgerblue4") +
  scale_colour_identity(name = NULL,
                        labels = c("green4"=name1, "dodgerblue4"=name2),
                        guide = "legend") +
  xlab("Mean diameter of tubule / nm") +
  ylab("Relative frequency / %") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

#Boxplot comparing mean tubule diameter
ggplot(tubule_stats_comb, aes(x=sample, y=mean)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  xlab(NULL) +
  scale_y_continuous(breaks = scales::breaks_pretty(10)) +
  ylab("Mean diameter of tubule / nm")

#Scatter plot comparing mean tubule diameters
ggplot(tubule_stats_comb) +
  geom_point(aes(x=sample, y=mean)) +
  geom_point(aes(x=name1, y=mean(tubule_stats1$mean)), colour="darkred", size=2) +
  geom_point(aes(x=name2, y=mean(tubule_stats2$mean)), colour="darkred", size=2) +
  geom_errorbar(data=stats, colour="darkred",
                aes(x=name1,
                    ymin=mean(tubule_stats1$mean)-mean.error(tubule_stats1$mean),
                    ymax=mean(tubule_stats1$mean)+mean.error(tubule_stats1$mean))) +
  geom_errorbar(data=stats, colour="darkred",
                aes(x=name2,
                    ymin=mean(tubule_stats2$mean)-mean.error(tubule_stats2$mean),
                    ymax=mean(tubule_stats2$mean)+mean.error(tubule_stats2$mean))) +
  xlab(NULL) +
  ylab("Mean diameter of tubule / nm") +
  scale_y_continuous(breaks = scales::breaks_pretty(10))

#Frequency plot comparing median tubule diameters
ggplot(NULL, aes(x=`50%`, y = (..count..)*100/sum(..count..))) +
  geom_freqpoly(data=tubule_stats1, binwidth = scale1, aes(colour = "green4")) +
  geom_freqpoly(data=tubule_stats2, binwidth = scale2, aes(colour = "dodgerblue4")) +
  geom_vline(xintercept = mean(tubule_stats1$`50%`),
             linetype = "dashed",
             colour = "green4") +
  geom_vline(xintercept = mean(tubule_stats2$`50%`),
             linetype = "dashed",
             colour = "dodgerblue4") +
  scale_colour_identity(name = NULL,
                        labels = c("green4"=name1, "dodgerblue4"=name2),
                        guide = "legend") +
  xlab("Median diameter of tubule / nm") +
  ylab("Relative frequency / %") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

#Boxplot comparing median tubule diameter
ggplot(tubule_stats_comb, aes(x=sample, y=`50%`)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  xlab(NULL) +
  scale_y_continuous(breaks = scales::breaks_pretty(10)) +
  ylab("Median diameter of tubule / nm")

#Frequency plot comparing IQR of tubule diameters
ggplot(NULL, aes(x=`75%`-`25%`, y = (..count..)*100/sum(..count..))) +
  geom_freqpoly(data=tubule_stats1, binwidth = scale1, aes(colour = "green4")) +
  geom_freqpoly(data=tubule_stats2, binwidth = scale2, aes(colour = "dodgerblue4")) +
  geom_vline(xintercept = mean(tubule_stats1$`75%`-tubule_stats1$`25%`),
             linetype = "dashed",
             colour = "green4") +
  geom_vline(xintercept = mean(tubule_stats2$`75%`-tubule_stats2$`25%`),
             linetype = "dashed",
             colour = "dodgerblue4") +
  scale_colour_identity(name = NULL,
                        labels = c("green4"=name1, "dodgerblue4"=name2),
                        guide = "legend") +
  xlab("IQR of diameter of tubule / nm") +
  ylab("Relative frequency / %") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

#Boxplot comparing IQR of tubule diameters
ggplot(tubule_stats_comb, aes(x=sample, y=`75%`-`25%`)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  scale_y_continuous(breaks = scales::breaks_pretty(10)) +
  xlab(NULL) +
  ylab("IQR of diameter of tubule / nm")

#Frequency plot comparing standard deviation of tubule diameters
ggplot(NULL, aes(x=sd, y = (..count..)*100/sum(..count..))) +
  geom_freqpoly(data=tubule_stats1, binwidth = scale1, aes(colour = "green4")) +
  geom_freqpoly(data=tubule_stats2, binwidth = scale2, aes(colour = "dodgerblue4")) +
  geom_vline(xintercept = mean(tubule_stats1$sd),
             linetype = "dashed",
             colour = "green4") +
  geom_vline(xintercept = mean(tubule_stats2$sd),
             linetype = "dashed",
             colour = "dodgerblue4") +
  scale_colour_identity(name = NULL,
                        labels = c("green4"=name1, "dodgerblue4"=name2),
                        guide = "legend") +
  xlab("SD of diameter of tubule / nm") +
  ylab("Relative frequency / %") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

#Boxplot comparing standard deviation of tubule diameters
ggplot(tubule_stats_comb, aes(x=sample, y=sd)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  scale_y_continuous(breaks = scales::breaks_pretty(10)) +
  xlab(NULL) +
  ylab("SD of diameter of tubule / nm")


# Quantile-Quantile plots ----
#Q-Q plot of all diameters to test for normal distribution
ggplot(res_combined_long, aes(sample=feret_nm)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  facet_wrap(~sample) +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of tubule lengths in nm to test for normal distribution
ggplot(res_lengths_combined, aes(sample=length_nm)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  facet_wrap(~sample) +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of tubule lengths in slices to test for normal distribution
ggplot(res_lengths_combined, aes(sample=length_slices)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  facet_wrap(~sample) +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of mean tubule diameters to test for normal distribution
ggplot(tubule_stats_comb, aes(sample=mean)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  facet_wrap(~sample) +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of IQR of tubule diameters to test for normal distribution
ggplot(tubule_stats_comb, aes(sample=`75%`-`25%`)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  facet_wrap(~sample) +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of SD of tubule diameters to test for normal distribution
ggplot(tubule_stats_comb, aes(sample=sd)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  facet_wrap(~sample) +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")



# Statistical tests ----
#Mann-Whitney tests to check for statistical significance
wilcox.test(results1_long$feret_nm, results2_long$feret_nm)
wilcox.test(res_lengths1$length_nm, res_lengths2$length_nm)
wilcox.test(res_lengths1$length_slices, res_lengths2$length_slices)
wilcox.test(tubule_stats1$mean, tubule_stats2$mean)
wilcox.test(tubule_stats1$`75%`-tubule_stats1$`25%`, tubule_stats2$`75%`-tubule_stats2$`25%`)
wilcox.test(tubule_stats1$sd, tubule_stats2$sd)



#Q-Q plot of all diameters to test for normal distribution
ggplot(results2_long, aes(sample=feret_nm)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of tubule lengths in nm to test for normal distribution
ggplot(res_lengths2, aes(sample=length_nm)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of tubule lengths in slices to test for normal distribution
ggplot(res_lengths2, aes(sample=length_slices)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of mean tubule diameters to test for normal distribution
ggplot(tubule_stats2, aes(sample=mean)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of IQR of tubule diameters to test for normal distribution
ggplot(tubule_stats2, aes(sample=`75%`-`25%`)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")

#Q-Q plot of SD of tubule diameters to test for normal distribution
ggplot(tubule_stats2, aes(sample=sd)) +
  stat_qq_band(conf=0.95) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles")
