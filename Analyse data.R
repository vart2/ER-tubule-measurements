##Load packages. Use install.packages("pkg_name") in the console if first time. 
library(tidyverse)
library(ggplot2)
library(imputeTS)

#Set file path to results file
file <- readline("Enter path to results file:")

#Set scale of images and slice thickness
scale <- as.numeric(readline("Enter scale of images in nm/pixel:"))
thickness <- as.numeric(readline("Enter slice thickness in nm:"))

# Data  ----

## Load results ----

#Read in results
results <- read.csv(file, check.names = FALSE)[-1]

#Convert results to nm
results_nm <- results
results_nm[2:ncol(results_nm)] <- results_nm[2:ncol(results_nm)]*scale

## Create dataframe of tubule lengths ----
lengths <- c()
for (i in 2:ncol(results)) {
  nonNAIndex <- which(!is.na(results[i]))
  lengths[i-1] <- max(nonNAIndex)-min(nonNAIndex)+1
}
res_lengths <- results[1, 2:ncol(results)]
res_lengths <- rbind(res_lengths, lengths)
res_lengths <- res_lengths[-1, ]
res_lengths <- res_lengths %>%
  pivot_longer(cols = everything(),
               names_to = "tubule",
               values_to = "length_px")
res_lengths<- res_lengths %>%
  mutate(length_nm=length_px*thickness)

## Get longest 10 tubules ----

#Order results by track length
res_ordered <- results[-1]
res_ordered <- res_ordered[, order(tibble::deframe(res_lengths[1:2]), decreasing = TRUE)]
res_ordered <- cbind(slice = c(seq(1:nrow(results))), res_ordered)

results_long <- res_ordered %>%      #convert to long format
  pivot_longer(!slice, 
               names_to = "tubule", 
               values_to = "feret_px",
               values_drop_na = TRUE)
results_long<- results_long %>%
  mutate(feret_nm=feret_px*scale)

#Selects top 10 tubules by track length
longest10 <- res_ordered[1:11]

longest10_long <- longest10 %>%      #convert to long format
  pivot_longer(!slice, 
               names_to = "tubule", 
               values_to = "feret_px",
               values_drop_na = TRUE)
longest10_long <- longest10_long %>%      #Appends column for diameter in nm
  mutate(feret_nm=feret_px*scale)

longest10[2:ncol(longest10)] <- longest10[2:ncol(longest10)]*scale

#Removes NAs of top 10 tubules
longest10_na.rm <- data.frame(number = 1:nrow(results))
temp_res <- rbind(c(1), longest10)
for (i in 2:ncol(longest10)) {
  nonNAIndex <- which(!is.na(longest10[i]))
  remove <- nonNAIndex[1]
  tubule <- data.frame(temp_res[-1:-remove, i])
  colnames(tubule) <- colnames(longest10[i])
  tubule$number <- seq(nrow(tubule))
  longest10_na.rm <- merge(longest10_na.rm, tubule, by="number", all=TRUE)
}
longest10_na.rm <- longest10_na.rm[rowSums(is.na(longest10_na.rm)) != ncol(longest10_na.rm)-1, ]
for (i in 2:ncol(longest10_na.rm)) {
  na_check <- which(!is.na(longest10_na.rm[i]))
  if (na_check[length(na_check)] < nrow(longest10_na.rm)) {
    remove_start <- na_check[length(na_check)] + 1
    remove_end <- nrow(longest10_na.rm)
    longest10_na.rm[remove_start:remove_end, i] <- 0
  }
}
longest10_na.rm <- na_ma(longest10_na.rm, k=1)
longest10_na.rm[longest10_na.rm==0] <- NA
colnames(longest10_na.rm) <- c("number",9,8,7,6,5,4,3,2,1,0)
longest10_na.rm <- longest10_na.rm %>%
  pivot_longer(!number,
               names_to = "tubule",
               values_to = "diameter_nm",
               values_drop_na = TRUE)

## Calculate stats ----
quantiles <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
o <- function(x) {
  subset(x, x < quantile(x, 0.05) | quantile(x, 0.95) < x)
}
mean.error <- function(x) sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))
mean_feret <- mean(results_long$feret_nm)
sd_feret <- sd(results_long$feret_nm)
sem_feret <- mean.error(results_long$feret_nm)
percentiles_feret <- quantile(results_long$feret_nm, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
stats <- data.frame(cbind(mean_feret,
                          sd_feret,
                          sem_feret,
                          percentiles_feret[1],
                          percentiles_feret[2],
                          percentiles_feret[3],
                          percentiles_feret[4],
                          percentiles_feret[5],
                          percentiles_feret[6],
                          percentiles_feret[7]))
rownames(stats) <- "1"
colnames(stats) <- c("mean", "sd", "sem", "5%", "10%", "25%", "50%", "75%", "90%", "95%")

## Calculate individual tubule statistics ----
tubule_stats <- results_nm[1, 2:ncol(results_nm)]
tubule_mean <- c()
for (i in 2:ncol(results_nm)) {
  tubule_mean[i-1] <- mean(results_nm[-1, i], na.rm = TRUE)
}
tubule_sd <- c()
for (i in 2:ncol(results_nm)) {
  tubule_sd[i-1] <- sd(results_nm[-1, i], na.rm = TRUE)
}
tubule_sem <- c()
for (i in 2:ncol(results_nm)) {
  tubule_sem[i-1] <- mean.error(results_nm[-1, i])
}
"tubule_5%" <- c()
for (i in 2:ncol(results_nm)) {
  `tubule_5%`[i-1] <- quantile(results_nm[-1, i], 0.05, na.rm = TRUE)
}
"tubule_10%" <- c()
for (i in 2:ncol(results_nm)) {
  `tubule_10%`[i-1] <- quantile(results_nm[-1, i], 0.1, na.rm = TRUE)
}
"tubule_25%" <- c()
for (i in 2:ncol(results_nm)) {
  `tubule_25%`[i-1] <- quantile(results_nm[-1, i], 0.25, na.rm = TRUE)
}
"tubule_50%" <- c()
for (i in 2:ncol(results_nm)) {
  `tubule_50%`[i-1] <- quantile(results_nm[-1, i], 0.5, na.rm = TRUE)
}
"tubule_75%" <- c()
for (i in 2:ncol(results_nm)) {
  `tubule_75%`[i-1] <- quantile(results_nm[-1, i], 0.75, na.rm = TRUE)
}
"tubule_90%" <- c()
for (i in 2:ncol(results_nm)) {
  `tubule_90%`[i-1] <- quantile(results_nm[-1, i], 0.9, na.rm = TRUE)
}
"tubule_95%" <- c()
for (i in 2:ncol(results_nm)) {
  `tubule_95%`[i-1] <- quantile(results_nm[-1, i], 0.95, na.rm = TRUE)
}
tubule_stats <- rbind(tubule_stats, tubule_mean, tubule_sd, tubule_sem, `tubule_5%`, `tubule_10%`, `tubule_25%`, `tubule_50%`, `tubule_75%`, `tubule_90%`, `tubule_95%`)
tubule_stats <- tubule_stats[-1, ]
rownames(tubule_stats) <- c("mean", "sd", "sem", "5%", "10%", "25%", "50%", "75%", "90%", "95%")
tubule_stats <- as.data.frame(t(tubule_stats))



# Graphs ----

## Plot top 10 tubules ----
ggplot(longest10_long, aes(slice, feret_nm, col=tubule)) +
  geom_line() +
  facet_wrap(~tubule, scales = "free_x") +
  xlab("Slice") +
  ylab("Diameter / nm") +
  theme(legend.position = "none")

## Stacked area chart of 10 longest tubules ----
ggplot(longest10_na.rm, aes(x=number*thickness, y=diameter_nm, fill=tubule)) +
  geom_area(position="stack", stat="identity") +
  xlab("Length along tubule / nm") +
  ylab("Diameter / nm") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(thickness,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"), legend.position = "none")

## Boxplot comparing diameter of 10 longest tubules ----
ggplot(longest10_long, aes(x=tubule, y=feret_nm)) + 
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  xlab("Tubule") +
  ylab("Diameter / nm") +
  scale_y_continuous(breaks = scales::breaks_pretty(10))

## Frequency graph of diameter for whole dataset ----
ggplot(results_long, aes(x=feret_nm)) +
  geom_freqpoly(binwidth=scale, aes(y = (..count..)*100/sum(..count..))) +
  xlab("Diameter / nm") +
  ylab("Relative frequency / %") +
  geom_vline(xintercept = mean(results_long$feret_nm), linetype = "dashed") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

## Boxplot of diameters ----
ggplot(data=results_long, aes(y=feret_nm, x=1)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  #geom_boxplot(aes(y=feret_nm)) +
  #stat_summary(fun=mean, geom="point",colour="darkred", size=3,
  #             aes(x=1, y=feret_nm)) +
  #stat_summary(fun=mean, geom="text", vjust=-1,
  #             aes(x=1, y=feret_nm),
  #             label = c(paste("mean =", round(mean_feret),"nm"))) +
  #geom_errorbar(data=stats, colour="darkred",
  #              aes(x=1, ymin=mean_feret-sem_feret, ymax=mean_feret+sem_feret,width=0.4)) +
  xlab(NULL) +
  ylab("Diameter / nm") +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(breaks = scales::breaks_pretty(10))

## Frequency graph of tubule lengths ----
ggplot(res_lengths, aes(x=length_nm)) +
  geom_freqpoly(binwidth = thickness, aes(y = (..count..)*100/sum(..count..))) +
  xlab("Tubule length / nm") +
  ylab("Relative frequency / %") +
  geom_vline(xintercept = mean(res_lengths$length_nm), linetype = "dashed") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

## Boxplot of tubule lengths ----
ggplot(res_lengths, aes(y=length_nm, x=1)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(breaks = scales::breaks_pretty(10)) +
  xlab(NULL) +
  ylab("Tubule length / nm")

## Frequency plot of mean tubule diameter ----
ggplot(tubule_stats, aes(x=mean)) +
  geom_freqpoly(binwidth = scale, aes(y = (..count..)*100/sum(..count..))) +
  xlab("Mean diameter of tubule / nm") +
  ylab("Relative frequency / %") +
  geom_vline(xintercept = mean(tubule_stats$mean), linetype = "dashed") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

## Boxplot of mean tubule diameter ----
ggplot(tubule_stats, aes(y=mean, x=1)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(breaks = scales::breaks_pretty(10)) +
  xlab(NULL) +
  ylab("Mean diameter of tubule / nm")

## Frequency plot of median tubule diameter ----
ggplot(tubule_stats, aes(x=`50%`)) +
  geom_freqpoly(binwidth = scale, aes(y = (..count..)*100/sum(..count..))) +
  xlab("Median diameter of tubule / nm") +
  ylab("Relative frequency / %") +
  geom_vline(xintercept = mean(tubule_stats$`50%`), linetype = "dashed") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

## Boxplot of median tubule diameter ----
ggplot(tubule_stats, aes(y=`50%`, x=1)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(breaks = scales::breaks_pretty(10)) +
  xlab(NULL) +
  ylab("Median diameter of tubule / nm")

## Frequency plot of IQR of tubule diameter ----
ggplot(tubule_stats, aes(x=`75%`-`25%`)) +
  geom_freqpoly(binwidth = scale, aes(y = (..count..)*100/sum(..count..))) +
  xlab("IQR of diameter of tubule / nm") +
  ylab("Relative frequency / %") +
  geom_vline(xintercept = mean(tubule_stats$`75%`-tubule_stats$`25%`), linetype = "dashed") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

## Boxplot of IQR of tubule diameter ----
ggplot(tubule_stats, aes(y=`75%`-`25%`, x=1)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(breaks = scales::breaks_pretty(10)) +
  xlab(NULL) +
  ylab("IQR of diameter of tubule / nm")

## Frequency plot of standard deviation of tubule diameter ----
ggplot(tubule_stats, aes(x=sd)) +
  geom_freqpoly(binwidth = scale, aes(y = (..count..)*100/sum(..count..))) +
  xlab("SD of diameter of tubule / nm") +
  ylab("Relative frequency / %") +
  geom_vline(xintercept = mean(tubule_stats$sd), linetype = "dashed") +
  scale_x_continuous(breaks = scales::breaks_pretty(10), expand = expansion(mult=c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
  coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) +
  theme(axis.line = element_line(linetype = "solid"))

## Boxplot of SD of tubule diameter ----
ggplot(tubule_stats, aes(y=sd, x=1)) +
  stat_summary(geom="errorbar", width = 0.5, fun.data=quantiles) +
  stat_summary(geom="boxplot", fun.data=quantiles) +
  #stat_summary(fun = o, geom="point") +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(breaks = scales::breaks_pretty(10)) +
  xlab(NULL) +
  ylab("SD of diameter of tubule / nm")

ggplot(results, aes(x=slice*thickness - 280, y=`8, 89`)) +
  geom_line() +
  ylab("Diameter / nm") +
  xlab("Tubule length / nm") +
  coord_cartesian(xlim=c(0, 350))