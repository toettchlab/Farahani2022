# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
library(xlsx)
rm(list = ls())

##########################################################################

# enter parameters

t_norm <- 10; # time point for normalizing data
cond_norm <- "cd3e"; # condition for normalizing data

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/Western blots/2022_06_21_EGF_timecourse");

names <-  list.files(pattern="*.csv");
list <- lapply(names, read.csv);
names <- t(str_replace(names,".csv",""));

##########################################################################

# data analysis

data <- list[[1]]

# normalize by actin bands

actin <- filter(data, protein == "actin (pEGFR)" | protein == "actin (EGFR)")

norm <- data.frame(matrix(ncol = 1, nrow = nrow(data)));
colnames(norm) <- c("norm_actin");

for (j in 1:nrow(data)){
  if (data$condition[j] == "EGFR" | data$condition[j] == "actin (EGFR)"){
    norm[j,1] <- data$value[j] / filter(actin, time == data$time[j] & 
                                                date == data$date[j] & 
                                                condition == data$condition[j] & 
                                                protein == "actin (EGFR)")["value"]
  } else {
    norm[j,1] <- data$value[j] / filter(actin, time == data$time[j] & 
                                                  date == data$date[j] & 
                                                  condition == data$condition[j] & 
                                                  protein == "actin (pEGFR)")["value"]
  }
}

data <- cbind(data,norm)

# calculate pEGFR/EGFR

norm <- data.frame(matrix(ncol = 1, nrow = nrow(data)));
colnames(norm) <- c("norm_EGFR");

for (j in 1:nrow(data)){
  if (data$protein[j] == "pEGFR"){
    norm[j,1] <- data$norm_actin[j] / filter(data, time == data$time[j] & 
                                                date == data$date[j] & 
                                                condition == data$condition[j] & 
                                                protein == "EGFR")["norm_actin"]
  } else {
    norm[j,1] <- 0
  }
}

data <- cbind(data,norm)

# normalize to single time point for each replicate

norm <- data.frame(matrix(ncol = 2, nrow = nrow(data)));
colnames(norm) <- c("norm_actin_time","norm_EGFR_time");

for (j in 1:nrow(data)){
  norm[j,1] <- data[j,"norm_actin"] / filter(data, date == data$date[j] &
                                                          protein == data$protein[j] &
                                                          condition == cond_norm &
                                                          time == t_norm)["norm_actin"]
  norm[j,2] <- data[j,"norm_EGFR"] / filter(data, date == data$date[j] &
                                                          protein == data$protein[j] &
                                                          condition == cond_norm &
                                                          time == t_norm)["norm_EGFR"]
}

data <- cbind(data,norm)

output_stats <- data %>%
  group_by(time,condition,protein) %>% 
  summarise(mean_norm_actin = mean(norm_actin),
            mean_norm_EGFR = mean(norm_EGFR),
            sd_norm_actin = sd(norm_actin),
            sd_norm_EGFR = sd(norm_EGFR),
            sem_norm_actin = sd_norm_actin / sqrt(3),
            sem_norm_EGFR = sd_norm_EGFR / sqrt(3),
            mean_norm_actin_time = mean(norm_actin_time),
            mean_norm_EGFR_time = mean(norm_EGFR_time),
            sd_norm_actin_time = sd(norm_actin_time),
            sd_norm_EGFR_time = sd(norm_EGFR_time),
            sem_norm_actin_time = sd_norm_actin_time / sqrt(3),
            sem_norm_EGFR_time = sd_norm_EGFR_time / sqrt(3))

##########################################################################

# plot themes

theme_black <- theme_classic() +
  theme(legend.position="right",
        legend.background = element_rect(fill = 'black', colour = 'black'),
        legend.text = element_text(colour="white"),
        plot.background = element_rect(fill = 'black', colour = 'black'),
        panel.background = element_rect(fill = 'black', colour = 'black'),
        axis.text = element_text(colour="white"),
        axis.title = element_text(colour="white"),
        axis.line = element_line(colour="white")) + 
  theme(aspect.ratio=1)

theme_white <-   theme_classic() +
  theme(legend.position="right",
        legend.background = element_rect(fill = 'white', colour = 'white'),
        legend.text = element_text(colour="black"),
        plot.background = element_rect(fill = 'white', colour = 'white'), 
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.text = element_text(colour="black"),
        axis.title = element_text(colour="black"),
        axis.line = element_line(colour="black")) + 
  theme(aspect.ratio=1)

##########################################################################

p_pEGFR_actin <- ggplot() +
  geom_line(data = filter(output_stats,protein == "pEGFR"), aes(x = time,
                                                                y = mean_norm_actin_time,
                                                                group = condition,
                                                                color = condition)) +
  geom_errorbar(data = filter(output_stats,protein == "pEGFR"), aes(x = time,
                                                                    ymax = mean_norm_actin_time + sem_norm_actin_time,
                                                                    ymin = mean_norm_actin_time - sem_norm_actin_time,
                                                                    group = condition,
                                                                    color = condition)) + 
  ylab("pEGFR/actin")

p_pEGFR_EGFR <- ggplot() +
  geom_point(data = filter(output_stats,protein == "pEGFR"), aes(x = time,
                                                                y = mean_norm_EGFR_time,
                                                                group = condition,
                                                                color = condition,
                                                                shape = condition)) +
  geom_line(data = filter(output_stats,protein == "pEGFR"), aes(x = time,
                                                        y = mean_norm_EGFR_time,
                                                        group = condition,
                                                        color = condition)) +
  geom_errorbar(data = filter(output_stats,protein == "pEGFR"), aes(x = time,
                                                            ymax = mean_norm_EGFR_time + sem_norm_EGFR_time,
                                                            ymin = mean_norm_EGFR_time - sem_norm_EGFR_time,
                                                            group = condition,
                                                            color = condition)) + 
  ylab("pEGFR/EGFR")

p_EGFR <- ggplot() +
  geom_point(data = filter(output_stats,protein == "EGFR"), aes(x = time,
                                                               y = mean_norm_actin_time,
                                                               group = condition,
                                                               color = condition,
                                                               shape = condition)) +
  geom_line(data = filter(output_stats,protein == "EGFR"), aes(x = time,
                                                                y = mean_norm_actin_time,
                                                                group = condition,
                                                                color = condition)) +
  geom_errorbar(data = filter(output_stats,protein == "EGFR"), aes(x = time,
                                                                    ymax = mean_norm_actin_time + sem_norm_actin_time,
                                                                    ymin = mean_norm_actin_time - sem_norm_actin_time,
                                                                    group = condition,
                                                                    color = condition)) + 
  ylab("EGFR/actin")

p_pAKT <- ggplot() +
  geom_point(data = filter(output_stats,protein == "pAKT"), aes(x = time,
                                                               y = mean_norm_actin_time,
                                                               group = condition,
                                                               color = condition,
                                                               shape = condition)) +
  geom_line(data = filter(output_stats,protein == "pAKT"), aes(x = time,
                                                               y = mean_norm_actin_time,
                                                               group = condition,
                                                               color = condition)) +
  geom_errorbar(data = filter(output_stats,protein == "pAKT"), aes(x = time,
                                                                   ymax = mean_norm_actin_time + sem_norm_actin_time,
                                                                   ymin = mean_norm_actin_time - sem_norm_actin_time,
                                                                   group = condition,
                                                                   color = condition)) + 
  ylab("pAKT/actin")

p_ppERK <- ggplot() +
  geom_point(data = filter(output_stats,protein == "ppERK"), aes(x = time,
                                                                y = mean_norm_actin_time,
                                                                group = condition,
                                                                color = condition,
                                                                shape = condition)) +
  geom_line(data = filter(output_stats,protein == "ppERK"), aes(x = time,
                                                               y = mean_norm_actin_time,
                                                               group = condition,
                                                               color = condition)) +
  geom_errorbar(data = filter(output_stats,protein == "ppERK"), aes(x = time,
                                                                   ymax = mean_norm_actin_time + sem_norm_actin_time,
                                                                   ymin = mean_norm_actin_time - sem_norm_actin_time,
                                                                   group = condition,
                                                                   color = condition)) + 
  ylab("ppERK/actin")



grid.arrange(p_pEGFR_EGFR + theme_white,
             p_pAKT + theme_white,
             p_ppERK + theme_white,
             nrow = 1, ncol = 3)

