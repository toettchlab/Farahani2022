library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
library(xlsx)
rm(list = ls())

##########################################################################

# enter parameters
t_start <- 50

##########################################################################

# gather files
setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/ERBB2 EGFR multiplex/2022_06_21_ERBB2_EGFR_mp/")

ZtSH2_names <-  list.files(pattern="*ZtSH2.csv");
ZtSH2_list <- lapply(ZtSH2_names, read.csv);
ZtSH2_names <- t(str_replace(ZtSH2_names,".csv",""));

VISH2_names <-  list.files(pattern="*VISH2.csv");
VISH2_list <- lapply(VISH2_names, read.csv);
VISH2_names <- t(str_replace(VISH2_names,".csv",""));

# create lists for dataframes
list_raw <- list();
list_abs_decrease <- list();
list_rate_of_change <- list();
list_half_life <- list();
list <- list();

##########################################################################

# data analysis

for (i in seq_along(ZtSH2_names)) {
  
  # raw values
  matrix_ZtSH2 <- ZtSH2_list[[i]];
  colnames(matrix_ZtSH2) <- str_replace_all(colnames(matrix_ZtSH2), "cell.", "");
  
  matrix_VISH2 <- VISH2_list[[i]];
  colnames(matrix_VISH2) <- str_replace_all(colnames(matrix_VISH2), "cell.", "");
  
  # melt raw reporter values matrix and classify observations by condition & date
  matrix_ZtSH2_melt <- melt(matrix_ZtSH2,id = c("time"));
  colnames(matrix_ZtSH2_melt) <- c("time", "cell", "reporter_raw");
  
  matrix_ZtSH2_melt <- mutate(matrix_ZtSH2_melt, 
                              condition = str_replace(str_extract(ZtSH2_names[[i]],"[:alnum:]*_"),"_",""), 
                              date = str_extract(ZtSH2_names[[i]],"[:digit:]{8}"), 
                              reporter = str_replace_all(str_extract(ZtSH2_names[[i]],"_[:alnum:]{4}2"),"_",""));
  
  matrix_VISH2_melt <- melt(matrix_VISH2,id = c("time"));
  colnames(matrix_VISH2_melt) <- c("time", "cell", "reporter_raw");
  
  matrix_VISH2_melt <- mutate(matrix_VISH2_melt, 
                              condition = str_replace(str_extract(VISH2_names[[i]],"[:alnum:]*_"),"_",""), 
                              date = str_extract(VISH2_names[[i]],"[:digit:]{8}"), 
                              reporter = str_replace_all(str_extract(VISH2_names[[i]],"_[:alnum:]{4}2"),"_",""));
  
  # combine matrices for ZtSH2 and VISH2
  matrix <- bind_rows(matrix_ZtSH2_melt,matrix_VISH2_melt, .id = NULL)
  
  # normalize time to addition of ligand
  matrix["time"] <- (matrix["time"] - t_start)/60;
  
  # normalize ZtSH2 readout
  
  matrix <- matrix %>%
    group_by(cell,reporter) %>%
    mutate(reporter_norm_max = reporter_raw / mean(reporter_raw[time <= 0]),
           reporter_abs_decrease = (reporter_norm_max - 1) * -100);
  
  list[[i]] <- matrix;
  
}

output <- bind_rows(list);

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

# mean trajectories vs time (mean across replicates)

dat <- output %>% 
  group_by(time,date,reporter) %>% 
  summarise(mean_reporter_abs_decrease = mean(reporter_abs_decrease),
            sd_reporter_abs_decrease = sd(reporter_abs_decrease))
dat <- dat %>% 
  group_by(time,reporter) %>% 
  summarise(mean = mean(mean_reporter_abs_decrease),
            sd = sd(mean_reporter_abs_decrease))


p1 <- ggplot() + 
  geom_ribbon(data=dat, aes(x = time, 
                            ymin = mean - sd,
                            ymax = mean + sd, 
                            group = reporter, 
                            fill = reporter), alpha = 0.2) +
  geom_line(data=dat, aes(x = time, y = mean, group = reporter,color = reporter)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(min(dat["time"]), 15)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-10,60)) + 
  labs(x = "time (min)", y = "RTK activity") 

p1 + theme_white

##########################################################################

# mean trajectories vs time (min/max normalized, mean across replicates)

dat <- filter(output, time <= 15) %>% 
  group_by(time,date,reporter) %>% 
  summarise(mean_reporter_abs_decrease = mean(reporter_abs_decrease),
            sd_reporter_abs_decrease = sd(reporter_abs_decrease))
dat <- dat %>% 
  group_by(date,reporter) %>% 
  mutate(norm_reporter_abs_decrease = (mean_reporter_abs_decrease - min(mean_reporter_abs_decrease)) / (max(mean_reporter_abs_decrease) - min(mean_reporter_abs_decrease)))
dat <- dat %>% 
  group_by(time,reporter) %>% 
  summarise(mean = mean(norm_reporter_abs_decrease),
            sd = sd(norm_reporter_abs_decrease))


p1 <- ggplot() + 
  geom_ribbon(data=dat, aes(x = time, 
                            ymin = mean - sd,
                            ymax = mean + sd, 
                            group = reporter, 
                            fill = reporter), alpha = 0.2) +
  geom_line(data=dat, aes(x = time, y = mean, group = reporter,color = reporter)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(min(dat["time"]), 15)) +
  labs(x = "time (min)", y = "RTK activity") 

p1 + theme_white

##########################################################################

# time to half of maximal response

dat <- filter(output, time <= 13 & time >= 0) %>% 
  group_by(cell,date,reporter) %>% 
  filter(reporter_abs_decrease > 0.5 * (max(reporter_abs_decrease) - min(reporter_abs_decrease))) %>% 
  mutate(t_half = min(time))

p1 <- ggplot(data=dat, aes(x = reporter, y = t_half, fill = reporter)) +
  geom_boxplot(colour = "black") +
  labs(x = "reporter", y = "half-life (min)")

p1 + theme_white

# run K-S tests
ks.test(filter(dat, reporter == "VISH2")[["t_half"]],
        filter(dat, reporter == "ZtSH2")[["t_half"]])





