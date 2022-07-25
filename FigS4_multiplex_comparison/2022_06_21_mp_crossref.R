library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
library(xlsx)
rm(list = ls())

##########################################################################

# gather files
setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/ERBB2 EGFR multiplex cross reference/2022_06_21/")

SH2_names <-  list.files(pattern="*SH2.csv");
SH2_list <- lapply(SH2_names, read.csv);
SH2_names <- t(str_replace(SH2_names,".csv",""));

# create lists for dataframes
list_raw <- list();
list_abs_decrease <- list();
list_rate_of_change <- list();
list_half_life <- list();
list <- list();

##########################################################################

# data analysis

for (i in seq_along(SH2_names)) {
  
  # raw values
  matrix_SH2 <- SH2_list[[i]];
  colnames(matrix_SH2) <- str_replace_all(colnames(matrix_SH2), "cell.", "");
  
  # melt raw reporter values matrix and classify observations by condition & date
  matrix_SH2_melt <- melt(matrix_SH2,id = c("time"));
  colnames(matrix_SH2_melt) <- c("time", "cell", "reporter_raw");
  
  matrix_SH2_melt <- mutate(matrix_SH2_melt,
                              experiment = str_replace(str_extract(SH2_names[[i]],"[:alnum:]*_"),"_",""),
                              receptor = str_replace_all(str_extract(SH2_names[[i]],"_[:alnum:]{4,5}_"),"_",""),
                              reporter = str_replace_all(str_extract(SH2_names[[i]],"_[:alnum:]{2}SH2"),"_",""),
                              date = str_extract(SH2_names[[i]],"[:digit:]{8}"));
  
  # normalize time to addition of ligand
  matrix_SH2_melt["time"] <- matrix_SH2_melt["time"]/60;
  
  # normalize ZtSH2 readout

  matrix_SH2_melt <- matrix_SH2_melt %>% 
    group_by(cell,experiment,receptor,reporter) %>% 
    mutate(reporter_norm_max = reporter_raw / mean(reporter_raw[time <= 0]), 
           reporter_abs_decrease = (reporter_norm_max - 1) * -100);
  
  list[[i]] <- matrix_SH2_melt;
  
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

# reporter_abs_decrease

dat <- filter(output, time <= 13) %>% 
  group_by(time,experiment,reporter,receptor) %>% 
  mutate(mean_reporter_abs_decrease = mean(reporter_abs_decrease),
         sd_reporter_abs_decrease = sd(reporter_abs_decrease))

p1 <- ggplot() + 
  geom_line(data=filter(dat, cell == "1"), aes(x = time, 
                                               y = mean_reporter_abs_decrease, 
                                               group = interaction(experiment,reporter,receptor), 
                                               color = interaction(experiment,reporter,receptor))) + 
  scale_x_continuous(expand = c(0, 0), limits = c(min(dat["time"]), 13)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(-0.1,0.6)) + 
  labs(x = "time (min)", y = "RTK activity") 

p1 + theme_white
##########################################################################

# mean trajectories vs time (mean across all cells)

dat <- filter(output, time <= 13) %>% 
  group_by(time,experiment,reporter,receptor) %>% 
  mutate(mean_reporter_abs_decrease = mean(reporter_abs_decrease),
            sd_reporter_abs_decrease = sd(reporter_abs_decrease))
dat <- dat %>% 
  group_by(experiment,reporter,receptor) %>% 
  mutate(norm_abs_decrease = (mean_reporter_abs_decrease - min(mean_reporter_abs_decrease)) / (max(mean_reporter_abs_decrease) - min(mean_reporter_abs_decrease)))

p1 <- ggplot() + 
  geom_line(data=filter(dat, cell == "1"), aes(x = time, 
                                               y = norm_abs_decrease, 
                                               group = interaction(experiment,reporter,receptor), 
                                               color = interaction(experiment,reporter,receptor))) + 
  scale_x_continuous(expand = c(0, 0), limits = c(min(dat["time"]), 13)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(-0.1,0.6)) + 
  labs(x = "time (min)", y = "RTK activity") 

p1 + theme_white

##########################################################################

# time to half of maximal response

dat <- filter(output, time <= 13 & time >= 0) %>% 
  group_by(cell,date,experiment,receptor,reporter) %>% 
  filter(reporter_abs_decrease > 0.5 * (max(reporter_abs_decrease) - min(reporter_abs_decrease))) %>% 
  mutate(t_half = min(time))

p1 <- ggplot(data=dat, aes(x = interaction(experiment,receptor,reporter), y = t_half, fill = interaction(experiment,receptor,reporter))) +
  geom_boxplot(colour = "black") +
  labs(x = "reporter", y = "half-life (min)")

p1 + theme_white





