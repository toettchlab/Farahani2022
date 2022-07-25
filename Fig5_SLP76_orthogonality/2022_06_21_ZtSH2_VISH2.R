library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
library(xlsx)
rm(list = ls())

##########################################################################

# enter parameters
frame_pre <- 3; # final frame before ligand addition
time_ligand <- 40; # ligand added after this time point

##########################################################################

# gather files
setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/SLP76 orthogonality/2022_06_21_ZtSH2_VISH2/")

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
  
  # normalize time to point of ligand addition
  matrix["time"] <- (matrix["time"] - time_ligand)/60;
  
  # calculate ZtSH2 response
  
  matrix <- matrix %>% 
    group_by(cell,reporter) %>% 
    mutate(reporter_norm_max = reporter_raw / mean(head(reporter_raw,frame_pre)), 
           reporter_abs_decrease = (reporter_norm_max - 1) * -100);

  list[[i]] <- matrix;

  print(str_c(ZtSH2_names[[i]]," done"))
  
}

output <- rbind_list(list);

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
        axis.line = element_line(colour="white"), 
        text = element_text(family = "Arial")) + 
  theme(aspect.ratio=1)

theme_white <-   theme_classic() +
  theme(legend.position="right",
        legend.background = element_rect(fill = 'white', colour = 'white'),
        legend.text = element_text(colour="black"),
        plot.background = element_rect(fill = 'white', colour = 'white'), 
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.text = element_text(colour="black"),
        axis.title = element_text(colour="black"),
        axis.line = element_line(colour="black"), 
        text = element_text(family = "Arial")) + 
  theme(aspect.ratio=1)

##########################################################################

# cd3e responses

dat <- filter(output, condition == "cd3e" & time <=10) %>%
  group_by(time,condition,reporter,date) %>%
  summarise(mean = mean(reporter_abs_decrease))
dat <- dat %>% 
  group_by(time,condition,reporter) %>% 
  summarise(mean_rep = mean(mean),
            sd = sd(mean))

p1 <- ggplot() +
  geom_ribbon(data=dat, aes(x = time, 
                            ymin = mean_rep - sd, 
                            ymax = mean_rep + sd,
                            group = reporter,
                            fill = reporter),
              alpha=0.2) + 
  geom_line(data=dat, aes(x = time, y = mean_rep,group = reporter, color=reporter)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(min(dat["time"]), 10)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-10, 60)) +
  labs(x = "time (min)", y = "RTK activity", title = "EGFR-CD3ε")

p1 + theme_white

##########################################################################

# SLP76 responses

dat <- filter(output, condition == "slp76" & time <=10) %>%
  group_by(time,condition,reporter,date) %>%
  summarise(mean = mean(reporter_abs_decrease))
dat <- dat %>% 
  group_by(time,condition,reporter) %>% 
  summarise(mean_rep = mean(mean),
            sd = sd(mean))

p1 <- ggplot() +
  geom_ribbon(data=dat, aes(x = time, 
                            ymin = mean_rep - sd, 
                            ymax = mean_rep + sd,
                            group = reporter,
                            fill = reporter),
              alpha=0.2) + 
  geom_line(data=dat, aes(x = time, y = mean_rep,group = reporter, color=reporter)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(min(dat["time"]), 10)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-10, 60)) +
  labs(x = "time (min)", y = "RTK activity", title = "EGFR-SLP76")

p1 + theme_white

##########################################################################

# boxplot
dat <- filter(output,time == 10)
x = c("slp76", "cd3e")
dat$condition <- factor(dat$condition,
                       levels = x)  

p1 <- ggplot(dat, aes(x = interaction(condition,reporter), y = reporter_abs_decrease, fill = reporter)) +
  geom_boxplot(colour = "black") +
  labs(x = "reporter", y = "norm. cyt. intensity") +
  scale_y_continuous(expand = c(0, 0), limits = c(-20, 75))

p1 + theme_white

