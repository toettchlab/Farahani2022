# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
rm(list = ls())

##########################################################################

# enter parameters
frame_pre <- 3; # final frame before ligand addition
time_ligand <- 40; # ligand added after this time point

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/CRISPR/2022_07_12/");

ZtSH2_names <-  list.files(pattern="*2.csv");
ZtSH2_list <- lapply(ZtSH2_names, read.csv);
ZtSH2_names <- t(str_replace(ZtSH2_names,"_ZtSH2.csv",""));

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
  matrix <- ZtSH2_list[[i]];
  colnames(matrix) <- str_replace_all(colnames(matrix), "cell.", "");

  # melt raw ZtSH2 values matrix and classify observations by condition & date
  matrix_melt <- melt(matrix,id = c("time"));
  colnames(matrix_melt) <- c("time", "cell", "ZtSH2_raw");
  
  matrix_melt <- mutate(matrix_melt, 
                        condition = str_replace(str_extract(ZtSH2_names[[i]],"[:alnum:]*_"),"_",""), 
                        date = str_extract(ZtSH2_names[[i]],"[:digit:]{8}"));

  # incorporate expression levels into data matrix
  
  expression <- data.frame(matrix(ncol = 1, nrow = nrow(matrix_melt)));
  colnames(expression) <- c("ZtSH2_expression");
  for (j in 1:nrow(matrix_melt)){
    expression[j,1] <- filter(matrix_melt, time == 0 & cell == matrix_melt[j,"cell"])[3];
  }

  # calculate ZtSH2 response
  
  matrix_melt <- cbind(matrix_melt, expression);
  
  matrix_melt <- matrix_melt %>% 
    group_by(cell) %>% 
    mutate(ZtSH2_norm_max = ZtSH2_raw / ZtSH2_raw[time <= 0], 
           ZtSH2_abs_decrease = (ZtSH2_norm_max - 1) * -100);
  
  # calculate integrated ZtSH2 response (AUC)
  
  matrix_melt <- matrix_melt %>% 
    group_by(cell) %>% 
    mutate(ZtSH2_AUC = cumsum(ZtSH2_abs_decrease))
  
  # normalize time to ligand addition 
  matrix_melt["time"] <- (matrix_melt["time"] - time_ligand)/60;

  list[[i]] <- matrix_melt;
  
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

# mean activity vs time

# mean across replicates
dat <- filter(output, time <= 6) %>%
  group_by(time,condition,date) %>%
  summarise(mean = mean(ZtSH2_abs_decrease))

dat <- dat %>%
  group_by(time,condition) %>%
  summarise(mean_rep = mean(mean),
            sd = sd(mean))

p1 <- ggplot() +
  geom_line(data = dat, aes(x = time,
                            y = mean_rep,
                            group = condition,
                            color = condition)) +
  geom_ribbon(data = dat, aes(x = time,
                              max = mean_rep + sd,
                              min = mean_rep - sd,
                              group = condition,
                              fill = condition),
              alpha = 1) +
  # scale_x_continuous(expand = c(0,0), limits = c(min(dat["time"]),10)) + 
  xlab("time (min)") + ylab("RTK activity")

p1 + theme_white

##########################################################################

# mean activity vs time (mean across all cells)

# mean across replicates
dat <- filter(output, time <= 6) %>%
  group_by(time,condition) %>%
  summarise(mean = mean(ZtSH2_abs_decrease),
            sd = sd(ZtSH2_abs_decrease))


p1 <- ggplot() +
  geom_line(data = dat, aes(x = time,
                            y = mean,
                            group = condition,
                            color = condition)) +
  geom_ribbon(data = dat, aes(x = time,
                              max = mean + sd,
                              min = mean - sd,
                              group = condition,
                              fill = condition),
              alpha = 0.2) +
  # scale_x_continuous(expand = c(0,0), limits = c(min(dat["time"]),10)) + 
  xlab("time (min)") + ylab("RTK activity")

# p1 + theme_white

##########################################################################

# boxplot
dat <- filter(output,time == 1)

x = c("WT", "mNG2")

dat$condition <- factor(dat$condition,
                       levels = x)  

p1 <- ggplot(dat, aes(x = condition, y = ZtSH2_abs_decrease, fill = condition)) +
  geom_boxplot(colour = "black") +
  labs(x = "reporter", y = "norm. cyt. intensity")

p1 + theme_white

# run K-S tests
ks.test(filter(dat, condition == "WT")[["ZtSH2_abs_decrease"]],
        filter(dat, condition == "mNG2")[["ZtSH2_abs_decrease"]])

