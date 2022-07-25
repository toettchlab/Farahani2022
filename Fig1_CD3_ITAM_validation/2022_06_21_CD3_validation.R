# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
rm(list = ls())

##########################################################################

# enter parameters
frame_pre <- 10; # final frame before ligand addition
time_ligand <- 9; # ligand added after this time point

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/CD3 ITAM Validation/2022_06_21_CD3_validation");

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
  
  condition <- data.frame(matrix(ncol = 1, nrow = nrow(matrix_melt)));
  colnames(condition) <- c("condition");
  condition$condition <- str_replace(str_extract(ZtSH2_names[[i]],"[:alnum:]*_"),"_","");

  date <- data.frame(matrix(ncol = 1, nrow = nrow(matrix_melt)));
  colnames(date) <- c("date");
  date$date <- str_extract(ZtSH2_names[[i]],"[:digit:]{8}");

  # normalize time to point of ligand addition
  matrix_melt["time"] <- matrix_melt["time"] - time_ligand;
  
  # incorporate expression levels into data matrix
  
  expression <- data.frame(matrix(ncol = 1, nrow = nrow(matrix_melt)));
  colnames(expression) <- c("ZtSH2_expression");
  for (j in 1:nrow(matrix_melt)){
    expression[j,1] <- filter(matrix_melt, time == 0 & cell == matrix_melt[j,"cell"])[3];
  }

  matrix_melt <- cbind(matrix_melt, expression);

  # calculate ZtSH2 response
  
  matrix_melt <- matrix_melt %>% 
    group_by(cell) %>% 
    mutate(ZtSH2_norm_max = ZtSH2_raw / mean(head(ZtSH2_raw,frame_pre)), 
           ZtSH2_abs_decrease = (ZtSH2_norm_max - 1) * -100);

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

# mean EGFR activity vs. time (heatmap)

dat <- output %>% 
  group_by(time,condition) %>% 
  mutate(mean_ZtSH2_abs_decrease = mean(ZtSH2_abs_decrease))
x = c("epsilon","delta", "gamma", "zeta3", "zeta2","zeta1", "EGFR", "none")
dat$condition <- factor(dat$condition,
                        levels = x)   

p1 <- ggplot() +
  geom_tile(data = dat, aes(x = time,
                            y = condition,
                            fill = mean_ZtSH2_abs_decrease)) + 
  scale_fill_viridis_c()

p1 + theme_white

##########################################################################

# boxplot

dat <- filter(output, time == "10" | time == "60") %>%
  unite("cond", c(condition,time), sep = "_", remove = FALSE)


# x = c("none_20", "none_60","EGFR_20", "EGFR_60","zeta1_20", "zeta1_60","zeta2_20", "zeta2_60","zeta3_20", "zeta3_60","gamma_20", "gamma_60","delta_20", "delta_60","epsilon_20", "epsilon_60")
x = c("none_10", "none_60","EGFR_10", "EGFR_60","zeta1_10", "zeta1_60","zeta2_10", "zeta2_60","zeta3_10", "zeta3_60","gamma_10", "gamma_60","delta_10", "delta_60","epsilon_10", "epsilon_60")
dat$cond <- factor(dat$cond,
                       levels = x)                        
                                        
p1 <- ggplot(data = dat, 
             aes(x = cond, 
                 y = ZtSH2_abs_decrease, 
                 fill = condition)) +
  geom_boxplot(colour = "black") +
  labs(x = "reporter", y = "RTK activity (t = 1 h)")
  
p1 + theme_classic() +
  theme(legend.position="right",
        legend.background = element_rect(fill = 'white', colour = 'white'),
        legend.text = element_text(colour="black"),
        plot.background = element_rect(fill = 'white', colour = 'white'), 
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.text = element_text(colour="black"),
        axis.title = element_text(colour="black"),
        axis.line = element_line(colour="black"), 
        text = element_text(family = "Arial"))

# run K-S tests
ks.test(filter(output, time == "10" & condition == "none")[["ZtSH2_abs_decrease"]],
        filter(output, time == "60" & condition == "epsilon")[["ZtSH2_abs_decrease"]])

# check number of cells analyzed for each independent experiment
num <- output %>% 
  group_by(condition,date) %>% 
  filter(time == 0 & cell == max(cell))
