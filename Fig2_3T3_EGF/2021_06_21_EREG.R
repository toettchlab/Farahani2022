# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
rm(list = ls())

##########################################################################

# enter parameters
frame_pre <- 2; # final frame before ligand addition
time_ligand <- 80; # ligand added after this time point

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/ERBB2:EGFR kinetics comparison/2022_06_21_EREG/");

ZtSH2_names <-  list.files(pattern="*2.csv");
ZtSH2_list <- lapply(ZtSH2_names, read.csv);
ZtSH2_names <- t(str_replace(ZtSH2_names,"_ZtSH2.csv",""));

RTK_names <- list.files(pattern = "*RTK.csv");
RTK_names <- t(RTK_names);
RTK_list <- lapply(RTK_names, read.csv);

EGFR_names <- list.files(pattern = "*EGFRcit.csv");
EGFR_names <- t(EGFR_names);
EGFR_list <- lapply(EGFR_names, read.csv);

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
  
  # ITAM-RTK expression levels
  matrix_RTK <- RTK_list[[i]];
  colnames(matrix_RTK) <- str_replace_all(colnames(matrix_RTK), "cell.", "");
  matrix_RTK <- melt(matrix_RTK);
  
  # EGFR-Citrine expression levels
  matrix_EGFR <- EGFR_list[[i]];
  colnames(matrix_EGFR) <- str_replace_all(colnames(matrix_EGFR), "cell.", "");
  matrix_EGFR <- melt(matrix_EGFR);
  
  # melt raw ZtSH2 values matrix and classify observations by condition & date
  matrix_melt <- melt(matrix,id = c("time"));
  colnames(matrix_melt) <- c("time", "cell", "ZtSH2_raw");
  
  matrix_melt <- mutate(matrix_melt, 
                        condition = str_replace(str_extract(ZtSH2_names[[i]],"[:alnum:]*_"),"_",""), 
                        date = str_extract(ZtSH2_names[[i]],"[:digit:]{8}"), 
                        ligand = str_replace_all(str_extract(ZtSH2_names[[i]],"_[:alpha:]*_"),"_",""));

  # incorporate expression levels into data matrix
  
  expression <- data.frame(matrix(ncol = 3, nrow = nrow(matrix_melt)));
  colnames(expression) <- c("RTK_expression", "EGFR_expression", "ZtSH2_expression");
  for (j in 1:nrow(matrix_melt)){
    expression[j,1] <- filter(matrix_RTK, variable == matrix_melt[j,"cell"])[2];
    expression[j,2] <- filter(matrix_EGFR, variable == matrix_melt[j,"cell"])[2];
    expression[j,3] <- filter(matrix_melt, time == 0 & cell == matrix_melt[j,"cell"])[3];
  }

  matrix_melt <- cbind(matrix_melt, expression);
  
  # calculate ZtSH2 response
  
  matrix_melt <- matrix_melt %>% 
    group_by(cell) %>% 
    mutate(ZtSH2_norm_max = ZtSH2_raw / mean(head(ZtSH2_raw,frame_pre)), 
           ZtSH2_abs_decrease = (ZtSH2_norm_max - 1) * -100);
  
  # normalize time to point of ligand addition
  matrix_melt["time"] <- (matrix_melt["time"] - time_ligand)/60;
  
  # add matrix to data frame
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

# mean cytoplasmic intensity vs. time

dat <- filter(output,
                condition == "EGFR" & 
                ligand == "EGF") %>%
  group_by(time,condition,ligand,date) %>%
  summarise(mean = mean(ZtSH2_abs_decrease))
dat <- dat %>% 
  group_by(time,condition,ligand) %>% 
  summarise(mean_rep = mean(mean),
            sd = sd(mean))


p1 <- ggplot() +
  geom_line(data = dat, aes(x = time, y = mean_rep, group = interaction(ligand,condition), color = interaction(ligand,condition))) +
  geom_ribbon(data = dat, aes(x = time,
                                       max = mean_rep + sd,
                                       min = mean_rep - sd,
                                       group = interaction(ligand,condition),
                                       fill = interaction(ligand,condition)),
              alpha = 1) + 
  scale_x_continuous(expand = c(0, 0), limits = c(min(dat["time"]),30)) + 
  ylab("RTK activity") + xlab("time (min)")

p1 + theme_white



             