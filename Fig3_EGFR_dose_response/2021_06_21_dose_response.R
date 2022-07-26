# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
library(xlsx)
rm(list = ls())

##########################################################################

# enter parameters
frame_pre <- 2; # final frame before ligand addition
time_ligand <- 2; # ligand added after this time point

MW_EGF <- 6; # kDa
MW_EREG <- 5.4;
MW_EPGN <- 6;

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/EGFR dose response/2022_06_21");

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
  
  # melt raw ZtSH2 values matrix and classify observations by condition & date
  matrix_melt <- melt(matrix,id = c("time"));
  colnames(matrix_melt) <- c("time", "cell", "ZtSH2_raw");
  
  matrix_melt <- mutate(matrix_melt, 
                        condition = str_replace(str_extract(ZtSH2_names[[i]],"[:alnum:]*_"),"_",""), 
                        date = str_extract(ZtSH2_names[[i]],"[:digit:]{8}"), 
                        dose = str_c(str_replace_all(str_extract(ZtSH2_names[[i]],"_[:digit:]{1,4}_"),"_",""),"ng/mL", sep = " "),
                        dose_nM = as.numeric(str_replace_all(str_extract(ZtSH2_names[[i]],"_[:digit:]{1,4}_"),"_","")));
  
  if (matrix_melt["condition"] == "EGF"){
    matrix_melt["dose_nM"] <- matrix_melt["dose_nM"] / MW_EGF;
  } else if (matrix_melt["condition"] == "EREG"){
    matrix_melt["dose_nM"] <- matrix_melt["dose_nM"] / MW_EREG;
  } else if (matrix_melt["condition"] == "EPGN"){
    matrix_melt["dose_nM"] <- matrix_melt["dose_nM"] / MW_EPGN;
  }

  # incorporate expression levels into data matrix
  
  expression <- data.frame(matrix(ncol = 2, nrow = nrow(matrix_melt)));
  colnames(expression) <- c("RTK_expression", "ZtSH2_expression");
  for (j in 1:nrow(matrix_melt)){
    expression[j,1] <- filter(matrix_RTK, variable == matrix_melt[j,"cell"])[2];
    expression[j,2] <- filter(matrix_melt, time == 0 & cell == matrix_melt[j,"cell"])[3];
  }

  matrix_melt <- cbind(matrix_melt, expression);
  
  # calculate ZtSH2 response
  
  matrix_melt <- matrix_melt %>% 
    group_by(cell) %>% 
    mutate(ZtSH2_norm_max = ZtSH2_raw / mean(head(ZtSH2_raw,frame_pre)), 
           ZtSH2_abs_decrease = (ZtSH2_norm_max - 1) * -100);
  
  # normalize time to point of ligand addition
  matrix_melt["time"] <- matrix_melt["time"] - time_ligand;

  # add matrix to data frame
  list[[i]] <- matrix_melt;
  
  print(str_c(ZtSH2_names[[i]]," done"))
  
}

output <- rbind_list(list);

##########################################################################

# statistics for output dataframe

# statistics of cells pooled together
output_stats <- output %>% 
  group_by(time,condition,dose) %>% 
  summarise(mean_ZtSH2_norm_max = mean(ZtSH2_norm_max),
            sd_ZtSH2_norm_max = sd(ZtSH2_norm_max),
            mean_ZtSH2_abs_decrease = mean(ZtSH2_abs_decrease),
            sd_ZtSH2_abs_decrease = sd(ZtSH2_abs_decrease))

# statistics of replicates
output_stats_rep <- output %>% 
  group_by(time,condition,dose,date) %>% 
  summarise(mean_ZtSH2_norm_max = mean(ZtSH2_norm_max),
            mean_ZtSH2_abs_decrease = mean(ZtSH2_abs_decrease))

output_stats_rep <- output_stats_rep %>% 
  group_by(time,condition,dose) %>% 
  summarise(mean_ZtSH2_abs_decrease_rep = mean(mean_ZtSH2_abs_decrease),
  sd_ZtSH2_abs_decrease_rep = sd(mean_ZtSH2_abs_decrease))

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

# EGF

dat <- filter(output_stats_rep, condition == "EGF" | (condition == "EGF" & dose == "0 ng/mL"))
x <- c("100 ng/mL","20 ng/mL","5 ng/mL","2 ng/mL","02 ng/mL","0 ng/mL")
dat$dose <- factor(dat$dose,levels = x)  

p1 <- ggplot() +
  geom_line(data = dat, aes(x = time,
                            y = mean_ZtSH2_abs_decrease_rep,
                            group = dose,
                            color = dose)) +
  geom_ribbon(data = dat, aes(x = time,
                              max = mean_ZtSH2_abs_decrease_rep + sd_ZtSH2_abs_decrease_rep,
                              min = mean_ZtSH2_abs_decrease_rep - sd_ZtSH2_abs_decrease_rep,
                              group = dose,
                              fill = dose),
              alpha = 0.2) +
  scale_x_continuous(expand = c(0,0), limits = c(-2,60)) +
  scale_y_continuous(expand = c(0,0), limits = c(-10,60)) +
  xlab("time (min)") + ylab("RTK activity")

p1 + theme_white

##########################################################################

# EREG

dat <- filter(output_stats_rep, condition == "EREG" | (condition == "EGF" & dose == "0 ng/mL"))
x <- c("1000 ng/mL","500 ng/mL","100 ng/mL","20 ng/mL","10 ng/mL","0 ng/mL")
dat$dose <- factor(dat$dose,levels = x)  

p1 <- ggplot() +
  geom_line(data = dat, aes(x = time,
                            y = mean_ZtSH2_abs_decrease_rep,
                            group = dose,
                            color = dose)) +
  geom_ribbon(data = dat, aes(x = time,
                              max = mean_ZtSH2_abs_decrease_rep + sd_ZtSH2_abs_decrease_rep,
                              min = mean_ZtSH2_abs_decrease_rep - sd_ZtSH2_abs_decrease_rep,
                              group = dose,
                              fill = dose),
              alpha = 0.2) +
  scale_x_continuous(expand = c(0,0), limits = c(-2,60)) +
  scale_y_continuous(expand = c(0,0), limits = c(-10,60)) +
  xlab("time (min)") + ylab("RTK activity")

p1 + theme_white

##########################################################################

# EPGN

dat <- filter(output_stats_rep, condition == "EPGN" | (condition == "EGF" & dose == "0 ng/mL"))
x <- c("1000 ng/mL","500 ng/mL","100 ng/mL","20 ng/mL","0 ng/mL")
dat$dose <- factor(dat$dose,levels = x)  

p1 <- ggplot() +
  geom_line(data = dat, aes(x = time,
                            y = mean_ZtSH2_abs_decrease_rep,
                            group = dose,
                            color = dose)) +
  geom_ribbon(data = dat, aes(x = time,
                              max = mean_ZtSH2_abs_decrease_rep + sd_ZtSH2_abs_decrease_rep,
                              min = mean_ZtSH2_abs_decrease_rep - sd_ZtSH2_abs_decrease_rep,
                              group = dose,
                              fill = dose),
              alpha = 0.2) +
  scale_x_continuous(expand = c(0,0), limits = c(-2,60)) +
  scale_y_continuous(expand = c(0,0), limits = c(-10,60)) +
  xlab("time (min)") + ylab("RTK activity")

p1 + theme_white

