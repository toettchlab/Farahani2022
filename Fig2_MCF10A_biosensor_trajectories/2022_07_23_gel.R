# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
rm(list = ls())

##########################################################################

# enter parameters
frame_pre <- 2; # final frame before ligand addition
time_ligand <- 30; # ligand added after this time point
time_end <- 1800; # final time point for analysis

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/MCF10A gel/2022_06_07/");

ZtSH2_names <-  list.files(pattern="*ZtSH2.csv");
ZtSH2_list <- lapply(ZtSH2_names, read.csv);
ZtSH2_names <- t(str_replace(ZtSH2_names,"_ZtSH2.csv",""));

KTRcyt_names <-  list.files(pattern="*ErkKTR_cyt.csv");
KTRcyt_list <- lapply(KTRcyt_names, read.csv);
KTRcyt_names <- t(str_replace(KTRcyt_names,"_ErkKTR_cyt.csv",""));

KTRnuc_names <-  list.files(pattern="*ErkKTR_nuc.csv");
KTRnuc_list <- lapply(KTRnuc_names, read.csv);
KTRnuc_names <- t(str_replace(KTRnuc_names,"_ErkKTR_nuc.csv",""));

class_names <- list.files(pattern = "*class.csv");
class_names <- t(class_names);
class_list <- lapply(class_names, read.csv);

# create lists for dataframes
list <- list();

##########################################################################

# data analysis
for (i in seq_along(ZtSH2_names)) {
  
  # raw ZtSH2 values
  matrix <- ZtSH2_list[[i]];
  colnames(matrix) <- str_replace_all(colnames(matrix), "cell.", "");
  
  # raw KTR values
  matrix_KTRnuc <- KTRnuc_list[[i]];
  colnames(matrix_KTRnuc) <- str_replace_all(colnames(matrix_KTRnuc), "cell.", "");
  
  matrix_KTRcyt <- KTRcyt_list[[i]];
  colnames(matrix_KTRcyt) <- str_replace_all(colnames(matrix_KTRcyt), "cell.", "");
  
  # cell classifications
  matrix_class <- class_list[[i]];
  colnames(matrix_class) <- str_replace_all(colnames(matrix_class), "cell.", "");
  matrix_class <- melt(matrix_class);
  
  # melt raw value matrices and classify observations by condition & date
  matrix_melt <- melt(matrix,id = c("time"));
  colnames(matrix_melt) <- c("time", "cell", "ZtSH2_raw");
  matrix_melt <- mutate(matrix_melt, 
                        condition = str_replace(str_extract(ZtSH2_names[[i]],"[:alnum:]*_"),"_",""), 
                        date = str_extract(ZtSH2_names[[i]],"[:digit:]{8}"));
  matrix_KTRcyt <- melt(matrix_KTRcyt,id = c("time"));
  colnames(matrix_KTRcyt) <- c("time", "cell", "KTRcyt");
  matrix_KTRcyt <- mutate(matrix_KTRcyt, 
                        condition = str_replace(str_extract(KTRcyt_names[[i]],"[:alnum:]*_"),"_",""), 
                        date = str_extract(KTRcyt_names[[i]],"[:digit:]{8}"));
  matrix_KTRnuc <- melt(matrix_KTRnuc,id = c("time"));
  colnames(matrix_KTRnuc) <- c("time", "cell", "KTRnuc");
  matrix_KTRnuc <- mutate(matrix_KTRnuc, 
                          condition = str_replace(str_extract(KTRnuc_names[[i]],"[:alnum:]*_"),"_",""), 
                          date = str_extract(KTRnuc_names[[i]],"[:digit:]{8}"));

  # incorporate cell classification into data matrix
  class <- data.frame(matrix(ncol = 1, nrow = nrow(matrix_melt)));
  colnames(class) <- c("class");
  for (j in 1:nrow(matrix_melt)){
    class[j,1] <- filter(matrix_class, variable == matrix_melt[j,"cell"])[2];
  }
  
  matrix_melt <- cbind(matrix_melt, class);
  
  matrix_melt <- full_join(matrix_melt, matrix_KTRcyt);
  matrix_melt <- full_join(matrix_melt,matrix_KTRnuc);
  
  class(matrix_melt$cell) = "Integer"
  
  # calculate KTR C/N ratio
  
  matrix_melt <- matrix_melt %>% 
    group_by(cell,time) %>% 
    mutate(KTR_CN = KTRcyt/KTRnuc);

  # calculate ZtSH2 response
  
  matrix_melt <- matrix_melt %>% 
    group_by(cell) %>% 
    mutate(ZtSH2_norm_max = ZtSH2_raw / mean(head(ZtSH2_raw,frame_pre)), 
           ZtSH2_abs_decrease = (ZtSH2_norm_max - 1) * -1);
  
  # calculate integrated ZtSH2 response (AUC)
  
  matrix_melt <- matrix_melt %>% 
    group_by(cell) %>% 
    mutate(ZtSH2_AUC = cumsum(ZtSH2_abs_decrease))

  # normalize time to ligand addition
  matrix_melt["time"] <- (matrix_melt["time"]-time_ligand) / 60;
  
  # add matrix to list
  list[[i]] <- matrix_melt;
  
  print(str_c(ZtSH2_names[[i]]," done"))
  
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

# single-cell traces (EGFR pYtag)

dat <- filter(output,time <= 16) %>% 
  group_by(condition,date,cell) %>% 
  mutate(ZtSH2_norm = ZtSH2_abs_decrease/max(ZtSH2_abs_decrease), 
         KTR_norm = (KTR_CN - min(KTR_CN)) / (max(KTR_CN) - min(KTR_CN)))
dat <- arrange(dat,class)

p1 <- ggplot() +
  geom_tile(data = dat, aes(x = time,
                            y = interaction(cell,date,condition,class),
                            fill = ZtSH2_norm)) + 
  scale_fill_viridis_c()

p1 + theme_white

##########################################################################

# single-cell traces (ErkKTR)

dat <- filter(output,time <= 16) %>% 
  group_by(condition,date,cell) %>% 
  mutate(ZtSH2_norm = ZtSH2_abs_decrease/max(ZtSH2_abs_decrease), 
         KTR_norm = (KTR_CN - min(KTR_CN)) / (max(KTR_CN) - min(KTR_CN)))
dat <- arrange(dat,class)

p1 <- ggplot() +
  geom_tile(data = dat, aes(x = time,
                            y = interaction(cell,date,condition,class),
                            fill = KTR_norm)) + 
  scale_fill_viridis_c()

p1 + theme_white

##########################################################################

# mean activity vs time

dat <- filter(output,time <= 16) %>% 
  group_by(class,time) %>% 
  mutate(mean_ZtSH2 = mean(ZtSH2_abs_decrease),mean_KTR = mean(KTR_CN))

p1 <- ggplot() +
  geom_tile(data = dat, aes(x = time,
                            y = class,
                            fill = mean_ZtSH2)) + 
  scale_fill_viridis_c()


p1 + theme_white

p1 <- ggplot() +
  geom_tile(data = dat, aes(x = time,
                            y = class,
                            fill = mean_KTR)) + 
  scale_fill_viridis_c()

p1 + theme_white

##########################################################################

# time-to-half-max

dat_ZtSH2 <- filter(output,time <= 16) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_norm = ZtSH2_abs_decrease/max(ZtSH2_abs_decrease), 
         reporter = "ZtSH2")
dat_ZtSH2 <- arrange(dat_ZtSH2,class)
dat_ZtSH2 <- filter(dat_ZtSH2, reporter_norm >= 0.5) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_half = min(time))

dat_KTR <- filter(output,time <= 16) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_norm = (KTR_CN - min(KTR_CN)) / (max(KTR_CN) - min(KTR_CN)), 
         reporter = "KTR")
dat_KTR <- arrange(dat_KTR,class)
dat_KTR <- filter(dat_KTR, reporter_norm >= 0.5) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_half = min(time))

dat <- bind_rows(dat_ZtSH2,dat_KTR)
x = c("ZtSH2", "KTR")

dat$reporter <- factor(dat$reporter,
                   levels = x) 

p1 <- ggplot() + 
  geom_boxplot(data = dat, aes(x = interaction(class,reporter), 
                               y = reporter_half, 
                               group = interaction(class,reporter), 
                               fill = interaction(class,reporter)))
  # scale_y_continuous(expand = c(0, 0), limits = c(0,17))

p1 + theme_white

# run K-S tests
ks.test(filter(dat, reporter == "KTR" & class == 1)[["reporter_half"]],
        filter(dat, reporter == "KTR" & class == 2)[["reporter_half"]])


##########################################################################

# calculate time delay between pYtag and ErkKTR (time-to-half-max boxplot)

dat_ZtSH2 <- filter(output,time <= 16) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_norm = ZtSH2_abs_decrease/max(ZtSH2_abs_decrease), 
         reporter = "ZtSH2")
dat_ZtSH2 <- arrange(dat_ZtSH2,class)
dat_ZtSH2 <- filter(dat_ZtSH2, reporter_norm >= 0.5) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_half = min(time))

dat_KTR <- filter(output,time <= 16) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_norm = (KTR_CN - min(KTR_CN)) / (max(KTR_CN) - min(KTR_CN)), 
         reporter = "KTR")
dat_KTR <- arrange(dat_KTR,class)
dat_KTR <- filter(dat_KTR, reporter_norm >= 0.5) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_half = min(time))

dat <- bind_rows(dat_ZtSH2,dat_KTR)
x = c("ZtSH2", "KTR")

dat$reporter <- factor(dat$reporter,
                       levels = x) 

dat <- dat %>% 
  group_by(condition,date,cell,reporter) %>% 
  filter(time == min(time))
dat <- dat %>% 
  group_by(condition,date,cell) %>% 
  mutate(delay = -reporter_half[reporter == "ZtSH2"] + reporter_half[reporter == "KTR"]);

p1 <- ggplot() + 
  geom_boxplot(data = dat, aes(x = 1, y = delay)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,10))

p1 + theme_white

##########################################################################

# calculate time delay between pYtag and ErkKTR (time-to-half-max dot plot)

dat_ZtSH2 <- filter(output,time <= 16) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_norm = ZtSH2_abs_decrease/max(ZtSH2_abs_decrease), 
         reporter = "ZtSH2")
dat_ZtSH2 <- arrange(dat_ZtSH2,class)
dat_ZtSH2 <- filter(dat_ZtSH2, reporter_norm >= 0.5) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_half = min(time))

dat_KTR <- filter(output,time <= 16) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_norm = (KTR_CN - min(KTR_CN)) / (max(KTR_CN) - min(KTR_CN)), 
         reporter = "KTR")
dat_KTR <- arrange(dat_KTR,class)
dat_KTR <- filter(dat_KTR, reporter_norm >= 0.5) %>% 
  group_by(condition,date,cell) %>% 
  mutate(reporter_half = min(time))

dat <- bind_rows(dat_ZtSH2,dat_KTR)
x = c("ZtSH2", "KTR")

dat$reporter <- factor(dat$reporter,
                       levels = x) 

p1 <- ggplot() + 
  geom_point(data = dat, aes(x = reporter, y = reporter_half, group = interaction(cell,condition,date,class))) + 
  geom_line(data = dat, aes(x = reporter, y = reporter_half, group = interaction(cell,condition,date,class)))

p1 + theme_white

# run K-S tests
ks.test(filter(dat, reporter == "KTR")[["reporter_half"]],
        filter(dat, reporter == "ZtSH2")[["reporter_half"]])



