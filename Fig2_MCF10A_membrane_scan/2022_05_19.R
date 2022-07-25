# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
rm(list = ls())

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/MCF10A gel (membranes)/2022_05_19");

names <-  list.files(pattern="quantification.csv");
list <- lapply(names, read.csv);
names <- t(str_replace(names,".csv",""));

output <- list[[1]]

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

dat <- filter(output, type == "apical" & time == -1)
p1 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,200)) + 
  theme_white
p1

dat <- filter(output, type == "lateral" & time == -1)
p2 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(20,160)) + 
  theme_white
p2

dat <- filter(output, type == "apical" & time == 6)
p3 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,200)) + 
  theme_white
p3

dat <- filter(output, type == "lateral" & time == 6)
p4 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(20,100)) + 
  theme_white
p4

dat <- filter(output, type == "apical" & time == 11)
p5 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,200)) + 
  theme_white
p5

dat <- filter(output, type == "lateral" & time == 11)
p6 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) +
  scale_y_continuous(expand = c(0, 0), limits = c(20,100)) + 
  theme_white
p6


# grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 3, ncol = 2)

##########################################################################

mean_ZtSH2 <- mean(filter(output,type == "lateral" & channel == "ZtSH2")[,"intensity"])
mean_EGFR <- mean(filter(output,type == "lateral" & channel == "EGFR")[,"intensity"])

data <- output

data[data$channel == "ZtSH2","intensity"] <- data[data$channel == "ZtSH2","intensity"] / mean_ZtSH2
data[data$channel == "EGFR","intensity"] <- data[data$channel == "EGFR","intensity"] / mean_EGFR

# y axis limits
y_ap <- scale_y_continuous(expand = c(0, 0), limits = c(0,3.5))
y_lat <- scale_y_continuous(expand = c(0, 0), limits = c(0.4,2.5))

dat <- filter(data, type == "apical" & time == -1)
p1 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  y_ap + 
  theme_white
dat <- filter(data, type == "lateral" & time == -1)
p2 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  y_lat + 
  theme_white
dat <- filter(data, type == "apical" & time == 6)
p3 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  y_ap + 
  theme_white
dat <- filter(data, type == "lateral" & time == 6)
p4 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  y_lat + 
  theme_white
dat <- filter(data, type == "apical" & time == 11)
p5 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) + 
  y_ap + 
  theme_white
dat <- filter(data, type == "lateral" & time == 11)
p6 <- ggplot() + 
  geom_line(data = dat, aes(x = distance, y = intensity, group = channel, color = channel)) +
  y_lat + 
  theme_white


grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 3, ncol = 2)



