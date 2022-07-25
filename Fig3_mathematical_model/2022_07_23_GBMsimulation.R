# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
rm(list = ls())

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/RTK_biosensor_modeling/2022_07_23/");

py_names <- list.files(pattern="data_GBMsimulation.csv")
py <- read_csv(py_names[[1]])

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

dat <- py
dat["time"] <- dat["time"]/60

p1 <- ggplot() +
  geom_line(data = dat, aes(x = time,
                            y = ZtSH2_clearance,
                            group = interaction(dimer_off_rate),
                            color = interaction(dimer_off_rate))) + 
  scale_x_continuous(expand = c(0,0), limits = c(-2,30)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.05,0.4)) +
  xlab("time (min)") + ylab("RTK activity")

p1 + theme_white




