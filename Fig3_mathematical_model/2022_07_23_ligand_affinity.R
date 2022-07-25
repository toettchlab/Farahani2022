# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
rm(list = ls())

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/RTK_biosensor_modeling/2022_07_23/");

py_names <- list.files(pattern="data_ligand_affinity.csv")
py <- read_csv(py_names[[1]])

##########################################################################

py["time"] <- py["time"] / 60;

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

dat <- filter(py,ligand_off_rate == 1 & dimer_off_rate == 100)

p1 <- ggplot() +
  geom_line(data = dat, aes(x = time,
                            y = ZtSH2_clearance,
                            group = ligand_dose,
                            color = ligand_dose)) + 
              # scale_colour_gradient(name = "count", trans = "log10",colours = terrain.colors(10)) + 
  scale_colour_viridis_c(trans = "log10") + 
  scale_x_continuous(expand = c(0,0), limits = c(-2,30)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.1,0.8)) +
  xlab("time (min)") + ylab("EGFR activity")

p1 + theme_white




