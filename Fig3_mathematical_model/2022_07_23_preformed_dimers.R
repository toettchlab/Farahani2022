# import libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(reshape2)
rm(list = ls())

##########################################################################

# gather files

setwd("/Users/payamfarahani/Documents/Python_R/RTK BIOSENSORS/RTK_biosensor_modeling/2022_07_23");

py_names <- list.files(pattern="data_preformed_dimers.csv")
py <- read_csv(py_names[[1]])

##########################################################################

py["time"] <- py["time"] / 60

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

dat <- filter(py,EGFR_per_cell == 250000)

p1 <- ggplot() +
  geom_line(data = dat, aes(x = time,
                            y = ZtSH2_clearance,
                            group = relative_off_rate,
                            color = relative_off_rate)) + 
              # scale_colour_gradient(name = "count", trans = "log10",colours = terrain.colors(10)) + 
  scale_colour_viridis_c(trans = "log10") + 
  scale_x_continuous(expand = c(0,0), limits = c(-2,30)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.05,0.9)) +
  xlab("time (min)") + ylab("RTK activity")

p1 + theme_white

p1 <- ggplot() +
  geom_point(data = filter(dat, time == 0), aes(x = relative_off_rate,
                            y = frac,
                            group = relative_off_rate,
                            color = relative_off_rate)) + 
  # scale_colour_gradient(name = "count", trans = "log10",colours = terrain.colors(10)) + 
  scale_colour_viridis_c(trans = "log10") +  
  # scale_x_continuous(expand = c(0,0), limits = c(0,1800)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.05,0.9)) +
  xlab("relative off rate") + ylab("relative pre-formed dimers")

p1 + theme_white

# dat <- filter(output, ligand == "EREG" & RTK_expression < 1500) %>% 
#   group_by(time,condition,ligand) %>% 
#   summarise(mean_ZtSH2_norm_max = mean(ZtSH2_norm_max),
#             sd_ZtSH2_norm_max = sd(ZtSH2_norm_max),
#             mean_ZtSH2_abs_decrease = mean(ZtSH2_abs_decrease),
#             sd_ZtSH2_abs_decrease = sd(ZtSH2_abs_decrease))
# 
# p1 <- ggplot() +
#   geom_line(data = dat, aes(x = time,
#                             y = mean_ZtSH2_abs_decrease,
#                             group = interaction(condition,ligand),
#                             color = interaction(condition,ligand))) +
#   scale_x_continuous(expand = c(0,0), limits = c(0,1800)) +
#   scale_y_continuous(expand = c(0,0), limits = c(-0.2,0.6)) +
#   xlab("time (min)") + ylab("RTK activity")
# 
# p1 + theme_white

