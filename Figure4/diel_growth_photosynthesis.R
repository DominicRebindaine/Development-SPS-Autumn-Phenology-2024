###########################################################################################################################################
# R Script for: 
# Figure 4 of the study:
# Development underpins the summer solstice reversal of the effects of climate warming on the autumn phenology of European beech
# Rebindaine et al. (2024)
# Dominic Rebindaine
# Last edited: 04/11/2024
###########################################################################################################################################

##### Using data extracted from Zweifel et al. (2021) "Why trees grow at night" Figure 2e to plot diel growth patterns of Fagus sylvatica
##### And using data extracted from Urban et al. (2014) "Impact of elevated CO2 concentration on dynamics of leaf photosynthesis in Fagus sylvatica is modulated by sky conditions" Figure 2a to plot diel photosynthesis patterns of Fagus sylvatica

library(tidyverse)
library(zoo)
library(slider)
library(patchwork)

input.path = '/Users/dominic/Documents/Crowther/r/buds/Data'

output.path = '/Users/dominic/Documents/Crowther/r/buds/figures/diurnal'

plotTheme1 = theme(
  legend.position   = "right",
  legend.background = element_rect(fill='white', linewidth=0.5, linetype="solid"),
  legend.text       = element_text(color="black", size = 14),
  panel.background  = element_blank(),
  axis.text         = element_text(colour = "black", size = 16),
  axis.title        = element_text(colour = "black", size = 18),
  panel.border      = element_rect(colour = "black", fill=NA),
  axis.line         = element_line(color = "black"),
  strip.background  = element_rect(fill=NA),
  strip.text        = element_text(colour = 'black'),
  plot.title        = element_text(face="bold",hjust = 0.5))


###### Hourly probability for growth
# data from Zweifel et al. (2021) figure 2e - percentage - number of hours with growth relative to total number of hours in the stem growth period
dat.e = read.csv(paste(input.path,"Fig2e_Zweifel_data.csv", sep = '/'))

# round x (hour)
dat.e$x = round(dat.e$x, 0)
# round y (%)
dat.e$y = round(dat.e$y, 2)

#Plot
dat.e %>% ggplot() + 
  geom_line(aes(x,y)) +
  xlab('Time of day (h)')+ylab('Probability for growth (%)') +
  plotTheme1

## Get relative values of y
dat.e$rel = dat.e$y/max(dat.e$y)*100

#Plot
dat.e %>% ggplot() + 
  geom_line(aes(x,rel)) +
  xlab('Time of day (h)')+ylab('Relative probability for growth (%)') +
  plotTheme1

  
###### Hourly carbon assimilation rate
# data from figure 4a of Urban et al 2014 
#"Impact of elevated CO2 concentration on dynamics of leaf photosynthesis in Fagus sylvatica is modulated by sky conditions"
dat.a = read.csv(paste(input.path,"Fig4a_Urban_data.csv", sep = '/'))

# if assimilation rate is -ve then force to 0
dat.a$y[dat.a$y < 0] = 0

# round x (hour)
dat.a$x = round(dat.a$x, 0)
# round y (%)
dat.a$y = round(dat.a$y, 2)

## Get relative values of y
dat.a$rel = dat.a$y/max(dat.a$y)*100

#### Plot relative rates of carbon assimilation and relative probability for growth on same figure (Fig. 4a)
# Step 1: Remove the 'y' column from both data frames
dat.a = dat.a[, c("x", "rel")]
dat.e = dat.e[, c("x", "rel")]

# Step 2: Merge the data frames by 'x', keeping all values from 'dat.e' (outer join)
dat <- merge(dat.e, dat.a, by = "x", all.x = TRUE, suffixes = c("_e", "_a"))

# Interpolate missing values in rel_a
dat$rel_a = na.approx(dat$rel_a, x = dat$x, na.rm = F)

# Make missing nighttime measurements = 0 in rel.a
dat = dat %>% mutate(rel_a = if_else(x %in% c(0,1,2,22,23),0,rel_a))

# Make the timeseries = 2 days

dat2 = data_frame(x = seq(24,47, by = 1), rel_e = dat$rel_e, rel_a = dat$rel_a)

dat3 = rbind(dat,dat2)

# +/- 1 rolling mean to smooth curves and add time col for 2 days
dat3 = dat3 %>% 
  mutate(roll1_rel_e = slider::slide_dbl(rel_e, mean, .before = 1, .after = 1),
         roll1_rel_a = slider::slide_dbl(rel_a, mean, .before = 1, .after = 1),
         time = seq(ISOdate(2024, 1, 1, 0), by = "hour", length.out = 48))

# Plot

plot.a = dat3 %>% 
  ggplot() +

  geom_line(aes(x = time, y = roll1_rel_e, color = 'Probability for growth'), size = 1) +
  geom_line(aes(x = time, y = roll1_rel_a, color = 'Carbon assimilation'), size = 1) +
  scale_color_manual(values = c('Probability for growth' = 'black', 'Carbon assimilation' = 'green')) +
  
  # Shade under curves
  geom_ribbon(aes(x = time, ymin = 0, ymax = roll1_rel_e), fill = 'black', alpha = 0.1) +
  geom_ribbon(aes(x = time, ymin = 0, ymax = roll1_rel_a), fill = 'green', alpha = 0.1) +
  
  
  # vertical lines at 7am and 7pm. The lines draws at time + 1?
  geom_vline(xintercept = as.POSIXct("2024-01-01 08:00:00"), color = "red", size = .5, alpha = .5) +
  geom_vline(xintercept = as.POSIXct("2024-01-02 08:00:00"), color = "red", size = .5, alpha = .5) +
  geom_vline(xintercept = as.POSIXct("2024-01-01 20:00:00"), color = "blue", size = .5, alpha = .5) +
  geom_vline(xintercept = as.POSIXct("2024-01-02 20:00:00"), color = "blue", size = .5, alpha = .5) +

  xlab('') +
  ylab('Relative value (%)') +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = '4 hours') +
  coord_cartesian(ylim = c(0, 100), expand = FALSE) +
  plotTheme1 +
  theme(legend.position = 'none',
        axis.text.x = element_blank())

plot.a

#### Threshold to induce overwintering processes (Fig. 4b)

dat.chill = data.frame(time = seq(ISOdate(2024, 1, 1, 0), by = "hour", length.out = 48), temp = 
                   rep(c(2.25,2,2,2.25,2.5,2.75,3.25,4.5,5.75,7,8.25,9,9.5,9.75,9.75,9.5,9,8.25,7,5.75,4.5,3.25,2.75,2.5),2))

plot.b = dat.chill %>% ggplot()+
  geom_line(aes(x = time, y = temp), size = 1, colour = 'darkorchid4') +
  geom_ribbon(aes(x = time, ymin = 0, ymax = temp), fill = 'darkorchid4', alpha = 0.1) +
  
  # vertical lines at 7am and 7pm. The lines draws at time + 1?
  geom_vline(xintercept = as.POSIXct("2024-01-01 08:00:00"), color = "red", size = .5, alpha = .5) +
  geom_vline(xintercept = as.POSIXct("2024-01-02 08:00:00"), color = "red", size = .5, alpha = .5) +
  geom_vline(xintercept = as.POSIXct("2024-01-01 20:00:00"), color = "blue", size = .5, alpha = .5) +
  geom_vline(xintercept = as.POSIXct("2024-01-02 20:00:00"), color = "blue", size = .5, alpha = .5) +
  
  xlab('Time') +
  ylab('Temperature (°C)') +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = '4 hours') +
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
  coord_cartesian(ylim = c(0,10), expand = FALSE) +
  plotTheme1 +
  theme(axis.text.x = element_text(margin = margin(t = 6)))

plot.b

# Merge a and b
#######################################

#define plot layout
layout <- "
AAA
BBB"

#Merge plots
CombinedPlot =  plot.a + plot.b +
  plot_layout(design = layout, tag_level = 'new') + plot_annotation(tag_levels = list(c('a','b')))&
  theme(plot.tag = element_text(face = 'bold', size = 22), axis.title = element_text(size = 24), 
        axis.text = element_text(size = 22))

CombinedPlot

#Save PDF
#pdf(paste(output.path,"CombinedPlot.pdf",sep="/"), width=12, height=10, useDingbats=FALSE)
#CombinedPlot 
#dev.off()
