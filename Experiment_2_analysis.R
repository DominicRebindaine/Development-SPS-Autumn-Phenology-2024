###########################################################################################################################################
# R Script for: 
# Experiment 2 of the study:
# Development underpins the summer solstice reversal of the effects of climate warming on the autumn phenology of European beech
# Rebindaine et al. (2024)
# Dominic Rebindaine
# Last edited: 18/09/2024
###########################################################################################################################################

#####################
# Required packages #
#####################

require(ggplot2)
require(tidyverse)
require(broom.mixed)
require(patchwork)
require(lubridate)
require(data.table)
require(tm)
require(lme4)
require(parameters)

##############################################################################################################################################
##############################################################################################################################################



#############################
## Set directory and paths ##
#############################


data.dir = '/Users/Dominic/Documents/Crowther/r/buds/Data'
output.dir = '/Users/Dominic/Documents/Crowther/r/buds/figures'


#####################
# Load bud set data #
#####################


bud.df = read.table(paste(data.dir,"bud_measurements_rm_2022.csv",sep="/"), header=T, sep=",") %>% 
  
  #create ID and treatment columns
  mutate(ID=readr::parse_number(TreeID), #keep only numbers in string
         TreatmentLetter=removeNumbers(TreeID),#remove numbers in string
         date = as.Date(date, format="%d.%m.%Y"),
         #merge controls
         Treatment = ifelse(TreatmentLetter %in% c("A","C", "D","H"),"Control",
                            ifelse(TreatmentLetter %in% c("E"),"Pre_day", 
                                ifelse(TreatmentLetter %in% c("F"),"Pre_night",
                                   ifelse(TreatmentLetter %in% c("G"),"Pre_full", 
                                          ifelse(TreatmentLetter %in% c("I"),"Post_day",
                                                 ifelse(TreatmentLetter %in% c("J"),"Post_night",
                                   ifelse(TreatmentLetter %in% c("K"), "Post_full", "Other")))))))) 
  
  #long format
  bud.df = pivot_longer(bud.df, -c(TreeID, ID, Treatment, TreatmentLetter, week, doy, date), 
               names_to = "bud_type", values_to = "bud_length") %>%
    # remove .length from bud type
  mutate(bud_type = gsub(".length","", bud_type)) %>% 
    # remove B treatment, and NAs
  filter(!Treatment %in% c('Other'),
         !is.na(bud_length)) %>% 
    # group to individ bud
  group_by(ID, bud_type) %>% 
   # only keep buds with more than 8 measurements
  filter(n() >= 8) %>% 
  ungroup()


  ##############################################################################################################################################
  ##############################################################################################################################################
  
# Plot theme
#####################################
  plotTheme1 = theme(
    legend.position   = "none",
    legend.background = element_rect(fill=NA, size=0.5, linetype="solid"),
    legend.text       = element_text(color="black"),
    panel.background  = element_blank(),
    axis.text         = element_text(colour = "black"),
    panel.border      = element_rect(colour = "black", fill=NA),
    axis.line         = element_line(color = "black"),
    strip.background  = element_rect(fill=NA),
    strip.text        = element_text(colour = 'black'),
    plot.title        = element_text(face="bold",hjust = 0.5))


##############################################################################################################################################
######## Sample sizes
##############################################################################################################################################

# Observed buds per treatment
table(bud.df[bud.df$date=="2022-09-01",]$Treatment, bud.df[bud.df$date=="2022-09-01",]$bud_type)

# Number of individual trees
length(unique(bud.df$ID))

# Number of individual trees per treatment
length(unique(bud.df[bud.df$Treatment=="Control",]$ID))

length(unique(bud.df[bud.df$Treatment=="Pre_day",]$ID))

length(unique(bud.df[bud.df$Treatment=="Pre_night",]$ID))

length(unique(bud.df[bud.df$Treatment=="Pre_full",]$ID))

length(unique(bud.df[bud.df$Treatment=="Post_day",]$ID))

length(unique(bud.df[bud.df$Treatment=="Post_night",]$ID))

length(unique(bud.df[bud.df$Treatment=="Post_full",]$ID))


##############################################################################################################################################
######## Data Preparation
##############################################################################################################################################


####################################################
## Get absolute and relative bud growth
####################################################


bud.df = bud.df %>%
  group_by(ID, bud_type)%>%
  mutate(
    #get relative bud length per individual
    bud.rel = bud_length / max(bud_length),
    #get bud growth rates per individual
    across(bud_length, ~c(NA, diff(.)), .names = "bud_growth"),
    #get relative bud growth rates per individual
    across(bud.rel, ~c(NA, diff(.)), .names = "relbud_growth"),
    #get date of maximum bud length per individual
    date.max.bud = date[which.max(bud_length)]) %>%
  ungroup()


######################################################################
## Mean abs and rel bud growth per treatment
######################################################################

minbud.df = bud.df %>% 
  group_by(ID, bud_type) %>% 
  filter(row_number() <= which.min(bud_length)) %>% 
  select(ID, Treatment, bud_type, bud_length, bud.rel, date.max.bud)

maxbud.df = bud.df %>% 
  group_by(ID, bud_type) %>% 
  filter(row_number() == min(which(bud.rel > 0.9))) %>% 
  select(ID, Treatment, bud_type, bud_length, bud.rel, date.max.bud)

minmaxbud.df = merge(minbud.df, maxbud.df, by = c('ID', 'bud_type'))

minmaxbud.df = minmaxbud.df %>% 
  group_by(ID, bud_type) %>% 
  mutate(diff = bud_length.y - bud_length.x,
         total.rel.growth = (1-bud_length.x/bud_length.y)*100)

averages = minmaxbud.df %>% 
  group_by(Treatment.x) %>% 
  summarise(meandiff = round(mean(diff),2),
            mean.total.rel = round(mean(total.rel.growth),2))


######################################################################
## Get buds that are >= 90% their total length at first measurement
######################################################################

over90start = bud.df %>% filter(doy == 237) %>% filter(bud.rel >= 0.9)


######################################################################
## Interpolate 90% bud set dates from seasonal bud length measurements
######################################################################


bud.data.90 = bud.df %>%
  group_by(ID, Treatment, bud_type) %>%
  # remove buds that were >= 90% relative bud length on day 1
  filter(!(ID %in% c(116) & bud_type == 'apical'),
         !(ID %in% c(148) & bud_type == 'apical'),
         !(ID %in% c(149) & bud_type == 'apical'),
         !(ID %in% c(156) & bud_type == 'apical'),
         !(ID %in% c(158) & bud_type == 'apical'),
         !(ID %in% c(58) & bud_type == 'lateral'),
         !(ID %in% c(80) & bud_type == 'lateral'),
         !(ID %in% c(82) & bud_type == 'lateral'),
         !(ID %in% c(99) & bud_type == 'lateral'),
         !(ID %in% c(147) & bud_type == 'lateral'),
         !(ID %in% c(153) & bud_type == 'lateral'),
         !(ID %in% c(154) & bud_type == 'lateral'),
         !(ID %in% c(156) & bud_type == 'lateral'),
         !(ID %in% c(157) & bud_type == 'lateral'),
         !(ID %in% c(158) & bud_type == 'lateral'),
         !(ID %in% c(159) & bud_type == 'lateral')) %>% 
  # Find the last row with <90% bud length in each group, keep only this row and next
   filter(row_number() <= min(which(bud.rel > 0.9)),
          row_number() >= min(which(bud.rel > 0.9))-1) %>%
  #Linear interpolation
  summarize(EOS     = as.Date(approx(bud.rel, date, .9,  ties=min)$y, origin="1970-01-01")) %>%
  mutate(EOS.DOY    = yday(EOS)) %>%
  ungroup()

# Observed buds per treatment after removal
table(bud.data.90$Treatment, bud.data.90$bud_type)


# Buds analysis data summary
##############################
BudsSummary = merge(bud.data.90, minmaxbud.df, by = c('ID', 'bud_type')) %>% 
  select(ID, bud_type, Treatment, EOS, EOS.DOY, bud_length.x, bud_length.y, diff, total.rel.growth)

AverageSummary = BudsSummary %>% 
  group_by(Treatment) %>% 
  summarise(
    meanEOS.DOY = round(mean(EOS.DOY),1),
            ci.lower.EOS = round(meanEOS.DOY - 1.96*standard_error(EOS.DOY),1),
            ci.upper.EOS = round(meanEOS.DOY + 1.96*standard_error(EOS.DOY),1),
    meandiff = round(mean(diff),2),
            ci.lower.diff = round(meandiff - 1.96*standard_error(diff),2),
            ci.upper.diff = round(meandiff + 1.96*standard_error(diff),2),
    mean.total.rel = round(mean(total.rel.growth),2),
            ci.lower.rel = round(mean.total.rel - 1.96*standard_error(total.rel.growth),2),
            ci.upper.rel = round(mean.total.rel + 1.96*standard_error(total.rel.growth),2))



##############################################################################################################################################
######## Bud set analysis
##############################################################################################################################################


############################
## Mixed effects models
############################


# Bud set dates
###############

resultsLM = bud.data.90 %>% 
  do({model = lmer(EOS.DOY ~ Treatment + (1|bud_type), data=.)  # create the model
  data.frame(tidy(model, effect="fixed"))}) %>%
  filter(!term %in% c("(Intercept)")) %>% 
  mutate(term = factor(term, levels = c("TreatmentPre_day", "TreatmentPre_night", "TreatmentPre_full",
                                        "TreatmentPost_day", "TreatmentPost_night", "TreatmentPost_full")),
         SEx2 = std.error*2)
       
#Show results
as.data.frame(resultsLM) %>%
  dplyr::select(-c(statistic))%>%
  mutate_if(is.numeric, round, digits=2)


# Autumn bud growth rates
#########################

resultsRelGrowth = bud.df %>% 
  do({model = lmer(relbud_growth*100/7 ~ Treatment + (1|bud_type) + (1|date), data=.)  # create your model
  data.frame(tidy(model, effect="fixed"))}) %>% 
  filter(!term %in% c("(Intercept)")) %>% 
  mutate(term = factor(term, levels = c("TreatmentPre_day", "TreatmentPre_night", "TreatmentPre_full",
                                        "TreatmentPost_day", "TreatmentPost_night", "TreatmentPost_full")))
      
#Show results
as.data.frame(resultsRelGrowth) %>%
  dplyr::select(-c(statistic))%>%
  mutate_if(is.numeric, round, digits=2)

# Relative Autumn bud growth
#############################

BudGrowth.df = bud.df %>% 
  group_by(Treatment, ID, bud_type) %>% 
  summarize(RelAutGrowth = (1-(min(bud_length)/max(bud_length)))*100) %>% 
  ungroup()

resultsRelAutGrowth = BudGrowth.df %>% 
  do({model = lmer(RelAutGrowth ~ Treatment + (1|bud_type), data=.)
  data.frame(tidy(model, effect="fixed"))}) %>% 
filter(!term %in% c("(Intercept)")) %>% 
  mutate(term = factor(term, levels = c("TreatmentPre_day", "TreatmentPre_night", "TreatmentPre_full",
                                        "TreatmentPost_day", "TreatmentPost_night", "TreatmentPost_full")),
         SEx2 = std.error*2)


#Show results
Intercept = resultsRelAutGrowth[1,]$estimate
as.data.frame(resultsRelAutGrowth) %>%
  #get relative autumn growth change
  mutate(PercentGrowthChange = estimate/Intercept*100,
         PercentGrowth.se = std.error/Intercept*100) %>% 
  dplyr::select(-c(statistic))%>%
  mutate_if(is.numeric, round, digits=1)

# Absolute Autumn bud growth
#############################

resultsAbsAutGrowth = BudsSummary %>% 
  do({model = lmer(diff ~ Treatment + (1|bud_type), data=.)
  data.frame(tidy(model, effect="fixed"))}) %>% 
  filter(!term %in% c("(Intercept)")) %>% 
  mutate(term = factor(term, levels = c("TreatmentPre_day", "TreatmentPre_night", "TreatmentPre_full",
                                        "TreatmentPost_day", "TreatmentPost_night", "TreatmentPost_full")),
         SEx2 = std.error*2)


##############################################################################################################################################
######## Figures
##############################################################################################################################################


# Bud set dates
###############

## all

LMplot90 = ggplot() + 
  scale_color_manual(values = c(rep('#F21A00',3), rep('#3B9AB2',3)))+
  geom_hline(yintercept=0)+
  geom_point(data=resultsLM, aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=resultsLM, 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=3.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-10,10))+
  xlab("")+ylab("Bud set anomaly (days)")+
  scale_x_discrete(labels = c('Pre_day','Pre_night', 'Pre_full','Post_day', 'Post_night', 'Post_full'))+
  plotTheme1

LMplot90

## full cooling

LMplot90.full = ggplot() + 
  scale_color_manual(values = c('#F21A00', '#3B9AB2'))+
  geom_hline(yintercept=0)+
  geom_point(data= filter(resultsLM, term %in% c('TreatmentPre_full', 'TreatmentPost_full')),
             aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=filter(resultsLM, term %in% c('TreatmentPre_full', 'TreatmentPost_full')), 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=1.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-10,10))+
  xlab("")+ylab("Bud set anomaly (days)")+
  scale_x_discrete(labels = c('Pre','Post'))+
  plotTheme1 +
  ggtitle('Full cooling')

LMplot90.full

## Day cooling

LMplot90.day = ggplot() + 
  scale_color_manual(values = c('#F21A00', '#3B9AB2'))+
  geom_hline(yintercept=0)+
  geom_point(data= filter(resultsLM, term %in% c('TreatmentPre_day', 'TreatmentPost_day')), aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=filter(resultsLM, term %in% c('TreatmentPre_day', 'TreatmentPost_day')), 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=1.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-10,10))+
  xlab("")+ylab("Bud set anomaly (days)")+
  scale_x_discrete(labels = c('Pre','Post'))+
  plotTheme1 +
  ggtitle('Day cooling') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())

LMplot90.day

## Night cooling

LMplot90.night = ggplot() + 
  scale_color_manual(values = c('#F21A00', '#3B9AB2'))+
  geom_hline(yintercept=0)+
  geom_point(data= filter(resultsLM, term %in% c('TreatmentPre_night', 'TreatmentPost_night')), aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=filter(resultsLM, term %in% c('TreatmentPre_night', 'TreatmentPost_night')), 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=1.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-10,10))+
  xlab("")+ylab("Bud set anomaly (days)")+
  scale_x_discrete(labels = c('Pre','Post'))+
  plotTheme1 +
  ggtitle('Night cooling') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())

LMplot90.night

## Pre-solstice cooling

LMplot90.pre = ggplot() + 
  scale_color_manual(values = c('#F21A00', '#F21A00','#F21A00'))+
  geom_hline(yintercept=0)+
  geom_point(data= filter(resultsLM, term %in% c('TreatmentPre_day', 'TreatmentPre_night','TreatmentPre_full')), 
             aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=filter(resultsLM, term %in% c('TreatmentPre_day', 'TreatmentPre_night','TreatmentPre_full')), 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=1.5, size=2, alpha=0.4)+
  geom_vline(xintercept=2.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-10,10))+
  xlab("")+ylab("Bud set anomaly (days)")+
  scale_x_discrete(labels = c('Day','Night', 'Full'))+
  plotTheme1 +
  ggtitle('Pre-solstice cooling')

LMplot90.pre

## Post-solstice cooling

LMplot90.post = ggplot() + 
  scale_color_manual(values = c('#3B9AB2', '#3B9AB2','#3B9AB2'))+
  geom_hline(yintercept=0)+
  geom_point(data= filter(resultsLM, term %in% c('TreatmentPost_day', 'TreatmentPost_night','TreatmentPost_full')), 
             aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=filter(resultsLM, term %in% c('TreatmentPost_day', 'TreatmentPost_night','TreatmentPost_full')), 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=1.5, size=2, alpha=0.4)+
  geom_vline(xintercept=2.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-10,10))+
  xlab("")+ylab("Bud set anomaly (days)")+
  scale_x_discrete(labels = c('Day','Night', 'Full'))+
  plotTheme1 +
  ggtitle('Post-solstice cooling') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())
LMplot90.post

# Relative bud growth anomaly
#############################

## all

LMplotRelAutGrowth = ggplot() + 
  scale_color_manual(values = c(rep('#F21A00',3), rep('#3B9AB2',3)))+
  geom_hline(yintercept=0)+
  geom_point(data=resultsRelAutGrowth, aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=resultsRelAutGrowth, aes(x=term, ymin=estimate-1.96*std.error, ymax=estimate+1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=3.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-15,10))+
  xlab("")+ylab("Relative bud growth anomaly (%)")+
  scale_x_discrete(label =c('Pre_day','Pre_night', 'Pre_full','Post_day', 'Post_night', 'Post_full'))+
  plotTheme1

LMplotRelAutGrowth 

#ggsave(filename = "RelBudGrowthExp2.pdf", path = output.dir, plot = LMplotRelAutGrowth)

# Absolute bud growth anomaly
#############################

## all

LMplotAbsAutGrowth = ggplot() + 
  scale_color_manual(values = c(rep('#F21A00',3), rep('#3B9AB2',3)))+
  geom_hline(yintercept=0)+
  geom_point(data=resultsAbsAutGrowth, aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=resultsAbsAutGrowth, aes(x=term, ymin=estimate-1.96*std.error, ymax=estimate+1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=3.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-5,5))+
  xlab("")+ylab("Absolute bud growth anomaly (mm)")+
  scale_x_discrete(label =c('Pre_day','Pre_night', 'Pre_full','Post_day', 'Post_night', 'Post_full'))+
  plotTheme1

LMplotAbsAutGrowth 

# Relative bud growth rates
#############################

## all

LMplotRelGrowth = ggplot() + 
  scale_color_manual(values = c(rep('#F21A00',3), rep('#3B9AB2',3)))+
  geom_hline(yintercept=0)+
  geom_point(data=resultsRelGrowth, aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=resultsRelGrowth, 
                aes(x=term, ymin=estimate-1.96*std.error, ymax=estimate+1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=3.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-0.2,.2))+
  xlab("")+ylab("Bud growth rate anomaly (% per day)")+
  scale_x_discrete(c('Pre_day','Pre_night', 'Pre_full','Post_day', 'Post_night', 'Post_full'))+
  plotTheme1

LMplotRelGrowth 

## full

LMplotRelGrowth.full = ggplot() + 
  scale_color_manual(values = c('#F21A00', '#3B9AB2'))+
  geom_hline(yintercept=0)+
  geom_point(data= filter(resultsRelGrowth, term %in% c('TreatmentPre_full', 'TreatmentPost_full')),
             aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=filter(resultsRelGrowth, term %in% c('TreatmentPre_full', 'TreatmentPost_full')), 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=1.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-0.2,0.2))+
  xlab("")+ylab("Bud growth rate anomaly (% per day)")+
  scale_x_discrete(labels = c('Pre','Post'))+
  plotTheme1 +
  ggtitle('Full cooling')

LMplotRelGrowth.full

## day

LMplotRelGrowth.day = ggplot() + 
  scale_color_manual(values = c('#F21A00', '#3B9AB2'))+
  geom_hline(yintercept=0)+
  geom_point(data= filter(resultsRelGrowth, term %in% c('TreatmentPre_day', 'TreatmentPost_day')),
             aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=filter(resultsRelGrowth, term %in% c('TreatmentPre_day', 'TreatmentPost_day')), 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=1.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-0.2,0.2))+
  xlab("")+ylab("Bud growth rate anomaly (% per day)")+
  scale_x_discrete(labels = c('Pre','Post'))+
  plotTheme1 +
  ggtitle('Day cooling') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())
  

LMplotRelGrowth.day


## Night

LMplotRelGrowth.night = ggplot() + 
  scale_color_manual(values = c('#F21A00', '#3B9AB2'))+
  geom_hline(yintercept=0)+
  geom_point(data= filter(resultsRelGrowth, term %in% c('TreatmentPre_night', 'TreatmentPost_night')),
             aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=filter(resultsRelGrowth, term %in% c('TreatmentPre_night', 'TreatmentPost_night')), 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  
  geom_vline(xintercept=1.5, size=2, alpha=0.4)+
  coord_cartesian(ylim=c(-0.2,0.2))+
  xlab("")+ylab("Bud growth rate anomaly (% per night)")+
  scale_x_discrete(labels = c('Pre','Post'))+
  plotTheme1 +
  ggtitle('Night cooling') +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())


LMplotRelGrowth.night


# Bud length curves
###################

budpal = c('black', 'magenta3', 'red', 'brown3', 'steelblue3', 'blue4', 'cyan2')


#move points .05 to the left and right
pd = position_dodge(2) 


# Relative bud length
RelBudLength_curve = bud.df %>%
  mutate(Treatment = factor(Treatment, 
  levels=c("Control",'Pre_day','Pre_night', 'Pre_full','Post_day', 'Post_night', 'Post_full'))) %>% 
  ggplot(aes(x=date, y=bud.rel*100, colour=Treatment, group=Treatment)) +
  
  geom_hline(yintercept = 60, colour="lightgrey")+
  geom_hline(yintercept = 70, colour="lightgrey")+
  geom_hline(yintercept = 80, colour="lightgrey")+
  geom_hline(yintercept = 90, colour="lightgrey")+
  geom_hline(yintercept = 100, colour="lightgrey")+
  
  geom_vline(xintercept = as.Date('2022-11-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-10-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-09-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-08-01'), colour="lightgrey")+
  
  stat_summary(fun.data = "mean_se", geom="line", size = 1.2, position=pd, alpha=0.7) +
  stat_summary(fun.data = "mean_se", geom="errorbar", size = 0.8, width=0, position=pd) +
  stat_summary(fun.data = "mean_se", geom="point", size = 1.2, position=pd) +
  
  scale_color_manual(values = budpal)+
  
  coord_cartesian(ylim=c(65,100),xlim=c(as.Date(c('2022-08-25','2022-11-05'))))+
  
  xlab("") +
  ylab("Relative bud length (%)") +
  
  #ggtitle("Relative bud growth")+
  
  plotTheme1 +
  theme(legend.position       = c(0.8, 0.5),
        legend.key            = element_rect(fill = "transparent"),
        legend.background     = element_rect(fill='white'),
        legend.box.background = element_rect(fill='transparent'),
        axis.title.x          = element_blank(),
        axis.text.x           = element_blank())

RelBudLength_curve


# Absolute bud growth
#####################


AbsBudLength_curve = bud.df %>%
  mutate(Treatment = factor(Treatment, 
    levels=c("Control",'Pre_day','Pre_night', 'Pre_full','Post_day', 'Post_night', 'Post_full'))) %>% 
  ggplot(aes(x=date, y=bud_length, colour=Treatment, group=Treatment)) +
  
  geom_hline(yintercept = 12, colour="lightgrey")+
  geom_hline(yintercept = 16, colour="lightgrey")+
  #geom_hline(yintercept = 17.5, colour="lightgrey")+
  geom_hline(yintercept = 20, colour="lightgrey")+
  geom_hline(yintercept = 24, colour="lightgrey")+
  #geom_hline(yintercept = 26, colour="lightgrey")+
  
  geom_vline(xintercept = as.Date('2022-11-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-10-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-09-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-08-01'), colour="lightgrey")+
  
  stat_summary(fun.data = "mean_se", geom="line", size = 1.2, position=pd, alpha=0.7) +
  stat_summary(fun.data = "mean_se", geom="errorbar", size = 0.8, width=0, position=pd) +
  stat_summary(fun.data = "mean_se", geom="point", size = 1.2, position=pd) +
  
  scale_color_manual(values = budpal)+
  
  coord_cartesian(ylim=c(12.5,26),xlim=c(as.Date(c('2022-08-25','2022-11-05'))))+
  
  xlab("") +
  ylab("Bud length (mm)") +
  
  #ggtitle("Absolute bud growth") +
  
  plotTheme1 

AbsBudLength_curve


# Relative bud growth rate
##########################


RelBudGrowth_curve = bud.df %>%
  ggplot(aes(x=date, y=relbud_growth*100/7, colour=Treatment, group=Treatment)) +
  
  geom_hline(yintercept = 0.25, colour="lightgrey")+
  geom_hline(yintercept = 0.5, colour="lightgrey")+
  geom_hline(yintercept = 0.75, colour="lightgrey")+
  geom_hline(yintercept = 1, colour="lightgrey")+
  
  geom_vline(xintercept = as.Date('2022-11-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-10-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-09-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-08-01'), colour="lightgrey")+
  
  stat_summary(fun.data = "mean_se", geom="line", size = 1.2, position=pd, alpha=0.7) +
  stat_summary(fun.data = "mean_se", geom="errorbar", size = 0.8, width=0, position=pd) +
  stat_summary(fun.data = "mean_se", geom="point", size = 1.2, position=pd) +
  
  scale_color_manual(values = budpal)+
  
  coord_cartesian(ylim=c(0.056,1.2),xlim=c(as.Date(c('2022-08-25','2022-11-05'))))+
  
  xlab("") +
  ylab("Relative growth per day (%)") +
  
  #ggtitle("Bud growth rate")+
  
  plotTheme1 
 

RelBudGrowth_curve


# Absolute bud growth
#####################


AbsBudGrowth_curve = bud.df %>%
  ggplot(aes(x=date, y=bud_growth/7, colour=Treatment, group=Treatment)) +
  
  geom_hline(yintercept = .1, colour="lightgrey")+
  geom_hline(yintercept = .2, colour="lightgrey")+
  geom_hline(yintercept = .3, colour="lightgrey")+
  
  geom_vline(xintercept = as.Date('2022-11-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-10-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-09-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2022-08-01'), colour="lightgrey")+
  
  stat_summary(fun.data = "mean_se", geom="line", size = 1.2, position=pd, alpha=0.7) +
  stat_summary(fun.data = "mean_se", geom="errorbar", size = 0.8, width=0, position=pd) +
  stat_summary(fun.data = "mean_se", geom="point", size = 1.2, position=pd) +
  
  scale_color_manual(values = budpal)+
  
  coord_cartesian(ylim=c(0.011,0.25),xlim=c(as.Date(c('2022-08-25','2022-11-05'))))+
  
  xlab("") +
  ylab("mm per day") +
  
  #ggtitle("Absolute bud growth")+
  
  plotTheme1 

AbsBudGrowth_curve

## Bud set anomaly Pre only full-night-day order
LMplot90.AC = ggplot() + 
  scale_color_manual(values = c(rep('#F21A00',3)))+
  geom_hline(yintercept=0)+
  geom_point(data=resultsLM %>% filter(term %in% c('TreatmentPre_day','TreatmentPre_night','TreatmentPre_full')), aes(x=term, y=estimate, color=term))+
  geom_errorbar(data=resultsLM%>% filter(term %in% c('TreatmentPre_day','TreatmentPre_night','TreatmentPre_full')), 
                aes(x=term, ymin=estimate+1.96*std.error, ymax=estimate-1.96*std.error, color=term), 
                width=.2, position=position_dodge(.9)) +
  coord_cartesian(ylim=c(-10,10))+
  xlab("")+ylab("Bud set anomaly (days)")+
  scale_x_discrete(labels = c('Day','Night', 'Full'))+
  plotTheme1 +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13))

LMplot90.AC

#ggsave(filename = 'Pre_only.png', path = '/Users/dominic/Documents/Crowther/admin/PhD/AC', plot = LMplot90.AC, height = 5, width = 7)

# Merge plots
#####################

#define plot layout
layout <- "
ACC
BDD"

#Merge plots
CombinedPlot =  LMplotRelGrowth + LMplot90+ RelBudLength_curve + AbsBudLength_curve +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('A')) &
  theme(plot.tag = element_text(face = 'bold'))

#Save PDF
#pdf(paste(output.dir,"CombinedPlot_Experiment2.pdf",sep="/"), width=7, height=7, useDingbats=FALSE)
#CombinedPlot 
#dev.off()


# Merge separated bud set anomaly plots
#######################################

#define plot layout
layout2 <- "
ABC"

#Merge plots
CombinedBudSepPlot =  LMplot90.full + LMplot90.day + LMplot90.night +
  plot_layout(design = layout2, tag_level = 'new') + plot_annotation(tag_levels = list(c('','','')))&
  theme(plot.tag = element_text(face = 'bold'), axis.title = element_text(size = 22), 
        axis.text = element_text(size = 18))

#Save PDF
#pdf(paste(output.dir,"CombinedBudSepPlot.pdf",sep="/"), width=7, height=7, useDingbats=FALSE)
#CombinedBudSepPlot 
#dev.off()

# Merge pre- and post plots
############################

# #define plot layout
layout3 <- "
AB"

#Merge plots
CombinedBudSepPlotPre_Post =  LMplot90.pre + LMplot90.post +
  plot_layout(design = layout3, tag_level = 'new') + plot_annotation(tag_levels = list(c('','','')))&
  theme(plot.tag = element_text(face = 'bold'), axis.title = element_text(size = 22), 
        axis.text = element_text(size = 18))

#Save PDF
#pdf(paste(output.dir,"CombinedBudSepPlotPre_Post.pdf",sep="/"), width=7, height=7, useDingbats=FALSE)
#CombinedBudSepPlotPre_Post 
#dev.off()

# Merge separated bud growth anomaly plots
##########################################

#Merge plots
CombinedRelBudGroSepPlot =  LMplotRelGrowth.full + LMplotRelGrowth.day + LMplotRelGrowth.night +
  plot_layout(design = layout2, tag_level = 'new') + plot_annotation(tag_levels = list(c('D','',''))) &
  theme(plot.tag = element_text(face = 'bold'))

#Save PDF
#pdf(paste(output.dir,"CombinedRelBudGroSepPlot.pdf",sep="/"), width=7, height=7, useDingbats=FALSE)
#CombinedRelBudGroSepPlot 
#dev.off()


# Merge bud length curves
#################################################
layout4 = "
A
B"

#Merge plots
CombinedBudLengthCurves = RelBudGrowth_curve + AbsBudGrowth_curve +
  plot_layout(design = layout4, tag_level = 'new') +  plot_annotation(tag_levels = list(c('E', 'F'))) &
  theme(plot.tag = element_text(face = 'bold'))


 
#Save PDF
#pdf(paste(output.dir,"CombinedBudLengthCurves.pdf",sep="/"), width=7, height=7, useDingbats=FALSE)
#CombinedBudLengthCurves 
#dev.off()


#######################
# Individual Bud length
#######################


budlength_individuals = bud.df %>%
  ggplot(aes(x=date, y=bud_length, 
             colour=Treatment, fill=Treatment)) +
  
  geom_vline(xintercept = as.Date('2021-11-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2021-10-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2021-09-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2021-08-01'), colour="lightgrey")+
  geom_vline(xintercept = as.Date('2021-07-01'), colour="lightgrey")+
  
  geom_hline(yintercept = 10, colour="lightgrey")+
  geom_hline(yintercept = 20, colour="lightgrey")+
  geom_hline(yintercept = 30, colour="lightgrey")+
  
  geom_line(aes(x=date, y=bud_length), color="black")+
  geom_point()+
  geom_area(aes(x=date, y=bud_length), alpha=0.5) +
  
  geom_vline(data=bud.data.90, 
             aes(xintercept = EOS), linetype = "dashed")+
  
  xlab("Date") +
  ylab("Bud length (mm)") +
  
  facet_grid(ID~bud_type)+
  plotTheme1

#Save PDF
#pdf(paste(output.dir,"IndividualPlot_Budlength.pdf",sep="/"), width=6, height=120, useDingbats=FALSE)
#budlength_individuals
#dev.off()


