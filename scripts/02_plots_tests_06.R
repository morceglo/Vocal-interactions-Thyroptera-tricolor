# Does kinship predict vocal interaction rates in Thyroptera tricolor?
# Gerald Carter, gcarter1640@gmail.com

# clear workspace
rm(list=ls())

# set ggplot theme
theme_set(theme_bw(base_size = 14))

# set directory to get data
setwd("/Users/gerry/Dropbox/Dropbox/_working/_ACTIVE/_collaborators/Gloriana/thyroptera-vocal-interaction/code")

### OR manually choose this file: 'calling-association-kinship-data.csv'
setwd(dirname(file.choose()))

# get data
raw <- read.csv('calling-association-kinship-data.csv') 

# set directory to store results
setwd("/Users/gerry/Dropbox/Dropbox/_working/_ACTIVE/_collaborators/Gloriana/thyroptera-vocal-interaction/code/results-2023-09-14")

# load packages
library(tidyverse)
library(lme4)
library(glmmTMB)
library(broom.mixed)
library(boot)
library(patchwork)
library(performance)
library(MASS)

###### Create functions------------------

# function to get mean and 95% CI of values x via bootstrapping
boot_ci <- function(x, perms=5000, bca=F) {
  get_mean <- function(x, d) {
    return(mean(x[d]))
  }
  x <- as.vector(na.omit(x))
  mean <- mean(x)
  if(bca){
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="bca")
    low <- boot$bca[1,4]
    high <- boot$bca[1,5] 
  }else{
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="perc")
    low <- boot$perc[1,4]
    high <- boot$perc[1,5] 
  }
  c(low=low,mean=mean,high=high, N=round(length(x)))
}


# function to get mean and 95% CI via bootstrapping of values y within grouping variable x
boot_ci2 <- function(d=d, y=d$y, x=d$x, perms=5000, bca=F){
  df <- data.frame(effect=unique(x))
  df$low <- NA
  df$mean <- NA
  df$high <- NA
  df$n.obs <- NA
  for (i in 1:nrow(df)) {
    ys <- y[which(x==df$effect[i])]
    if (length(ys)>1 & var(ys)>0 ){
      b <- boot_ci(y[which(x==df$effect[i])], perms=perms, bca=bca) 
      df$low[i] <- b[1]
      df$mean[i] <- b[2]
      df$high[i] <- b[3]
      df$n.obs[i] <- b[4]
    }else{
      df$low[i] <- min(ys)
      df$mean[i] <- mean(ys)
      df$high[i] <- max(ys)
      df$n.obs[i] <- length(ys)
    }
  }
  df
}

# look at data
head(raw)
# group is the year and social group
# dyad is the undirected pair (flying and roosting bat) in alphanumeric order
# sri is simple ratio index of association
# flying.sex is sex of flying bat
# roosting.sex is sex of roosting bat
# dyad.sex is the sex of the flying and roosting bat
# udyad.sex is females, males, or mixed

# add a few variables to the data
d <- 
  raw %>% 
  # get log count of inquiry calls
  mutate(log_inquiry = log(n_inquiry)) %>% 
  # rescale kinship so 1 unit is 0.5 
  mutate(kinship2 = kinship*2) %>% 
  # if the roosting bat leaves the roost (escapes) then flight time is < 300 s
  mutate(escape= as.numeric(flight_time<300)) %>% 
  # convert flight time from seconds to minutes (for easier interpretation)
  mutate(flight_time2 = flight_time/60)


###### Describe the data-------------

# how many trials?
d %>% nrow()
#446

# how many undirected pairs have vocal data? 
d %>% 
  pull(dyad) %>% 
  unique() %>% 
  length
#139 undirected pairs

d %>% 
  group_by(bat_flying, bat_roosting) %>% 
  summarize(n=n()) 
# 254 directed pairs

# how many groups? 
n_distinct(d$group) 
# 23

# how many of these have association data
d %>% 
  group_by(dyad) %>% 
  summarize(sri= mean(sri)) %>% 
  filter(!is.na(sri))
# 133 pairs have association

# how many of these have kinship data
d %>% 
  group_by(dyad) %>% 
  summarize(kinship= mean(kinship)) %>% 
  filter(!is.na(kinship))
# 82 pairs have kinship data

# how many related undirected pairs
d %>% 
  group_by(dyad) %>% 
  summarize(kinship= mean(kinship)) %>% 
  filter(kinship>0)
# 32 kin pairs

# how many unrelated undirected pairs
d %>% 
  group_by(dyad) %>% 
  summarize(kinship= mean(kinship)) %>% 
  filter(kinship==0)
# 50 nonkin pairs

# what is mean kinship in group?
(grp.kin <- 
  d %>% 
  group_by(dyad, group) %>% 
  summarize(kinship= mean(kinship, na.rm=T)) %>% 
  group_by(group) %>% 
  summarize(kinship= mean(kinship, na.rm=T), n=n()) %>% 
  as.data.frame())

set.seed(123)
grp.kin %>% pull(kinship) %>% boot_ci(bca=T)
# 0.23, 95% CI = [0.16, 0.31]

# how many groups have kin?
grp.kin %>% filter(kinship>0)
# 17

# how many groups have no kin?
grp.kin %>% filter(kinship==0)
# 4
  
# plot distribution of kinship
d %>% 
  ggplot(aes(x=kinship))+
  facet_wrap(~udyad.sex, ncol=1, scales= 'free_y')+
  geom_histogram(fill="light blue", color="black")+
  ggtitle("pairwise kinship")

# check categories
d %>%
  filter(!is.na(kinship)) %>% 
  group_by(dyad) %>% 
  summarize(kinship= mean(kinship)) %>% 
  group_by(kinship) %>% 
  summarize(n=n()) 

##### Effect of kinship on association-------------

# plot distribution of assoc
d %>% 
  ggplot(aes(x=sri))+
  facet_wrap(~udyad.sex, ncol=1, scale ="free_y")+
  geom_histogram(fill="light blue", color="black")+
  xlab("association rate")+
  ggtitle("pairwise SRI values")
# looks ok

# what is mean and 95% CI of assoc among flying and roosting bats?
set.seed(123)
d %>% 
  group_by(dyad) %>% 
  summarize(sri = mean(sri)) %>% 
  pull(sri) %>% 
  boot_ci(bca=T)
# 0.51, 95% CI = [0.47, 0.55]
# a bit low, expected from past work is 0.76

# now only nonkin
set.seed(123)
k1 <- 
  d %>% 
  filter(kinship==0) %>% 
  group_by(dyad) %>% 
  summarize(sri = mean(sri)) %>% 
  pull(sri) %>% 
  boot_ci(bca=T) %>% 
  c(kinship= 0)
k1
# 0.572, 95% CI = [0.500, 0.636]

# now only kin
set.seed(123)
d %>% 
  filter(kinship>0) %>% 
  group_by(dyad) %>% 
  summarize(sri = mean(sri)) %>% 
  pull(sri) %>% 
  boot_ci(bca=T)
# 0.686, 95% CI = [0.631, 0.741]

# now only close kin
set.seed(123)
k2 <- 
  d %>% 
  filter(kinship==0.5) %>% 
  group_by(dyad) %>% 
  summarize(sri = mean(sri)) %>% 
  pull(sri) %>% 
  boot_ci(bca=T)%>% 
  c(kinship= 0.5)
k2
# 0.687, 95% CI = [0.628, 0.750]

# now only 0.25 kin
set.seed(123)
k3 <- 
  d %>% 
  filter(kinship==0.25) %>% 
  group_by(dyad) %>% 
  summarize(sri = mean(sri)) %>% 
  pull(sri) %>% 
  boot_ci(bca=T)%>% 
  c(kinship= 0.25)
k3
# 0.684, 95% CI = [0.542, 0.783]

# compile means
k <- 
  rbind(k1,k2,k3) %>% 
  data.frame() %>% 
  mutate(kin= as.character(kinship)) %>% 
  mutate(kin2= kinship>0) %>% 
  mutate(assoc= mean)

# plot association by kinship
(kin.dyads.plot<- 
    d %>% 
    filter(!is.na(kinship)) %>% 
    group_by(dyad) %>% 
    summarize(kinship= mean(kinship), assoc= mean(sri), udyad.sex= first(udyad.sex)) %>% 
    mutate(kin2= kinship>0) %>%  
    mutate(kin= as.character(kinship)) %>% 
    ggplot(aes(x=kin, y=assoc, group= kin, color=kin2))+
    geom_jitter(size=2, width= 0.1, height=0, aes(shape= udyad.sex))+
    geom_point(data= k, position = position_nudge(x = 0.25), size=3)+
    geom_errorbar(data=k, aes(ymin=low, ymax=high, width=.1), position = position_nudge(x = 0.25), size=1)+
    xlab("kinship")+
    ylab("association rate (simple ratio index)")+
    scale_color_manual(values= c("darkgrey", "darkred"))+
    theme(legend.position= "none", legend.title = element_blank()))

# save as PDF
ggsave(
  "kin_association.pdf",
  plot = kin.dyads.plot,
  scale = 1,
  width = 3,
  height = 5,
  units = "in",
  dpi = 300)


##### Effect of flight time on calling ----------

### How does flight time and count of inquiry calls predict count of response calls?

# plot number of response calls by flight time
t1 <- 
  d %>% 
  mutate(kinship = kinship>0.1) %>% 
  ggplot(aes(x=flight_time, y=n_response))+
  geom_point(size=2, alpha=0.3, aes(color= kinship, shape= kinship))+ 
  geom_smooth(method="glm.nb")+
  scale_color_manual(values= c("darkgrey", "darkred"))+
  xlab("flight time (seconds)")+
  ylab("response call count")+
  theme(legend.position= "none")

# plot number of inquiry calls by flight time
t2 <- 
  d %>% 
  mutate(kinship = kinship>0.1) %>% 
  ggplot(aes(x=flight_time, y=n_inquiry))+
  geom_point(size=2, alpha=0.3, aes(color= kinship, shape=kinship))+ 
  geom_smooth(method = "glm.nb")+
  scale_color_manual(values= c("darkgrey", "darkred"))+
  xlab("flight time (seconds)")+
  ylab("inquiry call count")+
  theme(legend.position= "none")

# plot together
t <- t1+t2
ggsave(
  filename= "flight_time.pdf",
  plot = t,
  scale = 1,
  width = 7,
  height = 7,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# try poisson model for effect of flight time on inquiry calling
t <- glmer (n_inquiry ~ flight_time2 + (1|bat_flying) + (1|bat_roosting) + (1|group), data= d, family= poisson)
check_overdispersion(t)

# negative binomial model (NBM)
t <- glmmTMB(n_inquiry ~ flight_time2 + (1|bat_flying) + (1|bat_roosting) + (1|group),
                data=d,
                ziformula=~0,
                family=nbinom2)
check_overdispersion(t)
t2 <- tidy(t,conf.int=TRUE,exponentiate=T,effects="fixed", conf.method="profile")
t2
# for every minute of flight time, the inquiry call count increases by a factor of 1.43 [1.36, 1.49]

# try poisson model for effect of flight time on response calling
t <- glmer (n_response ~ flight_time2 + log_inquiry + (1|bat_flying) + (1|bat_roosting) + (1|group), data= d, family= poisson)
check_overdispersion(t)

# try NBM
t <- glmmTMB(n_response ~ flight_time2 + log_inquiry+ (1|bat_flying) + (1|bat_roosting) + (1|group),
             data=d,
             ziformula=~1,             
             family=nbinom2)
check_overdispersion(t)
check_zeroinflation(t)
t2 <- tidy(t,conf.int=TRUE,exponentiate=T,effects="fixed", conf.method="profile")
t2
# for every minute of flight time, the response call count increases by a factor of 1.58 [1.26, 1.95]

##### Effect of inquiry calling on response calling------------

# plot with log counts
d %>% 
  ggplot(aes(x=log_inquiry, y=log(n_response+1)))+
  geom_point(size=2)+ 
  geom_smooth(method= "lm")+
  xlab("log inquiry call count")+
  ylab("log (x+1) response call count")

# plot negative binomial curve
# plot all bats
(t1 <- 
  d %>% 
  mutate(label= "all pairs") %>% 
  mutate(escape= as.logical(escape)) %>% 
  ggplot(aes(x=log_inquiry, y=n_response))+
    facet_wrap(~label)+
  geom_point(size=2, alpha=0.5)+ 
  geom_smooth(method="glm.nb")+
  xlab("log of inquiry call count")+
  ylab("response call count"))+
  theme(legend.position= 'none')

# plot with kinship
(t2 <- 
  d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(escape= as.logical(escape)) %>% 
  mutate(kin= if_else(kinship>0.1, "kin", "nonkin")) %>% 
  mutate(kin= factor(kin, levels= c("nonkin", "kin"))) %>% 
  ggplot(aes(x=log_inquiry, y=(n_response), color=kin, group=kin))+
  facet_wrap(~kin)+
  geom_point(size=2, alpha=0.5, aes(shape= escape))+ 
  geom_smooth(method="glm.nb")+
  xlab("log of inquiry call count")+
  ylab("response call count")+
  scale_color_manual(values= c("darkgrey", "darkred"))+
  theme(legend.position= "none")+
  theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none"))

# plot as linear model
(t3 <- 
    d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kin= if_else(kinship>0.1, "kin", "nonkin")) %>% 
  mutate(kin= factor(kin, levels= c("nonkin", "kin"))) %>% 
  ggplot(aes(x=log_inquiry, y=log(n_response+1), color=kin, group=kin))+
  facet_wrap(~kin)+
  geom_point(size=2, alpha=0.5)+ 
  geom_smooth(method="lm")+
  xlab("log of inquiry call count")+
  ylab("log of response call count")+
  scale_color_manual(values= c("darkgrey", "darkred"))+
  theme(legend.position= "none"))

(inquiry.response <- t1+t2+plot_layout(ncol = 2) + plot_layout(widths = c(1, 2)))

# save as PDF
ggsave(
  "inquiry_response.pdf",
  plot = inquiry.response,
  width = 14,
  height = 7,
  units = "in",
  dpi = 300)

# for plots, get residuals from negative binomial mixed effect model predicting number of responses by number of inquiry calls
fit <- glmmTMB(n_response ~ log_inquiry + offset(log(flight_time)) + (1|bat_flying) + (1|bat_roosting) + (1|group),
                data=d,
                ziformula=~1,
                family=nbinom2)
# get difference between observed response count and predicted response count
# add this response score to data for plotting
# this value represents response calling more than expected from inquiry calls
d$resid_response <- d$n_response - predict(fit, type= "response", newdata= d, allow.new.levels = T)


# plot number of inquiry calls by kinship with partner
d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kin= as.character(kinship)) %>% 
  ggplot(aes(x=kin, y=n_inquiry, group= kin))+
  geom_violin(fill=NA, width=1)+
  geom_jitter(size=2, width= 0.1, height=0, aes(color= udyad.sex))+
  xlab("kinship")+
  ylab("number of inquiry calls")+
  theme(legend.position= "top", legend.title = element_blank())

# plot number of response calls by kinship with partner
d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kin= as.character(kinship)) %>% 
  ggplot(aes(x=kin, y=n_response, group= kin))+
  geom_violin(fill=NA, width=1)+
  geom_jitter(size=2, width= 0.1, height=0, aes(color= udyad.sex))+
  xlab("kinship")+
  ylab("number of response calls")+
  theme(legend.position= "top", legend.title = element_blank())

# plot residual response call variation 
d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kin= as.character(kinship)) %>% 
  ggplot(aes(x=kin, y=resid_response, group= kin))+
  geom_violin(fill=NA, width=1)+
  geom_jitter(size=2, width= 0.1, height=0, aes(color= udyad.sex))+
  xlab("kinship")+
  ylab("residual response call variation")+
  theme(legend.position= "top", legend.title = element_blank())

# create function to plot variable by kinship with means and 95% CIs
plot_kinship_effect <- function(d=d, y= d$sri, label= 'label') {
  # plot association by kinship
  set.seed(123)
  means <- 
    d %>% 
    mutate(y=y) %>% 
    filter(!is.na(kinship)) %>% 
    group_by(dyad) %>% 
    summarize(kinship= mean(kinship), y= mean(y), udyad.sex= first(udyad.sex)) %>% 
    mutate(kin= ifelse(kinship>0, "kin", "nonkin")) %>% 
    dplyr::select(kin, y) %>% 
    boot_ci2(y= .$y, x=.$kin, bca=T) %>% 
    rename(kin= effect) %>% 
    mutate(kin= factor(kin, levels= c("nonkin", "kin")))
  points <- 
    d %>% 
    mutate(y=y) %>% 
    filter(!is.na(kinship)) %>% 
    group_by(dyad) %>% 
    summarize(kinship= mean(kinship), y= mean(y), udyad.sex= first(udyad.sex)) %>% 
    mutate(kin= ifelse(kinship>0, "kin", "nonkin")) %>% 
    mutate(kin= factor(kin, levels= c("nonkin", "kin")))
  # plot means and 95% CI
  means %>% 
    ggplot(aes(x=kin, y=mean, color= kin))+
    geom_jitter(data= points, aes(y= y), size=1, alpha=0.5, width=0.1)+
    geom_point(position = position_nudge(x = 0.25), size=3)+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), position = position_nudge(x = 0.25), size=1)+
    xlab("")+
    ylab(label) +
    scale_color_manual(values= c("darkgrey", "darkred"))+
    coord_flip()+                     
    theme_bw()+
    theme(legend.position= "none")
}

# plot effect of kinship on SRI
(p1 <- plot_kinship_effect(d,d$sri, label= 'association (SRI)'))

# plot effect of kinship on inquiry calling
(p2 <- plot_kinship_effect(d,d$n_inquiry/(d$flight_time/60), label= 'inquiry calls per min of flight time'))

# plot effect of kinship on response calling
(p3 <- plot_kinship_effect(d,d$n_response/(d$flight_time/60), label= 'response calls per min of flight time'))

# plot effect of kinship on response calling
(p4 <- plot_kinship_effect(d,d$resid_response, label= 'residual response call variation'))

# plot all
(kinship.effects <- p1+p2+p3+p4+plot_layout(ncol = 1))

# save as PDF
ggsave(
  "kinship_effects.pdf",
  plot = kinship.effects,
  scale = 1,
  width = 3.5,
  height = 7,
  units = "in",
  dpi = 300)

# plot number of inquiry calls by association and kinship (by dyad type)
t1 <- 
  d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kin= ifelse(kinship>0.1, "kin", "nonkin")) %>% 
  ggplot(aes(x=sri, y=log_inquiry))+
  facet_wrap(~dyad.sex, scales= "free_y")+
  geom_point(aes(color=kin), size=2, alpha=0.6)+
  geom_smooth(method= "lm")+
  xlab("association rate (simple ratio index)")+
  ylab("inquiry call count (log transformed)")+
  scale_colour_manual(values= c("darkred", "grey"))+
  theme_bw()+
  theme(legend.position= "none", legend.title = element_blank())

# plot number of response calls by association and kinship (by dyad type)
t2 <- 
  d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kin= ifelse(kinship>0.1, "kin", "nonkin")) %>% 
  ggplot(aes(x=sri, y=log(n_response+1)))+
  facet_wrap(~dyad.sex, scales= "free_y")+
  geom_point(aes(color=kin), size=2, alpha=0.6)+
  geom_smooth(method= "lm")+
  xlab("association rate (simple ratio index)")+
  ylab("response call count (log transformed)")+
  scale_colour_manual(values= c("darkred", "grey"))+
  theme_bw()+
  theme(legend.position= "none", legend.title = element_blank())
 
# plot residual response calling by association (by dyad type)- controlling for inquiry calls
t3 <- 
  d %>% 
   filter(!is.na(kinship)) %>% 
   mutate(kin= ifelse(kinship>0.1, "kin", "nonkin")) %>% 
   ggplot(aes(x=sri, y=resid_response))+
   facet_wrap(~dyad.sex, scales= "free_y")+
   geom_point(aes(color=kin), size=2, alpha=0.6)+
   geom_smooth(method= "lm")+
   xlab("association rate (simple ratio index)")+
   ylab("residual response call variation")+
   scale_colour_manual(values= c("darkred", "grey"))+
  theme_bw()+
   theme(legend.position= "none", legend.title = element_blank())

# combine plots
(by.sex.plot <- t1+t2+t3+plot_layout(ncol = 1)+ plot_annotation(tag_levels = 'A'))

# save as PDF
ggsave(
  "by_sex_plot.pdf",
  plot = by.sex.plot,
  scale = 1,
  width = 7,
  height = 14,
  units = "in",
  dpi = 300)

######## Effect of kinship and association on response calling---------

# poisson model is overdispersed
t <- glmer(n_response ~ kinship*sri + log_inquiry + offset(log(flight_time)) +  (1|bat_flying) + (1|bat_roosting) + (1|group), family= poisson, data = d)
check_overdispersion(t)

# fit negative binomial model (NBM)
t <- glmmTMB(n_response ~ kinship*sri + log_inquiry + offset(log(flight_time)) + (1|bat_flying) + (1|bat_roosting) + (1|group),
               data=d,
               ziformula=~0,
               family=nbinom2)
check_overdispersion(t)
check_zeroinflation(t, tolerance = 0.05)
AIC(t) # 2700
BIC(t) # 2734

# fit zero-inflated NBM
fit <- glmmTMB(n_response ~ kinship*sri + log_inquiry + offset(log(flight_time)) + (1|bat_flying) + (1|bat_roosting) + (1|group),
               data=d,
               ziformula=~1,
               family=nbinom2)
AIC(fit) #2649
BIC(fit) # 2687
summary(fit)

# fit zero-inflated NBM without interaction
fit <- glmmTMB(n_response ~ kinship + sri + log_inquiry + offset(log(flight_time)) + (1|bat_flying) + (1|bat_roosting) + (1|group),
                data=d,
                ziformula=~1,
                family=nbinom2)
AIC(fit) #2650
BIC(fit) # 2684
summary(fit)

# get predicted values
d$n_response_predicted <- predict(fit, type= "response", newdata= d, allow.new.levels = T)

# plot model performance
d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kinship= kinship>0) %>% 
  ggplot(aes(x=n_response, y=n_response_predicted))+
    geom_point(size=2, aes(color=kinship))+
    geom_smooth(method= "lm")+
    xlab("observed count of response calls")+
  ylab("predicted count of response calls")+
  scale_color_manual(values= c("darkgrey", "darkred"))+
  theme(legend.position= "top")
# get relationship between predicted and observed
summary(lm(n_response_predicted~n_response, data=d)) # r-squared= 0.46

# get model coefficients with 95% CIs
coefs <- 
  tidy(fit,conf.int=TRUE,exponentiate=F,effects="fixed", conf.method="profile") %>% 
  filter(component=="cond") %>% 
  dplyr::select(-effect, -component) %>% 
  mutate(type= "full model")
coefs

# save model results
write.csv(coefs, file="response_model_results.csv")

# effect of SRI with kin only
fitk <- glmmTMB(n_response ~ sri+log_inquiry + offset(log(flight_time)) + (1|bat_flying) + (1|bat_roosting) + (1|group),
                data=d %>% filter(kinship>0),
                ziformula=~1,
                family=nbinom2)
coefs.k <- 
  tidy(fitk,conf.int=TRUE,exponentiate=F,effects="fixed", conf.method="profile") %>% 
  filter(component=="cond") %>% 
  dplyr::select(-effect, -component) %>% 
  mutate(type= "kin only")
coefs.k

# effect of SRI with nonkin only
fitn <- glmmTMB(n_response ~ sri + log_inquiry + offset(log(flight_time))  + (1|bat_flying) + (1|bat_roosting) + (1|group),
                data=d %>% filter(kinship==0),
                ziformula=~1,
                family=nbinom2)
# get 95% CI
coefs.n <- 
  tidy(fitn,conf.int=TRUE,exponentiate=F,effects="fixed", conf.method="profile") %>% 
  filter(component=="cond") %>% 
  dplyr::select(-effect, -component) %>% 
  mutate(type= "nonkin only")
coefs.n

# plot model results
theme_set(theme_bw(base_size = 12))
plot1 <- 
  coefs %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term= case_when(
    term == 'sri' ~ "association",
    term == 'kinship' ~ "kinship",
    term == 'log_inquiry' ~ "log count of inquiry calls",
    term == 'kinship:sri' ~ "association x kinship interaction")) %>% 
  mutate(term= factor(term, levels =c("log count of inquiry calls",
                                      "association x kinship interaction",
                                      "association", 
                                      "kinship"))) %>% 
  mutate(term = fct_rev(term)) %>% 
  ggplot(aes(x=estimate, y=term))+
  geom_point(size=2)+
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high, height=0.2), size=1)+
  geom_vline(xintercept = 0, linetype= "dashed")+
  ylab("")+
  xlab("coefficient")+
  coord_cartesian(xlim= c(-1.5, 1.5))+
  theme(axis.text=element_text(size=12), strip.text = element_text(size=12))
plot1

# save as PDF
ggsave(
  "response_model_results.pdf",
  plot = plot1,
  scale = 1,
  width = 5,
  height = 2,
  units = "in",
  dpi = 300)


# check that sri and kinship do not predict response calls as single predictors 
# fit zero-inflated negative binomial model with only sri
t <- glmmTMB(n_response ~ sri + log_inquiry+ offset(log(flight_time))  + (1|bat_flying) + (1|bat_roosting) + (1|group),
               data=d,
               ziformula=~1,
               family=nbinom2)
summary(t)

# fit zero-inflated negative binomial model with only kinship
t <- glmmTMB(n_response ~ kinship +  log_inquiry+  offset(log(flight_time)) + (1|bat_flying) + (1|bat_roosting) + (1|group),
               data=d,
               ziformula=~1,
               family=nbinom2)
summary(t)

######## Effect of kinship and association on inquiry calling---------

# poisson model is overdispersed
t <- glmer(n_inquiry ~ kinship*sri + log(n_response+1) + offset(log(flight_time)) +  (1|bat_flying) + (1|bat_roosting) + (1|group), family= poisson, data = d)
check_overdispersion(t)

# fit negative binomial model (NBM) with log responses
t <- glmmTMB(n_inquiry ~ kinship*sri + log(n_response+1) + offset(log(flight_time)) +  (1|bat_flying) + (1|bat_roosting) + (1|group),
             data=d,
             ziformula=~0,
             family=nbinom2)
check_overdispersion(t)
AIC(t) # 3130
BIC(t) # 3163

# fit NBM with responses
t<- glmmTMB(n_inquiry ~ kinship*sri + n_response + offset(log(flight_time)) +  (1|bat_flying) + (1|bat_roosting) + (1|group),
             data=d,
             ziformula=~0,
             family=nbinom2)
check_overdispersion(t)
AIC(t) #3130
BIC(t) # 3164
summary(t)

# fit NBM with responses without interaction
fit2 <- glmmTMB(n_inquiry ~ kinship + sri + n_response + offset(log(flight_time)) +  (1|bat_flying) + (1|bat_roosting) + (1|group),
                data=d,
                ziformula=~0,
                family=nbinom2)
AIC(fit2) #3129
BIC(fit2) # 3159
summary(fit2)
# get predicted values
d$n_inquiry_predicted <- predict(fit2, type= "response", newdata= d, allow.new.levels = T)

# plot model performance
d %>% 
  filter(!is.na(kinship)) %>% 
  mutate(kinship= kinship>0) %>% 
  ggplot(aes(x=n_inquiry, y=n_inquiry_predicted))+
  geom_point(size=2, aes(color=kinship))+
  geom_smooth(method= "lm")+
  xlab("observed count of inquiry calls")+
  ylab("predicted count of inquiry calls")+
  scale_color_manual(values= c("darkgrey", "darkred"))+
  theme(legend.position= "top")
# get relationship between predicted and observed
summary(lm(n_inquiry_predicted~n_inquiry, data=d)) # r-squared= 0.71

# get model coefficients
summary(fit2)
coefs2 <- 
  tidy(fit2,conf.int=TRUE,exponentiate=F,effects="fixed", conf.method="profile") %>% 
  filter(component=="cond") %>% 
  dplyr::select(-effect, -component) %>% 
  mutate(type= "full model")
coefs2

# save model results
write.csv(coefs2, file="inquiry_model_results.csv")

# plot model results
theme_set(theme_bw(base_size = 12))
plot2 <- 
  coefs2 %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term= case_when(
    term == 'sri' ~ "association",
    term == 'kinship' ~ "kinship",
    term == 'n_response' ~ "response calls")) %>% 
  mutate(term= factor(term, levels =c("response calls",
                                      "association", 
                                      "kinship"))) %>% 
  mutate(term = fct_rev(term)) %>% 
  ggplot(aes(x=estimate, y=term))+
  geom_point(size=2)+
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high, height=0.2), size=1)+
  geom_vline(xintercept = 0, linetype= "dashed")+
  ylab("")+
  xlab("estimate (log odds)")+
  theme(axis.text=element_text(size=12), strip.text = element_text(size=12))
plot2


# save as PDF
ggsave(
  "inquiry_model_results.pdf",
  plot = plot1,
  scale = 1,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300)

# get relative amount of variance explained by random intercept
tidy(fit2,conf.int=F,exponentiate=F) %>% 
  filter(effect== "ran_pars") %>% 
  mutate(variance= estimate^2) %>% 
  dplyr::select(effect, group, variance) %>% 
  group_by(effect) %>% 
  mutate(sum= sum(variance)) %>% 
  mutate(prop= variance/sum)

# confirm no effect with single predictors
# association
t <- glmmTMB(n_inquiry ~ sri + n_response + offset(log(flight_time)) +  (1|bat_flying) + (1|bat_roosting) + (1|group),
             data=d,
             ziformula=~0,
             family=nbinom2)
summary(t)

# kinship
t <- glmmTMB(n_inquiry ~ kinship + n_response + offset(log(flight_time)) +  (1|bat_flying) + (1|bat_roosting) + (1|group),
                data=d,
                ziformula=~0,
                family=nbinom2)
summary(t)
  

