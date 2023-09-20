# Does kinship predict vocal interaction rates in Thyroptera tricolor?
# Gerald Carter, gcarter1640@gmail.com

# set working directory
setwd("/Users/gerry/Dropbox/Dropbox/_working/_ACTIVE/_collaborators/Gloriana/thyroptera-vocal-interaction/code")

############## CLEAN DATA
# clear workspace
rm(list=ls())

# load packages
library(tidyverse)
library(asnipe)

# set ggplot theme
theme_set(theme_bw(base_size = 14))

### clean and wrangle data

# create function to convert matrix to data-frame
matrix_to_df <- function(m1){
  data.frame(dyad = paste(rownames(m1)[col(m1)], colnames(m1)[row(m1)], sep="_"),
             value = c(t(m1)), stringsAsFactors = FALSE)
}

# get kinship data
k <- as.matrix(read.csv('KinshipMatrix.Sept14.csv', sep= ",", row.names = 1))

# check symmetry
mean(k==t(k), na.rm=T)

# get sex data
sex1 <- read.csv('supp_21.csv', sep= ";") 
sex2 <- read.csv('supp_22.csv', sep= ";") 

# get association data
a21 <- read.csv('associations_2021.csv', sep= ";")
a22 <- read.csv('associations_2022.csv', sep= ";")

# get calling data
i21 <- read.csv('vocal_interactions_2021.csv', sep= ";")
i22 <- read.csv('vocal_interactions_2022.csv', sep= ";")

# tidy sex data
sex <- 
  rbind(sex1,sex2) %>% 
  mutate(bat= paste0('bat', bat_id)) %>% 
  group_by(bat, sex) %>% 
  summarize(n= n()) %>% 
  dplyr::select(bat, sex) %>% 
  ungroup()

# tidy 2021 association data
assoc21 <- 
  a21 %>% 
  pivot_longer(Bat1:Bat10, names_to = 'names', values_to = "bat") %>% 
  filter(!is.na(bat)) %>% 
  mutate(bat= paste0('bat', bat)) %>% 
  mutate(group= paste0('group', Group)) %>%  
  mutate(date= dmy(Date)) %>% 
  dplyr::select(obs= data_id, group, date, bat)

# tidy 2022 association data
assoc22 <- 
  a22 %>% 
  pivot_longer(Bat1:Bat8, names_to = 'names', values_to = "bat") %>% 
  filter(!is.na(bat)) %>% 
  mutate(bat= paste0('bat', bat)) %>% 
  mutate(group= paste0('group', Group)) %>%  
  mutate(date= dmy(Date)) %>% 
  dplyr::select(obs= data_id, group, date, bat)

# combine association data
assoc <- 
  rbind(assoc21, assoc22) %>% 
  mutate(obs= paste(obs,group,date, sep= "_")) %>% 
  dplyr::select(bat, obs) %>% 
  as.data.frame()

# get number of observations per bat
obs_per_bat <- 
  assoc %>% 
  group_by(bat) %>% 
  summarize(n.obs=n()) %>% 
  arrange(desc(n.obs))

# range is 1 to 51 times
range(obs_per_bat$n.obs)

# pick threshold of observations for including bats in analysis
# if they have fewer observations than we make their association rates NA
threshold <- 4

# plot number of observations per bat
assoc %>% 
  group_by(bat) %>% 
  summarize(n.obs=n()) %>% 
  ggplot(aes(x=n.obs))+
    geom_histogram(fill="light blue", color="black")+
    geom_vline(xintercept= threshold, color= 'red', size=1)+
    xlab("number of observations of bat")+
    ylab("count of bats")

# get bats seen fewer than threshold number
low.n.bats <- 
  assoc %>% 
  group_by(bat) %>% 
  summarize(n.obs=n()) %>% 
  filter(n.obs < threshold) %>% 
  pull(bat)
low.n.bats
length(low.n.bats)
# 14 bats were seen fewer than 4 times

### get SRI from asnipe

# get association rates as group-by-individual
gbi <- get_group_by_individual(assoc, data_format="individuals")

# get association network of SRI
assoc.net<- get_network(association_data=gbi, data_format = "GBI", association_index = "SRI")

# check symmetry
mean(assoc.net==t(assoc.net))

# get SRI values for every undirected pair as dataframe
SRI <- matrix_to_df(assoc.net) 

# tidy kinship
kinship <- 
  k %>% 
  matrix_to_df() %>% 
  separate(dyad, into= c("bat1", "bat2")) %>% 
  mutate(bat2= str_remove(bat2, "X")) %>% 
  mutate(bat1= paste0('bat', bat1), bat2= paste0('bat', bat2)) %>% 
  filter(bat1!=bat2) %>% 
  # label undirected pair
  mutate(dyad= ifelse(bat1<bat2, paste(bat1,bat2, sep="_"), paste(bat2,bat1, sep="_"))) %>% 
  group_by(dyad) %>% 
  summarize(kinship= first(value)) 

# count trials
nrow(i21)
nrow(i22)

# trials with missing data
rbind(i21,i22) %>% 
  as_tibble() %>% 
  filter(is.na(n_response)) %>% 
  nrow()
#2

# trials where neither bat called
rbind(i21,i22) %>% 
  as_tibble() %>% 
  filter(n_inquiry==0) %>% 
  nrow()
# 5

# fix and tidy calling data
calling <- 
  rbind(i21,i22) %>% 
  mutate(Date= as.character(mdy(Date))) %>% 
  mutate(group= paste(substr(Date, start=1, stop=4),Group, sep="_")) %>% 
  # need to recreate these labels because messed up by excel in raw data
  separate(Dyad, into=c("bat_flying", "bat_roosting"), remove=F) %>% 
  mutate(bat_flying= paste0('bat', bat_flying), bat_roosting= paste0('bat', bat_roosting)) %>% 
  # label undirected dyads
  mutate(dyad= ifelse(bat_flying<bat_roosting, paste(bat_flying, bat_roosting, sep= "_"), paste(bat_roosting, bat_flying, sep= "_"))) %>% 
  dplyr::select(date= Date, group, bat_flying, bat_roosting, dyad, flight_time, n_inquiry, n_response)

# combine all data
d <- 
  # start with calling data
  calling %>% 
  # add SRI 
  mutate(sri= SRI$value[match(.$dyad, SRI$dyad)]) %>% 
  # add kinship
  mutate(kinship= kinship$kinship[match(.$dyad, kinship$dyad)]) %>% 
  # add sexes of flying and roosting bats
  mutate(flying.sex= sex$sex[match(.$bat_flying, sex$bat)]) %>% 
  mutate(roosting.sex= sex$sex[match(.$bat_roosting, sex$bat)]) %>% 
  # label sex combination for flying-->roosting dyad
  mutate(dyad.sex= paste(flying.sex, roosting.sex, sep= "-->")) %>% 
  # label sex combination for dyad (male, female, both)
  mutate(udyad.sex = case_when(
    flying.sex == 'female' & roosting.sex == 'female' ~ "female-female",
    flying.sex == 'male' & roosting.sex == 'male' ~ "male-male",
    TRUE ~ "mixed-sex")) %>% 
  # replace zero sri (never previously seen together) with NA because association is not 0
  mutate(sri = if_else(sri==0, NA, sri)) %>% 
  # if flying bat seen fewer than 4 times, then SRI is NA
  mutate(sri = if_else(bat_flying %in% low.n.bats, NA, sri)) %>% 
  # if roosting bat seen fewer than 4 times, then SRI is NA
  mutate(sri = if_else(bat_roosting %in% low.n.bats, NA, sri))
           
  
# inspect and fix cases where flying bats in more than 2 groups
t <- 
  d %>% 
  group_by(bat_flying, date, group) %>%
  summarize(n=n()) %>% 
  arrange(date) %>% 
  group_by(bat_flying, group) %>% 
  summarize(n=n()) %>%   
  group_by(bat_flying) %>% 
  summarize(n=n()) %>% 
  filter(n>2) %>% 
  pull(bat_flying)
d %>% 
  dplyr::select(date, group, bat_flying) %>% 
  filter(bat_flying %in% t) %>% 
  arrange(bat_flying) %>% 
  as.data.frame()
# some of these seem impossible and must be errors

# fix impossible group labels
d$group[which(d$date=="2021-07-03" & d$bat_flying=="bat982126051278521")] <- "2021_9"
d$group[which(d$date=="2021-07-03" & d$bat_flying=="bat982126058484300")] <- "2021_9"

# fix cases where roosting bats in more than 2 groups
t <- 
  d %>% 
  group_by(bat_roosting, date, group) %>%
  summarize(n=n()) %>% 
  arrange(date) %>% 
  group_by(bat_roosting, group) %>% 
  summarize(n=n()) %>%   
  group_by(bat_roosting) %>% 
  summarize(n=n()) %>% 
  filter(n>2) %>% 
  pull(bat_roosting)
d %>% 
  dplyr::select(date, group, bat_roosting) %>% 
  filter(bat_roosting %in% t) %>% 
  arrange(bat_roosting) %>% 
  as.data.frame()

# fix group labels
d$group[which(d$date=="2021-07-03" & d$bat_roosting=="bat982126058484300")] <- "2021_9"

# group 9 split into groups 9,1 and 9,2 during 2022
# There are interesting movements between group 9 and 10 and between groups 9,1 and 9,2
# For this analysis, I'm going to consider group 9,1 and 9,2 as the same group
d$group[which(d$group=="2022_9,1" | d$group=="2022_9,2")] <- "2021_9"  

# remove trials without inquiry or response calls
d <- 
  d %>% # 453 rows
  as_tibble() %>% 
  filter(! (is.na(n_inquiry) & is.na(n_response))) %>%   # 451 rows
  filter(n_inquiry>0) %>% # 446 rows
  mutate(year= substr(date, 1,4)) # add year

# save data as csv
write.csv(d, file = paste0('./processed/calling-association-kinship-data.csv'), row.names= F)

