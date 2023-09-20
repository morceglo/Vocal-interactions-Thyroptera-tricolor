require(lme4)

#Data manipulation
## Joining the two data sets
d1<-read.csv("vocal_interactions_2021.csv", sep=";")
d2<-read.csv("vocal_interactions_2022.csv", sep=";")

colnames(d2)[6]<-"Dyad"
d<-rbind(d1,d2)

## Adding observation level random effect
d$ob<-1:nrow(d)

#Stats models
##Model for inquiry with all observations
mod_inquiry<-glmer(n_inquiry~ 1 + (1|Bat_roosting)+ (1|Bat_flying) + (1|Group) + (1|ob), data=d, family="poisson")

###Extracting variance components
VBF_I<-VarCorr(mod_inquiry)$Bat_flying[1,1] #Variance Bats Flying
VBR_I<-VarCorr(mod_inquiry)$Bat_roosting[1,1] #Variance Bats Roosting
VBG_I<-VarCorr(mod_inquiry)$Group[1,1] #Variance Group
VOD_I<-VarCorr(mod_inquiry)$ob[1,1] #Variance Overddispersion

PV_I<-log(1/exp(fixef(mod_inquiry)[1]) + 1)  #Variance poisson processs

####Repeatabilty roosting individual
VBR_I/(VBR_I+VOD_I +VBF_I+PV_I + VBG_I)

####Repeatabilty flying individual
VBF_I/(VBR_I+VOD_I +VBF_I+PV_I+ VBG_I)


##Model for inquiry only when someone responded
d2<-d[d$n_response>0,]
mod_inquiry2<-glmer(n_inquiry~ 1 + (1|Bat_roosting)+ (1|Bat_flying) + (1|ob), data=d2, family="poisson")

###Extracting variance components
VBF_I<-VarCorr(mod_inquiry2)$Bat_flying[1,1] #Variance Bats Flying
VBR_I<-VarCorr(mod_inquiry2)$Bat_roosting[1,1] #Variance Bats Roosting
VBG_I<-VarCorr(mod_inquiry)$Group[1,1]
VOD_I<-VarCorr(mod_inquiry2)$ob[1,1] #Variance Overddispersion
PV_I<-log(1/exp(fixef(mod_inquiry2)[1]) + 1)  #Variance poisson processs

####Repeatabilty roosting individual
VBR_I/(VBR_I+VOD_I +VBF_I+PV_I+ VBG_I)

####Repeatabilty flying individual
VBF_I/(VBR_I+VOD_I +VBF_I+PV_I+ VBG_I)

##Model for response 
mod_response<-glmer(n_response~ 1+ (1|Bat_roosting)+ (1|Bat_flying) + (1|Group) + (1|ob), data=d, family="poisson")
summary(mod_response)

###Extracting variance components
VBF_R<-VarCorr(mod_response)$Bat_flying[1,1] #Variance Bats Flying
VBR_R<-VarCorr(mod_response)$Bat_roosting[1,1] #Variance Bats Roosting
VBG_R<-VarCorr(mod_response)$Group[1,1]
VOD_R<-VarCorr(mod_response)$ob[1,1] #Variance Overddispersion
PV_R<-log(1/exp(fixef(mod_response)[1]) + 1) #Variance poisson processs

##Repeatabilty roosting individual
VBR_R/(VBR_R + VOD_R + VBF_R + PV_R + VBG_R)

##Repeatabilty flying individual
VBF_R/(VBR_R+ VOD_R + VBF_R + PV_R + VBG_R)



