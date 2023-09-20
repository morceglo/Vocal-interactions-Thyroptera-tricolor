require(lme4)
setwd()
#Data manipulation
## Joining the two data sets
d1<-read.csv("./data/raw/vocal_interactions_2021.csv", sep=";")
d2<-read.csv("./data/raw/vocal_interactions_2022.csv", sep=";")
colnames(d2)[6]<-"Dyad"
d<-rbind(d1,d2)

## Adding observation level random effect
d$ob<-1:nrow(d)

#Stats models
##Model for inquiry
mod_inquiry<-glmer.nb(n_inquiry~ 1 + (1|Bat_roosting) + (1|Bat_flying) +  (1|Group), data=d)
summary(mod_inquiry)

###Extracting variance components
mu_l<-fixef(mod_inquiry)[1]
mu<-exp(mu_l)
theta<-getME(mod_inquiry, "glmer.nb.theta")

VBF_I_l<-VarCorr(mod_inquiry)$Bat_flying[1,1] #Variance Bats Flying
VBF_I<- (exp(VBF_I_l)-1)*exp(2*mu_l + VBF_I_l)

VBR_I_l<-VarCorr(mod_inquiry)$Bat_roosting[1,1] #Variance Bats Roosting
VBR_I<-(exp(VBR_I_l)-1)* exp(2*mu_l + VBR_I_l)

VBG_I_l<-VarCorr(mod_inquiry)$Group[1,1] #Variance GRoup
VBG_I<-(exp(VBG_I_l)-1)* exp(2*mu_l + VBG_I_l)


NBV_I<-mu + mu^2/theta

####Repeatabilty roosting individual
VBR_I/(VBR_I+VBF_I+NBV_I + VBG_I)

####Repeatabilty flying individual
VBF_I/(VBR_I+VBF_I+NBV_I)

#Stats models
##Model for inquiry2
d2<-d[d$n_response>0,]
mod_inquiry<-glmer.nb(n_inquiry~ 1 + (1|Bat_roosting)+ (1|Bat_flying) +  (1|Group), data=d2)
summary(mod_inquiry)

###Extracting variance components
mu_l<-fixef(mod_inquiry)[1]
mu<-exp(mu_l)
theta<-getME(mod_inquiry, "glmer.nb.theta")

VBF_I_l<-VarCorr(mod_inquiry)$Bat_flying[1,1] #Variance Bats Flying
VBF_I<- (exp(VBF_I_l)-1)*exp(2*mu_l + VBF_I_l)

VBR_I_l<-VarCorr(mod_inquiry)$Bat_roosting[1,1] #Variance Bats Roosting
VBR_I<-(exp(VBR_I_l)-1)* exp(2*mu_l + VBR_I_l)

VBG_I_l<-VarCorr(mod_inquiry)$Group[1,1] #Variance GRoup
VBG_I<-(exp(VBG_I_l)-1)* exp(2*mu_l + VBG_I_l)


NBV_I<-mu + mu^2/theta

####Repeatabilty roosting individual
VBR_I/(VBR_I+VBF_I+NBV_I)

####Repeatabilty flying individual
VBF_I/(VBR_I+VBF_I+NBV_I)


##Model for response
mod_response<-glmer.nb(n_response~ 1 + (1|Bat_roosting)+ (1|Bat_flying) +  (1|Group), data=d)
summary(mod_response)

###Extracting variance components
mu_l<-fixef(mod_response)[1] #+ fixef(mod_response)[2]*mean(d$n_inquiry, na.rm=TRUE)
mu<-exp(mu_l)
theta<-getME(mod_response, "glmer.nb.theta")

VBF_R_l<-VarCorr(mod_response)$Bat_flying[1,1] #Variance Bats Flying
VBF_R<- (exp(VBF_R_l)-1)*exp(2*mu_l + VBF_R_l)

VBR_R_l<-VarCorr(mod_response)$Bat_roosting[1,1] #Variance Bats Roosting
VBR_R<-(exp(VBR_R_l)-1)* exp(2*mu_l + VBR_R_l)

VBG_R_l<-VarCorr(mod_response)$Group[1,1] #Variance GRoup
VBG_R<-(exp(VBG_R_l)-1)* exp(2*mu_l + VBG_R_l)


NBV_R<-mu + mu^2/theta

####Repeatabilty roosting individual
VBR_R/(VBR_R+VBF_R+NBV_R +VBG_R)

####Repeatabilty flying individual
VBF_R/(VBR_R+VBF_R+NBV_R+VBG_R)



