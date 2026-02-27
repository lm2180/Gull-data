# Gull-data
Uses general linear mixed models to analyse data collected on gull behaviour. Code for R studio.
#Install packages
install.packages("ggplot2")
install.packages("maps")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("psych")
install.packages("tidyr")
install.packages("lme4")
install.packages("ggeffects")
library(ggplot2)
library(tidyverse)
library(sf)
library(dplyr)
library(psych)
library(tidyr)
library(lme4)
#Citations
sessionInfo()
citation()
citation("ggplot2")
citation("tidyverse")
citation("ggspatial")
citation("dplyr")
citation("psych")
citation("tidyr")
citation("lme4")
#Set up working directory and retrieve raw data
setwd("C:/Users/jizza/OneDrive/Desktop/Uni/Year 3/Bio 350 Dissertaion project/Seagull/GLMM data")
glmmdata<-read.csv("gull_obs_GLMM.csv",header=T)
view(glmmdata)
#remove NA and OOS data points
table(is.na(glmmdata$Behaviour))
unique(glmmdata$Behaviour)
glmmdata$Behaviour[glmmdata$Behaviour==""]<-NA
glmmdata$Behaviour[glmmdata$Behaviour=="OOS"]<-NA
gd2<-glmmdata%>%
  filter(!is.na(Behaviour))
table(is.na(glmmdata$Behaviour))
unique(gd2$Behaviour)
str(gd2)
view(gd2)
#Create a binary response variable for Resting
gd2$Resting<-ifelse(gd2$Behaviour=="R",1,0)
table(gd2$Resting)
gd2$Individual.ID<-as.factor(gd2$Individual.ID)
gd2$Area<-as.factor(gd2$Area)
#GLMM for resting
m_rest<-glmer(Resting~People.Count*Area+(1|Individual.ID),family=binomial,data=gd2)
summary(m_rest)
#outliers
install.packages("performance")
install.packages("see")
library(performance)
library(see)
check_model(m_rest,check="outliers") # leverage plot within the contour lines.
check_outliers(m_rest)#No outliers
#Check for over-dispersion
overdisp_fun<-function(model){
  rdf<-df.residual(model)
  rp<-residuals(model,type="pearson")
  Pearson.chisq<-sum(rp^2)
  prat<-Pearson.chisq/rdf
  pval<-pchisq(Pearson.chisq, df=rdf,lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m_rest)#Ratio<1 = 0.952, p=0.98
#Create a binary response variable for Locomotion
gd2$Locomotion<-ifelse(gd2$Behaviour=="L",1,0)
table(gd2$Locomotion)
#GLMM for Locomotion
m_loc<-glmer(Locomotion~People.Count*Area+(1|Individual.ID),family=binomial,data=gd2)
summary(m_loc)
#outliers
check_model(m_loc,check="outliers") # leverage plot within the contour lines.
check_outliers(m_loc)#3 outliers
#Check for over-dispersion
overdisp_fun<-function(model){
  rdf<-df.residual(model)
  rp<-residuals(model,type="pearson")
  Pearson.chisq<-sum(rp^2)
  prat<-Pearson.chisq/rdf
  pval<-pchisq(Pearson.chisq, df=rdf,lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m_loc)#Ratio<1 = 0.911, p=1.0
#Create a binary response variable for Foraging
gd2$Foraging<-ifelse(gd2$Behaviour=="F",1,0)
table(gd2$Foraging)
#GLMM for Foraging
m_for<-glmer(Foraging~People.Count*Area+(1|Individual.ID),family = binomial,data = gd2)
summary(m_for)
#outliers
check_model(m_for,check="outliers") # leverage plot within the contour lines.
check_outliers(m_for)
#Check for over-dispersion
overdisp_fun<-function(model){
  rdf<-df.residual(model)
  rp<-residuals(model,type="pearson")
  Pearson.chisq<-sum(rp^2)
  prat<-Pearson.chisq/rdf
  pval<-pchisq(Pearson.chisq, df=rdf,lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m_for)#Ratio<1 = 0.671, p=1.0
#Create a binary response variable for Vigilance
gd2$Vigilance<-ifelse(gd2$Behaviour=="V",1,0)
table(gd2$Vigilance)
#GLMM for Vigilance
m_vig<-glmer(Vigilance~People.Count*Area+(1|Individual.ID),family=binomial,data=gd2)
summary(m_vig)
#model unidentifiable - re scale variables
table(gd2$Area,gd2$Vigilance)
gd2$People_sc<-scale(gd2$People.Count)
m_vig<-glmer(Vigilance~People_sc*Area+(1|Individual.ID),family=binomial,data=gd2,
             control=glmerControl(optimizer="bobyqa"))
summary(m_vig)
#Create a binary response variable for Maintenance
gd2$Maintenance<-ifelse(gd2$Behaviour=="G",1,0)
table(gd2$Maintenance)
#GLMM for Maintenance
m_mai<-glmer(Maintenance~People.Count*Area+(1|Individual.ID),family=binomial,data=gd2)
summary(m_mai)
#outliers
check_model(m_mai,check="outliers") # leverage plot within the contour lines.
check_outliers(m_mai)
#Check for over-dispersion
overdisp_fun<-function(model){
  rdf<-df.residual(model)
  rp<-residuals(model,type="pearson")
  Pearson.chisq<-sum(rp^2)
  prat<-Pearson.chisq/rdf
  pval<-pchisq(Pearson.chisq, df=rdf,lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m_mai)#Ratio<1 = 0.744, p=1.0
#Create a binary response variable for Social
gd2$Social<-ifelse(gd2$Behaviour=="S",1,0)
table(gd2$Social)
#GLMM for Social
m_soc<-glmer(Social~People.Count*Area+(1|Individual.ID),family=binomial,data=gd2)
summary(m_soc)
#outliers
check_model(m_soc,check="outliers") # leverage plot within the contour lines.
check_outliers(m_soc)
#Check for over-dispersion
overdisp_fun<-function(model){
  rdf<-df.residual(model)
  rp<-residuals(model,type="pearson")
  Pearson.chisq<-sum(rp^2)
  prat<-Pearson.chisq/rdf
  pval<-pchisq(Pearson.chisq, df=rdf,lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m_soc)#Ratio<1 = 0.296, p=1.0
#Create a binary response variable for Antagonistic
gd2$Antagonistic<-ifelse(gd2$Behaviour=="A",1,0)
table(gd2$Antagonistic)
#GLMM for Antagonistic
m_ant<-glmer(Antagonistic~People.Count*Area+(1|Individual.ID),family=binomial,data=gd2)
summary(m_ant)
#Needs rescaling
table(gd2$Area,gd2$Antagonistic)
#Load in data set for sub-category 
glmmdata2<-read.csv("gull data collection GLMM.csv",header=T)
#remove NA and OOS
table(is.na(glmmdata2$Behaviour))
unique(glmmdata2$Behaviour)
glmmdata2$Behaviour[glmmdata2$Behaviour==""]<-NA
glmmdata2$Behaviour[glmmdata2$Behaviour=="OOS"]<-NA
behavioursections<-glmmdata2%>%
  filter(!is.na(Behaviour))
table(is.na(glmmdata2$Behaviour))
unique(behavioursections$Behaviour)
str(behavioursections)
view(behavioursections)
#LOCOMOTION
#Subset data for walking
glmmdata2$Locomotion_walking<-ifelse(glmmdata2$Behaviour=="L1",1,0)
table(glmmdata2$Locomotion_walking)
glmmdata2$Individual.ID<-as.factor(glmmdata2$Individual.ID)
glmmdata2$Area<-as.factor(glmmdata2$Area)
#GLMM for walking
m_walk<-glmer(Locomotion_walking~People.Count*Area+(1|Individual.ID),family=binomial,data=glmmdata2)
summary(m_walk)
#Subset data for flying
glmmdata2$Locomotion_flying<-ifelse(glmmdata2$Behaviour=="L2",1,0)
table(glmmdata2$Locomotion_flying)
#GLMM for flying
m_fly<-glmer(Locomotion_flying~People.Count*Area+(1|Individual.ID),family=binomial,data=glmmdata2)
summary(m_fly)
#Subset data for gliding
glmmdata2$Locomotion_gliding<-ifelse(glmmdata2$Behaviour=="L3",1,0)
table(glmmdata2$Locomotion_gliding)
#GLMM for gliding
m_gli <- glmer(Locomotion_gliding~People.Count*Area+(1|Individual.ID),family=binomial,data=glmmdata2)
summary(m_gli)
table(glmmdata2$Area,glmmdata2$Locomotion_gliding)
#FORAGING
#Subset data for pecking
glmmdata2$Foraging_pecking<-ifelse(glmmdata2$Behaviour=="F1",1,0)
table(glmmdata2$Foraging_pecking)
#GLMM for pecking
m_pec<-glmer(Foraging_pecking~People.Count*Area+(1|Individual.ID),family=binomial,data=glmmdata2)
summary(m_pec)
#Subset data for foot paddle
glmmdata2$Foraging_foot<-ifelse(glmmdata2$Behaviour=="F2",1,0)
table(glmmdata2$Foraging_foot)
#GLMM for paddle
m_foot<-glmer(Foraging_foot~People.Count*Area+(1|Individual.ID),family=binomial,data=glmmdata2)
summary(m_foot)
table(glmmdata2$Area,glmmdata2$Foraging_foot)
glmmdata2$Area2<-ifelse(glmmdata2$Area=="Residential","Urban")
glmm_sub<-subset(glmmdata2,Area%in%c("Beach","Singleton"))
m_foot<-glmer(Foraging_foot~People.Count*Area+(1|Individual.ID),family=binomial,data=glmm_sub)
summary(m_foot)
#Subset data for kleptoparasitism
glmmdata2$Foraging_klepto<-ifelse(glmmdata2$Behaviour=="F3",1,0)
table(glmmdata2$Foraging_klepto)
#GLMM for klepto
m_kle<-glmer(Foraging_klepto~People.Count*Area+(1|Individual.ID),family=binomial,data=glmmdata2)
summary(m_kle)
#Data is unstable.
table(glmmdata2$Area,glmmdata2$Foraging_klepto)
#Some Areas had 0 activity
glmm_sub2<-subset(glmmdata2,Area%in%c("Beach","Singleton"))
m_kle<-glmer(Foraging_klepto~People.Count*Area+(1|Individual.ID),family=binomial,data=glmm_sub2)
summary(m_kle)
#Subset data for scavenging
glmmdata2$Foraging_scavenging<-ifelse(glmmdata2$Behaviour=="F4",1,0)
table(glmmdata2$Foraging_scavenging)
#GLMM for scavenging
m_sca<-glmer(Foraging_scavenging~People.Count*Area+(1|Individual.ID),family=binomial,data=glmmdata2)
summary(m_sca)
#Data is unstable
table(glmmdata2$Area,glmmdata2$Foraging_scavenging)
glmm_sub3<-subset(glmmdata2,Area%in%c("City Centre","Residential","Singleton"))
m_sca<-glmer(Foraging_scavenging~People.Count*Area+(1|Individual.ID),family=binomial,data=glmm_sub3)
summary(m_sca)
#RESTING
#Subset data for standing
glmmdata2$Resting_stand<-ifelse(glmmdata2$Behaviour=="R1",1,0)
table(glmmdata2$Resting_stand)
#GLMM for standing
m_sta<-glmer(Resting_stand~People.Count*Area+(1|Individual.ID),family=binomial,data=glmmdata2)
summary(m_sta)
#Subset data for sitting
glmmdata2$Resting_sit<-ifelse(glmmdata2$Behaviour=="R2",1,0)
table(glmmdata2$Resting_sit)
#GLMM for sitting
m_sit <- glmer(Resting_sit~People.Count*Area+(1|Individual.ID),family=binomial,data=glmmdata2)
summary(m_sit)
table(glmmdata2$Area, glmmdata2$Resting_sit)
glmm_sub4<-subset(glmmdata2,Area%in%c("City Centre","Beach"))
m_sit<-glmer(Resting_sit~People.Count*Area+(1|Individual.ID),family=binomial,data=glmm_sub4)
summary(m_sit)
#Subset data for Perch
glmmdata2$Resting_perch<-ifelse(glmmdata2$Behaviour=="R3",1,0)
table(glmmdata2$Resting_perch)
#GLMM for perching
m_per<-glmer(Resting_perch~People.Count*Area+(1|Individual.ID),family=binomial,data = glmmdata2)
summary(m_sta)
table(glmmdata2$Area,glmmdata2$Resting_perch)#No beach results
glmm_sub5<-subset(glmmdata2,Area%in%c("City Centre","Residential","Singleton"))
m_per<-glmer(Resting_perch~People.Count*Area+(1|Individual.ID),family=binomial,data=glmm_sub5)
summary(m_per)
#Descriptive stastistics
#People count
people_summary<-gd2%>%
  group_by(Area)%>%
  summarise(mean_people=mean(People.Count,na.rm=TRUE),
            sd_people=sd(People.Count,na.rm=TRUE),
            max_people=max(People.Count,na.rm=TRUE))
view(people_summary)
#Figures
#GML Behaviours
glm_plot_data<-gd2%>%
  pivot_longer(cols=c(Resting,Locomotion,Foraging,Maintenance,Vigilance,Antagonistic,Social),
               names_to="Behaviours",
               values_to="value")
ggplot(glm_plot_data,aes(x=People,y=value))+
  geom_jitter(height=0.05,width=0,alpha=0.2)+
  geom_smooth(method="glm",method.args=list(family="binomial"),
              se=TRUE)+facet_wrap(~Behaviours,scales="free_y")+
  labs(x="Number of people present",y="Probability of behaviour")+
  theme(strip.text=element_text(size=18,face="bold"),axis.title=element_text(size=18),
        axis.text=element_text(size=18),legend.title=element_text(size=18),
        legend.text=element_text(size=18))+
  theme_classic()
#People per area
ggplot(glmmdata,aes(x=Area,y=People.Count,fill=Area))+
  geom_boxplot()+
  scale_fill_manual(values=c("Beach"="steelblue","City Centre"="plum1","Residential"="seagreen2","Singleton"="red"))+
  labs(y="Number of People",x="Area")+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))
pred_rest<-ggpredict(m_rest,terms=c("People.Count [all]","Area"),bias_correction=TRUE)
y
pred_loc<-ggpredict(m_loc,terms=c("People.Count [all]","Area"),bias_correction=TRUE)
ggplot()+
  geom_jitter(data=gd2,aes(x=People.Count,y=Resting,colour=Area),alpha=0.15,height=0.03)+
  geom_line(data=pred_rest,aes(x=x,y=predicted,colour=group),linewidth=1.2)+
  labs(x="Number of People",y="Probability of Resting")+
  theme_classic(base_size=14)
ggplot()+
  geom_jitter(data=gd2,aes(x=People.Count,y=Locomotion,colour=Area),alpha=0.15,height=0.03)+
  geom_line(data=pred_loc,aes(x=x,y=predicted,colour=group),linewidth=1.2)+
  labs(x="Number of People",y="Probability of Locomotion")+
  theme_classic(base_size=14)


