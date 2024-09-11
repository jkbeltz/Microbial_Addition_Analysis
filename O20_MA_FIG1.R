#### FIGURE 1 ####

####Library Packages####
#install.packages("lubridate")
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(cowplot)
library(nlme)
library(lme4)
library(lubridate)
library(lsmeans)

####Import Files####
weather=read.csv("raw/O20_weather_raw.csv")
census=read.csv("raw/O20MA_CS_raw.csv")

##separating date and time##

weather <- separate(weather, TimeStamp, into = c("Date", "Time"), sep = "^\\S*\\K\\s+")
weather$Date = format(as.Date(weather$Date, "%m/%d/%Y"), "20%y-%m-%d")

### create summary df####
#View(weather)
weather.sum=weather %>%
  group_by(Date) %>%
  dplyr::summarise(N=n(),
                   MaxAirTC=max(AirTC_Avg,na.rm = TRUE),
                   MinAIRTC=min(AirTC_Avg,na.rm =TRUE),
                   AirTC=mean(AirTC_Avg,na.rm = TRUE),
                   AirTCsd=sd(AirTC_Avg,na.rm = TRUE),
                   AirTCse= AirTCsd/sqrt(N),
                   Temp_Cavg=mean(Temp_C_Avg,na.rm = TRUE),
                   Temp_Cavgsd=sd(Temp_C_Avg,na.rm = TRUE),
                   Temp_Cavgse= Temp_Cavgsd/sqrt(N),
                   PARden=mean(PAR_Den_Avg,na.rm = TRUE),
                   PARdensd=sd(PAR_Den_Avg,na.rm = TRUE),
                   PARdense= PARdensd/sqrt(N),
                   PARtot=mean(PAR_Tot_Tot,na.rm = TRUE),
                   PARtotsd=sd(PAR_Tot_Tot,na.rm = TRUE),
                   PARtotse= PARtotsd/sqrt(N),
                   WS=mean(WS_ms_Avg,na.rm =TRUE),
                   WSsd=sd(WS_ms_Avg,na.rm = TRUE),
                   WSse= WSsd/sqrt(N),
                   Rain=mean(Rain_mm_Tot,na.rm = TRUE),
                   Rainsd=sd(Rain_mm_Tot,na.rm = TRUE),
                   Rainse=Rainsd/sqrt(N)
  )
wdata=unique(weather.sum)
#View(wdata)

###graph min/max temperature####

AirTCMin <-ggplot(wdata, aes(x=as.factor(Date), y=MaxAirTC, group=1))+
  geom_line(color="Red")+
  #geom_linerange(aes(ymin=AirTC-AirTCse,ymax=AirTC+AirTCse))+
  ylab("Max/Min Daily Temperature C")+
  scale_x_discrete(name=" Experiment Date", 
                   breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
                   labels=c("Start","August", "September", "October", "November", "End")) +
  #xlab("Sample Timepoints")+
  #scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  #scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot()+
  theme(legend.position = "none")

AirTCMinMax<-AirTCMin+geom_line(linetype= "solid" , data=wdata, aes(x=as.factor(Date), y=MinAIRTC, group=1), color="Blue")
AirTCMinMax


####census data####
#View(census)

census$Date=as.Date(census$Date, format = "%m/%d/%y")

Census.sum=census %>%
  group_by(Date, Cage_treatment) %>%
  dplyr::summarise(N=n(),
                   countavg=mean(Picture_total,na.rm = TRUE),
                   countsd=sd(Picture_total,na.rm = TRUE),
                   countse= countsd/sqrt(N),
                   combestavg=mean(est_total,na.rm = TRUE),
                   combestsd=sd(est_total,na.rm = TRUE),
                   combestse=combestsd/sqrt(N^2), 
                   picavgest=mean(est_avg,na.rm = TRUE),
                   picavgsd=sd(est_avg,na.rm = TRUE),
                   picavgse=picavgsd/sqrt(N),
                   pictotavg=mean(est_tot,na.rm = TRUE),
                   pictotsd=sd(est_tot,na.rm = TRUE),
                   pictotse= pictotsd/sqrt(N),
  )
View(Census.sum)

##combine dfs####
Census.sum$Date=as.Date(Census.sum$Date)
wdata$Date=as.Date(wdata$Date)

dplrcombinedtreatment <- left_join(wdata, Census.sum, by="Date")
dplrcombinedtreatment <- dplrcombinedtreatment[! dplrcombinedtreatment$Cage_treatment %in% c('Bloom'), ]
View(dplrcombinedtreatment)
###make plot###
coeff=5000

test3=ggplot(dplrcombinedtreatment, aes(x=as.factor(Date))) +
  geom_line(aes(y=combestavg, group=dplrcombinedtreatment$Cage_treatment, color=dplrcombinedtreatment$Cage_treatment), size=1.3)+
  geom_errorbar(aes(ymin=combestavg-combestse,ymax=combestavg+combestse, color=dplrcombinedtreatment$Cage_treatment), size=.5)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "LB", "Apple"),
                      labels = c("+AT", "+LB", "Control"),
                      values = c("red","blue","darkgrey"))+
  geom_line( aes(y=MinAIRTC * coeff), linetype="dotted", size=1, group=1, color="black", alpha=.6) + 
  geom_line( aes(y=MaxAirTC * coeff), linetype="dashed",size=1, group=1, color="black", alpha=.6) +
  scale_x_discrete(
    name=" Experiment Date", 
    breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
    labels=c("Start","August 14th", "September 18th", "October 18th", "November 10th", "End"))+
  scale_y_continuous(
    # Features of the first axis
    name = "Population Size",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Min/ Max Daily Temperature (C°)"))+
  theme_cowplot()

test3+ theme(legend.position='bottom')

####linear regression analysis####
#
##All Timepoints##
CensusALL = census
names(census)[9] ="picavg"
census <- census[! census$Cage_treatment %in% c('Bloom'), ]

CensusALL$Date=as.Date(CensusALL$Date)

CensusALL=CensusALL %>% 
  mutate(TP = rank(Date)) 

CensusALL$TP=as.numeric(CensusALL$TP)
CensusALL$Date=as.factor(CensusALL$Date)
CensusALL$DateCage <- paste(CensusALL$Date,CensusALL$Cage_number)
View(CensusALL)

CensusALLanv<-lme(picavg ~ Date*Cage_treatment, random=~1|Cage_number/Date, data=CensusALL)
anova(CensusALLanv) ## sig Date, insig cage_treatment and interaction 


pairs=lsmeans(CensusALLanv, pairwise ~ Cage_treatment|Date)
summ=summary(pairs$contrasts)
summ=summ %>% pivot_wider(Date, names_from = contrast, values_from = c("t.ratio","p.value"))
View(summ)

CensusALL=CensusALL %>%
  mutate(
    treat = case_when(
      Cage_treatment == "Apple" ~ "C",
      Cage_treatment == "AT" ~ "MA",
      Cage_treatment == "LB" ~ "MA"
    )
  )

CensusALLTreatanv<-lme(picavg ~ Date*treat, random=~1|Cage_number/Date, data=CensusALL)
anova(CensusALLTreatanv) ## sig Date, insig cage_treatment and interaction 


Treatpairs=lsmeans(CensusALLTreatanv, pairwise ~ treat|Date)
Treatsumm=summary(Treatpairs$contrasts)
Treatsumm=Treatsumm %>% pivot_wider(Date, names_from = contrast, values_from = c("t.ratio","p.value"))
View(Treatsumm)
write.csv(Treatsumm, "Treat_MA_lsmeans-pairwise-fig1.csv")





