####Figure 2####
install.packages("patchwork")

library(plyr)
library(cowplot)
library(nlme)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
library(ggbreak)
library(patchwork)
library(ggpubr)

####Import####
LDM = read.csv("raw/O20MA_LDM_raw.csv")
LDF = read.csv("raw/O20MA_LDF_raw.csv")
StarvM = read.csv("raw/O20MA_StarvM_raw.csv")
StarvF = read.csv("raw/O20MA_StarvF_raw.csv")
Fecund = read.csv("raw/O20MA_Fecund_raw.csv")
Viab = read.csv("raw/O20MA_Viab_raw.csv")
BSLipid = read.csv("raw/O20MA_BSLipid_raw.csv")
Pig = read.csv("raw/O20MA_Pig_raw.csv")


####Make Summary dfs####

######LDM####
##by treatment / timepoint##
LDM$M.AVG.LD=as.numeric(LDM$M.AVG.LD)
LDM.sum=LDM %>%
  group_by(timepoint,cage.treat2) %>%
  reframe(
                   N=n(),
                   timepoint=as.factor(timepoint),
                   cage.treat=as.factor(cage.treat),
                   cage.treat2=as.factor(cage.treat2),
                   LDM=mean(M.AVG.LD,na.rm = TRUE),
                   LDMsd=sd(M.AVG.LD,na.rm = TRUE),
                   LDMse=LDMsd/sqrt(N),
  )
LDM.sum=unique(LDM.sum)
View(LDM.sum)

##by cage / timepoint / treatment##
LDM$M.AVG.LD=as.numeric(LDM$M.AVG.LD)

LDM.cage.sum=LDM %>%
  group_by(timepoint,cage) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    LDM=mean(M.AVG.LD,na.rm = TRUE),
    LDMsd=sd(M.AVG.LD,na.rm = TRUE),
    LDMse=LDMsd/sqrt(N),
  )
LDM.cage.sum=unique(LDM.cage.sum)
View(LDM.cage.sum)


######LDF####
##by treatment / timepoint##
LDF$F.AVG.LD=as.numeric(LDF$F.AVG.LD)
LDF.sum=LDF %>%
  group_by(timepoint,cage.treat2) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    LDF=mean(F.AVG.LD,na.rm = TRUE),
    LDFsd=sd(F.AVG.LD,na.rm = TRUE),
    LDFse=LDFsd/sqrt(N),
  )
LDF.sum=unique(LDF.sum)
View(LDF.sum)

##by cage / timepoint / treatment##
LDF$F.AVG.LD=as.numeric(LDF$F.AVG.LD)

LDF.cage.sum=LDF %>%
  group_by(timepoint,cage) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    LDF=mean(F.AVG.LD,na.rm = TRUE),
    LDFsd=sd(F.AVG.LD,na.rm = TRUE),
    LDFse=LDFsd/sqrt(N),
  )
LDF.cage.sum=unique(LDF.cage.sum)
View(LDF.cage.sum)

######SRM####
##by treatment / time point##
View(StarvM)
StarvM <- StarvM %>%
  mutate(mean.starv = (
                         (ifelse(`X20` - `X16` > 0, `X20` - `X16`, 0) * 20) + 
                         (ifelse(`X24` - `X20` > 0, `X24` - `X20`, 0) * 24) + 
                         (ifelse(`X40` - `X24` > 0, `X40` - `X24`, 0) * 40) + 
                         (ifelse(`X44` - `X40` > 0, `X44` - `X40`, 0) * 44) + 
                         (ifelse(`X48` - `X44` > 0, `X48` - `X44`, 0) * 48) + 
                         (ifelse(`X64` - `X48` > 0, `X64` - `X48`, 0) * 64) + 
                         (ifelse(`X68` - `X64` > 0, `X68` - `X64`, 0) * 68) + 
                         (ifelse(`X72` - `X68` > 0, `X72` - `X68`, 0) * 72) + 
                         (ifelse(`X88` - `X72` > 0, `X88` - `X72`, 0) * 88) + 
                         (ifelse(`X92` - `X88` > 0, `X92` - `X88`, 0) * 92) + 
                         (ifelse(`X96` - `X92` > 0, `X96` - `X92`, 0) * 96) + 
                         (ifelse(`X112` - `X96` > 0, `X112` - `X96`, 0) * 112) + 
                         (ifelse(`X116` - `X112` > 0, `X116` - `X112`, 0) * 116) + 
                         (ifelse(`X120` - `X116` > 0, `X120` - `X116`, 0) * 120) + 
                         (ifelse(`X136` - `X120` > 0, `X136` - `X120`, 0) * 136) + 
                         (ifelse(`X140` - `X136` > 0, `X140` - `X136`, 0) * 140) + 
                         (ifelse(`X144` - `X140` > 0, `X144` - `X140`, 0) * 144) + 
                         (ifelse(`X160` - `X144` > 0, `X160` - `X144`, 0) * 160) + 
                         (ifelse(`X164` - `X160` > 0, `X164` - `X160`, 0) * 164) + 
                         (ifelse(`X168` - `X164` > 0, `X168` - `X164`, 0) * 168) +
                         (ifelse(`X184` - `X168` > 0, `X184` - `X168`, 0) * 184) + 
                         (ifelse(`X188` - `X184` > 0, `X188` - `X184`, 0) * 188) + 
                         (ifelse(`X192` - `X188` > 0, `X192` - `X188`, 0) * 192)
                       
                       
                       ) / `X192`)

SRM.sum=StarvM %>%
  group_by(timepoint,cage.treat2) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    SRM=mean(mean.starv,na.rm = TRUE),
    SRMsd=sd(mean.starv,na.rm = TRUE),
    SRMse=SRMsd/sqrt(N),
  )
SRM.sum=unique(SRM.sum)
#View(SRM.sum)

##by cage / timepoint / treatment##
SRM.cage.sum=StarvM %>%
  group_by(timepoint, cage) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    SRM=mean(mean.starv,na.rm = TRUE),
    SRMsd=sd(mean.starv,na.rm = TRUE),
    SRMse=SRMsd/sqrt(N),
  )
SRM.cage.sum=unique(SRM.cage.sum)
View(SRM.cage.sum)

######SRF####
##by treatment / time point##
View(StarvF)

StarvF <- StarvF %>%
  mutate(mean.starv = (
    (ifelse(`X20` - `X16` > 0, `X20` - `X16`, 0) * 20) + 
      (ifelse(`X24` - `X20` > 0, `X24` - `X20`, 0) * 24) + 
      (ifelse(`X40` - `X24` > 0, `X40` - `X24`, 0) * 40) + 
      (ifelse(`X44` - `X40` > 0, `X44` - `X40`, 0) * 44) + 
      (ifelse(`X48` - `X44` > 0, `X48` - `X44`, 0) * 48) + 
      (ifelse(`X64` - `X48` > 0, `X64` - `X48`, 0) * 64) + 
      (ifelse(`X68` - `X64` > 0, `X68` - `X64`, 0) * 68) + 
      (ifelse(`X72` - `X68` > 0, `X72` - `X68`, 0) * 72) + 
      (ifelse(`X88` - `X72` > 0, `X88` - `X72`, 0) * 88) + 
      (ifelse(`X92` - `X88` > 0, `X92` - `X88`, 0) * 92) + 
      (ifelse(`X96` - `X92` > 0, `X96` - `X92`, 0) * 96) + 
      (ifelse(`X112` - `X96` > 0, `X112` - `X96`, 0) * 112) + 
      (ifelse(`X116` - `X112` > 0, `X116` - `X112`, 0) * 116) + 
      (ifelse(`X120` - `X116` > 0, `X120` - `X116`, 0) * 120) + 
      (ifelse(`X136` - `X120` > 0, `X136` - `X120`, 0) * 136) + 
      (ifelse(`X140` - `X136` > 0, `X140` - `X136`, 0) * 140) + 
      (ifelse(`X144` - `X140` > 0, `X144` - `X140`, 0) * 144) + 
      (ifelse(`X160` - `X144` > 0, `X160` - `X144`, 0) * 160) + 
      (ifelse(`X164` - `X160` > 0, `X164` - `X160`, 0) * 164) + 
      (ifelse(`X168` - `X164` > 0, `X168` - `X164`, 0) * 168) +
      (ifelse(`X184` - `X168` > 0, `X184` - `X168`, 0) * 184) + 
      (ifelse(`X188` - `X184` > 0, `X188` - `X184`, 0) * 188) + 
      (ifelse(`X192` - `X188` > 0, `X192` - `X188`, 0) * 192)
    
    
  ) / `X192`)

SRF.sum=StarvF %>%
  group_by(timepoint,cage.treat2) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    SRF=mean(mean.starv,na.rm = TRUE),
    SRFsd=sd(mean.starv,na.rm = TRUE),
    SRFse=SRFsd/sqrt(N),
  )
SRF.sum=unique(SRF.sum)
View(SRF.sum)

##by cage / timepoint / treatment##
SRF.cage.sum=StarvF %>%
  group_by(timepoint, cage) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    SRF=mean(mean.starv,na.rm = TRUE),
    SRFsd=sd(mean.starv,na.rm = TRUE),
    SRFse=SRFsd/sqrt(N),
  )
SRF.cage.sum=unique(SRF.cage.sum)
View(SRF.cage.sum)

######Fecundity####
View(Fecund)

Fec.sum=Fecund %>%
  group_by(timepoint,cage.treat2) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    Fec=mean(totalfecund,na.rm = TRUE),
    Fecsd=sd(totalfecund,na.rm = TRUE),
    Fecse=Fecsd/sqrt(N),
  )
Fec.sum=unique(Fec.sum)

View(Fec.sum)

##by cage/ timepoint/ treatment##

Fec.cage.sum=Fecund %>%
  group_by(timepoint, cage) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    Fec=mean(totalfecund,na.rm = TRUE),
    Fecsd=sd(totalfecund,na.rm = TRUE),
    Fecse=Fecsd/sqrt(N),
  )
Fec.cage.sum=unique(Fec.cage.sum)

#View(Fec.cage.sum)

######Viability####

View(Viab)
Via.sum=Viab %>%
  group_by(timepoint,cage.treat2) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    Viab=mean(viability,na.rm = TRUE),
    Viabsd=sd(viability,na.rm = TRUE),
    Viabse=Viabsd/sqrt(N),
  )
Via.sum=unique(Via.sum)
View(Via.sum)

##by cage/ timepoint/ treatment##

Via.cage.sum=Viab %>%
  group_by(timepoint,cage) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    Viab=mean(viability,na.rm = TRUE),
    Viabsd=sd(viability,na.rm = TRUE),
    Viabse=Viabsd/sqrt(N),
  )
Via.cage.sum=unique(Via.cage.sum)
View(Via.cage.sum)

######BS LIPID######

View(BSLipid)
BS.sum=BSLipid %>%
  group_by(timepoint,cage.treat2) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    DW=mean(dry.weight,na.rm = TRUE),
    DWsd=sd(dry.weight,na.rm = TRUE),
    DWse=DWsd/sqrt(N),
    LW=mean(lipid.weight,na.rm = TRUE),
    LWsd=sd(lipid.weight,na.rm = TRUE),
    LWse=LWsd/sqrt(N),
    PercLipid=(LW/DW)/100,
    PercLipidse=(LWse/DWse)/10000,
  )
BS.sum=unique(BS.sum)

View(BS.sum)

##by cage/ timepoint/ treatment##

BS.cage.sum=BSLipid %>%
  group_by(timepoint,Cage) %>%
  reframe(
    N=n(),
    timepoint=as.factor(timepoint),
    cage=as.factor(Cage),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    DW=mean(dry.weight,na.rm = TRUE),
    DWsd=sd(dry.weight,na.rm = TRUE),
    DWse=DWsd/sqrt(N),
    LW=mean(lipid.weight,na.rm = TRUE),
    LWsd=sd(lipid.weight,na.rm = TRUE),
    LWse=LWsd/sqrt(N),
    PercLipid=(LW/DW)/100,
    PercLipidse=(LWse/DWse)/10000,
  )
BS.cage.sum=unique(BS.cage.sum)
View(BS.cage.sum)


######PIG####
View(Pig)
Pig.sum=Pig %>%
  group_by(Time.Point,cage.treat2) %>%
  reframe(
    N=n(),
    timepoint=as.factor(Time.Point),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    Pig=mean(AVG,na.rm = TRUE),
    Pigsd=sd(AVG,na.rm = TRUE),
    Pigse=Pigsd/sqrt(N),
  )
Pig.sum=unique(Pig.sum)

View(Pig.sum)

##by cage/ timepoint/ treatment##

Pig.cage.sum=Pig %>%
  group_by(Time.Point,Cage) %>%
  reframe(
    N=n(),
    timepoint=as.factor(Time.Point),
    cage=as.factor(Cage),
    cage.treat=as.factor(cage.treat),
    cage.treat2=as.factor(cage.treat2),
    Pig=mean(AVG,na.rm = TRUE),
    Pigsd=sd(AVG,na.rm = TRUE),
    Pigse=Pig/100,
  )
Pig.cage.sum=unique(Pig.cage.sum)

View(Pig.cage.sum)



####Merging Dfs####
##by tp/treatment##

Phenos.list <- list(LDM.sum, LDF.sum, SRM.sum, SRF.sum, Fec.sum, Via.sum, BS.sum, Pig.sum)      
#View(Phenos.list)
#merge all data frames together
MAPhenos=Phenos.list %>% reduce(full_join, by=c('timepoint','cage.treat','cage.treat2'))
#View(MAPhenos)
##Remove N columns
MAPhenos <- MAPhenos[, -c(3,8,12,16,20,24,28,37,38)]
#View(MAPhenos)

write.csv (MAPhenos, "MAPhenos.csv", row.names=FALSE)

#View(Phenos) #### All phenotypic data in one table by treatment and time point

##by tp/cage/treatment##

Phenos.cage.list <- list(LDM.cage.sum, LDF.cage.sum, SRM.cage.sum, SRF.cage.sum, Fec.cage.sum, Via.cage.sum, BS.cage.sum, Pig.cage.sum)      

#merge all data frames together
MAPhenos.cage=Phenos.cage.list %>% reduce(full_join, by=c('cage','timepoint','cage.treat','cage.treat2'))
#View(MAPhenos.cage)

##Remove N columns
MAPhenos.cage <- MAPhenos.cage[, -c(3,9,13,17,21,25,29,30,39,40,41,43)]

#View(MAPhenos.cage) #### All phenotypic data in one table by cage, treatment and time point
write.csv (MAPhenos.cage, "MAPhenos.cage.csv", row.names=FALSE)

##Subsetting##
MAfounder= subset(MAPhenos, MAPhenos$timepoint == 0)
MAfounder.cage=subset(MAPhenos.cage, MAPhenos.cage$timepoint == 0)

MAPhenos.24.cage=subset(MAPhenos.cage, MAPhenos.cage$timepoint %in% c("2","4"))
View(MAPhenos.24.cage)

MAphenos=subset(MAPhenos, MAPhenos$timepoint %in% c("2","4","5"))
View(MAphenos)
MAphenos.cage=subset(MAPhenos.cage, MAPhenos.cage$timepoint %in% c("2","4","5"))
View(MAphenos.cage)
####PCA####
MApheno.cage.pca<-prcomp(MAPhenos.cage[c(5,8,11,14,17,20,31)],scale=TRUE)
MApheno.cage.pca.df<-data.frame(cage=MAPhenos.cage$cage,timepoint=MAPhenos.cage$timepoint,treat=MAPhenos.cage$cage.treat,treat2=MAPhenos.cage$cage.treat2,MApheno.cage.pca$x)
View(MApheno.cage.pca.df)

MApheno.24.cage.pca<-prcomp(MAPhenos.24.cage[c(5,8,11,14,17,20,23,26,31)],scale=TRUE)
MApheno.24.cage.pca.df<-data.frame(cage=MAPhenos.24.cage$cage,timepoint=MAPhenos.24.cage$timepoint,treat=MAPhenos.24.cage$cage.treat,treat2=MAPhenos.24.cage$cage.treat2,MApheno.24.cage.pca$x)
View(MApheno.24.cage.pca.df)

#### viz and manova
MAPheno.pca=ggplot(MApheno.cage.pca.df, aes(x=PC1, y=PC2, color = treat2, shape = timepoint, group=treat)) +
  geom_point(size=5)+
  stat_ellipse()+
  scale_colour_manual( 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey","blue"))+
  theme_bw()+
  theme(
    legend.position = c(.2, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(.5, .5, .5,.5), 
    legend.text=element_text(size=12),
    legend.title = element_blank()
  )
MAPheno.pca

MApheno.cage.pca.df.manov<-manova(cbind(PC1,PC2)~timepoint*treat,MApheno.cage.pca.df)
summary(MApheno.cage.pca.df.manov)
MApheno2.cage.pca.df.manov<-manova(cbind(PC1,PC2)~timepoint*treat2,MApheno.cage.pca.df)
summary(MApheno2.cage.pca.df.manov) ## C vs, AT vs LB Time and interaction sig 


MAPheno.24.pca=ggplot(MApheno.24.cage.pca.df, aes(x=PC1, y=PC2, color = treat2, shape = timepoint, group=treat)) +
  geom_point(size=5)+
  stat_ellipse()+
  scale_colour_manual( 
    breaks = c("AT", "C", "LB"),
    labels = c("+AT","Control", "+LB"),
    values = c("red","grey","blue"))+
  theme_bw()+
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(.5, .5, .5,.5), 
    legend.text=element_text(size=12),
    legend.title = element_blank()
  )
MAPheno.24.pca

MApheno.24.cage.pca.df.manov<-manova(cbind(PC1,PC2)~timepoint*treat,MApheno.24.cage.pca.df)
summary(MApheno.24.cage.pca.df.manov) ## C vs. MA VERY SIG TREAT, TIME and Interaction
MApheno2.24.cage.pca.df.manov<-manova(cbind(PC1,PC2)~timepoint*treat2,MApheno.24.cage.pca.df)
summary(MApheno2.24.cage.pca.df.manov) ## C Vs. AT VS. LB SIG TREAT AND TIME NO interaction 

Fig2a=ggarrange(MAPheno.pca,MAPheno.24.pca,
                labels = c("J", "K"),
                ncol = 2, nrow = 1, common.legend = FALSE)
Fig2a
####GRAPHING Fig 2 ####

## Male Development Time Fig##
#View(MAfounder)
#View(MAphenos)
LDM1 <- ggplot(MAphenos, aes(x=as.numeric(timepoint), y=LDM, group=cage.treat2))+
  geom_point(data = MAfounder, aes(x=as.numeric(timepoint), y=LDM,color = cage.treat2), size=5)+
  geom_linerange(data=MAfounder, aes(ymin=LDM-LDMse,ymax=LDM+LDMse, color= cage.treat2), size=1)+
  geom_point(aes(color=cage.treat2),size=5)+
  geom_line(aes(color=cage.treat2))+
  geom_linerange(aes(ymin=LDM-LDMse,ymax=LDM+LDMse, color= cage.treat2), size=1)+
  ylab("Male Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_continuous(limits = c(.9,4.5), breaks = c(1,2,3,4),labels=c("F","TP2","TP4","TP5"))+
  scale_x_break(c(1.2,2.0), scales = 13)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey", "blue"))+
  geom_line(data=MAphenos.cage, aes(x=as.numeric(timepoint), y=LDM, group=as.factor(cage), color=cage.treat2), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
LDM1

## Female Development Time Fig##

LDF1 <- ggplot(MAphenos, aes(x=as.numeric(timepoint), y=LDF, group=cage.treat2))+
  geom_point(data = MAfounder, aes(x=as.numeric(timepoint), y=LDF,color = cage.treat2), size=5)+
  geom_linerange(data=MAfounder, aes(ymin=LDF-LDFse,ymax=LDF+LDFse, color= cage.treat2), size=1)+
  geom_point(aes(color=cage.treat2),size=5)+
  geom_line(aes(color=cage.treat2))+
  geom_linerange(aes(ymin=LDF-LDFse,ymax=LDF+LDFse, color= cage.treat2), size=1)+
  ylab("Female Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_continuous(limits = c(.9,4.5), breaks = c(1,2,3,4),labels=c("F","TP2","TP4","TP5"))+
  scale_x_break(c(1.2,2.0), scales = 13)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey", "blue"))+
  geom_line(data=MAphenos.cage, aes(x=as.numeric(timepoint), y=LDF, group=as.factor(cage), color=cage.treat2), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
LDF1

##Starvation Resistance Female##

SRF1 <- ggplot(MAphenos, aes(x=as.numeric(timepoint), y=SRF, group=cage.treat2))+
  geom_point(data = MAfounder, aes(x=as.numeric(timepoint), y=SRF,color = cage.treat2), size=5)+
  geom_linerange(data=MAfounder, aes(ymin=SRF-SRFse,ymax=SRF+SRFse, color= cage.treat2), size=1)+
  geom_point(aes(color=cage.treat2),size=5)+
  geom_line(aes(color=cage.treat2))+
  geom_linerange(aes(ymin=SRF-SRFse,ymax=SRF+SRFse, color= cage.treat2), size=1)+
  ylab("Female Starvation Resistance (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_continuous(limits = c(.9,4.5), breaks = c(1,2,3,4),labels=c("F","TP2","TP4","TP5"))+
  scale_x_break(c(1.2,2.0), scales = 13)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey", "blue"))+
  geom_line(data=MAphenos.cage, aes(x=as.numeric(timepoint), y=SRF, group=as.factor(cage), color=cage.treat2), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
SRF1

##Starvation Resistance Male##

SRM1 <- ggplot(MAphenos, aes(x=as.numeric(timepoint), y=SRM, group=cage.treat2))+
  geom_point(data = MAfounder, aes(x=as.numeric(timepoint), y=SRM,color = cage.treat2), size=5)+
  geom_linerange(data=MAfounder, aes(ymin=SRM-SRMse,ymax=SRM+SRMse, color= cage.treat2), size=1)+
  geom_point(aes(color=cage.treat2),size=5)+
  geom_line(aes(color=cage.treat2))+
  geom_linerange(aes(ymin=SRM-SRMse,ymax=SRM+SRMse, color= cage.treat2), size=1)+
  ylab("Male Starvation Resistance (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_continuous(limits = c(.9,4.5), breaks = c(1,2,3,4),labels=c("F","TP2","TP4","TP5"))+
  scale_x_break(c(1.2,2.0), scales = 13)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey", "blue"))+
  geom_line(data=MAphenos.cage, aes(x=as.numeric(timepoint), y=SRM, group=as.factor(cage), color=cage.treat2), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
SRM1

##Fecundity per day##

Fec1 <- ggplot(MAphenos, aes(x=as.numeric(timepoint), y=Fec, group=cage.treat2))+
  geom_point(data = MAfounder, aes(x=as.numeric(timepoint), y=Fec,color = cage.treat2), size=5)+
  geom_linerange(data=MAfounder, aes(ymin=Fec-Fecse,ymax=Fec+Fecse, color= cage.treat2), size=1)+
  geom_point(aes(color=cage.treat2),size=5)+
  geom_line(aes(color=cage.treat2))+
  geom_linerange(aes(ymin=Fec-Fecse,ymax=Fec+Fecse, color= cage.treat2), size=1)+
  ylab("3 Day Total Fecundity ")+
  xlab("Sample Timepoints")+
  scale_x_continuous(limits = c(.9,4.5), breaks = c(1,2,3,4),labels=c("F","TP2","TP4","TP5"))+
  scale_x_break(c(1.2,2.0), scales = 13)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey", "blue"))+
  geom_line(data=MAphenos.cage, aes(x=as.numeric(timepoint), y=Fec, group=as.factor(cage), color=cage.treat2), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
Fec1


##Viability ##

Viab1 <- ggplot(MAphenos, aes(x=as.numeric(timepoint), y=Viab, group=cage.treat2))+
  geom_point(data = MAfounder, aes(x=as.numeric(timepoint), y=Viab,color = cage.treat2), size=5)+
  geom_linerange(data=MAfounder, aes(ymin=Viab-Viabse,ymax=Viab+Viabse, color= cage.treat2), size=1)+
  geom_point(aes(color=cage.treat2),size=5)+
  geom_line(aes(color=cage.treat2))+
  geom_linerange(aes(ymin=Viab-Viabse,ymax=Viab+Viabse, color= cage.treat2), size=1)+
  ylab("Egg to Adult Viability")+
  xlab("Sample Timepoints")+
  scale_x_continuous(limits = c(.9,4.5), breaks = c(1,2,3,4),labels=c("F","TP2","TP4","TP5"))+
  scale_x_break(c(1.2,2.0), scales = 13)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey", "blue"))+
  geom_line(data=MAphenos.cage, aes(x=as.numeric(timepoint), y=Viab, group=as.factor(cage), color=cage.treat2), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
Viab1

##Dry Weight ##

DW1 <- ggplot(MAphenos, aes(x=as.numeric(timepoint), y=DW, group=cage.treat2))+
  geom_point(data = MAfounder, aes(x=as.numeric(timepoint), y=DW,color = cage.treat2), size=5)+
  geom_linerange(data=MAfounder, aes(ymin=DW-DWse,ymax=DW+DWse, color= cage.treat2), size=1)+
  geom_point(aes(color=cage.treat2),size=5)+
  geom_line(aes(color=cage.treat2))+
  geom_linerange(aes(ymin=DW-DWse,ymax=DW+DWse, color= cage.treat2), size=1)+
  ylab("Pooled Dry Weight (g)")+
  xlab("Sample Timepoints")+
  scale_x_continuous(limits = c(.9,3.5), breaks = c(1,2,3,4),labels=c("F","TP2","TP4","TP5"))+
  scale_x_break(c(1.2,2.0), scales = 13)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey", "blue"))+
  geom_line(data=MAphenos.cage, aes(x=as.numeric(timepoint), y=DW, group=as.factor(cage), color=cage.treat2), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
DW1

##Lipid Weight ##

LW1 <- ggplot(MAphenos, aes(x=as.numeric(timepoint), y=LW, group=cage.treat2))+
  geom_point(data = MAfounder, aes(x=as.numeric(timepoint), y=LW,color = cage.treat2), size=5)+
  geom_linerange(data=MAfounder, aes(ymin=LW-LWse,ymax=LW+LWse, color= cage.treat2), size=1)+
  geom_point(aes(color=cage.treat2),size=5)+
  geom_line(aes(color=cage.treat2))+
  geom_linerange(aes(ymin=LW-LWse,ymax=LW+LWse, color= cage.treat2), size=1)+
  ylab("Pooled Lipid Weight (g)")+
  xlab("Sample Timepoints")+
  scale_x_continuous(limits = c(.9,3.5), breaks = c(1,2,3,4),labels=c("F","TP2","TP4","TP5"))+
  scale_x_break(c(1.2,2.0), scales = 13)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey", "blue"))+
  geom_line(data=MAphenos.cage, aes(x=as.numeric(timepoint), y=LW, group=as.factor(cage), color=cage.treat2), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
LW1



##Pigmentation##

Pig1 <- ggplot(MAphenos, aes(x=as.numeric(timepoint), y=Pig, group=cage.treat2))+
  geom_point(data = MAfounder, aes(x=as.numeric(timepoint), y=Pig,color = cage.treat2), size=5)+
  geom_linerange(data=MAfounder, aes(ymin=Pig-Pigse,ymax=Pig+Pigse, color= cage.treat2), size=1)+
  geom_point(aes(color=cage.treat2),size=5)+
  geom_line(aes(color=cage.treat2))+
  geom_linerange(aes(ymin=Pig-Pigse,ymax=Pig+Pigse, color= cage.treat2), size=1)+
  ylab("Percent Melanization")+
  xlab("Sample Timepoints")+
  scale_x_continuous(limits = c(.9,4.5), breaks = c(1,2,3,4),labels=c("F","TP2","TP4","TP5"))+
  scale_x_break(c(1.2,2.0), scales = 13)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AT", "C", "LB"),
                      labels = c("+AT","Control", "+LB"),
                      values = c("red","grey", "blue"))+
  geom_line(data=MAphenos.cage, aes(x=as.numeric(timepoint), y=Pig, group=as.factor(cage), color=cage.treat2), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
Pig1

######Combine####

library(ggpubr)



MAFig2 <- ggarrange(LDM1+ rremove("xlab"), LDF1+ rremove("xlab"), DW1+ rremove("xlab"), SRM1 + rremove("xlab"), SRF1 + rremove("xlab"), LW1 + rremove("xlab"), Viab1, Pig1, Fec1,
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                  ncol = 3, nrow = 3, common.legend = TRUE)
MAFig2
save(Fig2, file = "fig2.rdata")



####


######Linear Mixed Model Analysis F####
###Found in Table 2###

 ### TP 2,4,5
MALmmdata=MAphenos.cage
MALmmdata$timepoint=as.factor(MALmmdata$timepoint)
MALmmdata$cage.treat=as.factor(MALmmdata$cage.treat)
MALmmdata$cage.treat2=as.factor(MALmmdata$cage.treat2)
View(MALmmdata)
  
  
MALDm1 <- lmer(LDM ~ timepoint*cage.treat + (1|cage), data = MALmmdata)  
anova(MALDm1)
lsmeans(MALDm1, pairwise ~ cage.treat)

MALDm2 <- lmer(LDM ~ timepoint*cage.treat2 + (1|cage), data = MALmmdata)  
anova(MALDm2)
lsmeans(MALDm2, pairwise ~ cage.treat2)

MALDf1 <- lmer(LDF ~ timepoint*cage.treat + (1|cage), data = MALmmdata)  
anova(MALDf1)
lsmeans(MALDf1, pairwise ~ cage.treat)

MALDf2 <- lmer(LDF ~ timepoint*cage.treat2 + (1|cage), data = MALmmdata)  
anova(MALDf2)
lsmeans(MALDf2, pairwise ~ cage.treat2)

MASRm1 <- lmer(SRM ~ timepoint*cage.treat + (1|cage), data = MALmmdata)  
anova(MASRm1)
lsmeans(MASRm1, pairwise ~ cage.treat)

MASRm2 <- lmer(SRM ~ timepoint*cage.treat2 + (1|cage) , data = MALmmdata)  
anova(MASRm2)
lsmeans(MASRm2, pairwise ~ cage.treat2)

MASRf1 <- lmer(SRF ~ timepoint*cage.treat + (1|cage), data = MALmmdata)  
anova(MASRf1)
lsmeans(MASRf1, pairwise ~ cage.treat)

MASRf2 <- lmer(SRF ~ timepoint*cage.treat2 + (1|cage) , data = MALmmdata)  
anova(MASRf2)
lsmeans(MASRf2, pairwise ~ cage.treat2)

MAFec1 <- lmer(Fec ~ timepoint*cage.treat + (1|cage), data = MALmmdata)  
anova(MAFec1)
lsmeans(MAFec1, pairwise ~ cage.treat)

MAFec2 <- lmer(Fec ~ timepoint*cage.treat2 + (1|cage), data = MALmmdata)  
anova(MAFec2)
lsmeans(MAFec2, pairwise ~ cage.treat2)

MAV1 <- lmer(Viab ~ timepoint*cage.treat + (1|cage), data = MALmmdata)  
anova(MAV1)
lsmeans(MAV1, pairwise ~ cage.treat)

MAV2<- lmer(Viab ~ timepoint*cage.treat2 + (1|cage) , data = MALmmdata)  
anova(MAV2)
lsmeans(MAV2, pairwise ~ cage.treat2)

MALW1 <- lmer(LW ~ timepoint*cage.treat + (1|cage), data = MALmmdata)  
anova(MALW1)
lsmeans(MALW1, pairwise ~ cage.treat)

MALW2<- lmer(LW ~ timepoint*cage.treat2 + (1|cage) , data = MALmmdata)  
anova(MALW2)
lsmeans(MALW2, pairwise ~ cage.treat2)

MADW1 <- lmer(DW ~ timepoint*cage.treat + (1|cage), data = MALmmdata)  
anova(MADW1)
lsmeans(MADW1, pairwise ~ cage.treat)

MADW2<- lmer(DW ~ timepoint*cage.treat2 + (1|cage) , data = MALmmdata)  
anova(MADW2)
lsmeans(MADW2, pairwise ~ cage.treat2)

MAPIG1 <- lmer(Pig ~ timepoint*cage.treat + (1|cage), data = MALmmdata)  
anova(MAPIG1)
lsmeans(MAPIG1, pairwise ~ cage.treat)

MAPIG2<- lmer(Pig ~ timepoint*cage.treat2 + (1|cage) , data = MALmmdata)  
anova(MAPIG2)
lsmeans(MAPIG2, pairwise ~ cage.treat2)



######PHENO PCA####
library(reshape2)
library(plyr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(latex2exp)
library(ggrepel)
library(scales)
library(grid)
library(ggtext)
library(ggpubr)

# PREPARATION

MApca245 <- MAPhenos.cage ### this df is produced when generating figure 2, must complete that first. 
MApca245 <- MApca245[, -c(6,7,9,10,12,13,15,16,18,19,21,22,23:30,32)]## remove SE and SD
View(MApca245)

MApca24 <- MAPhenos.24.cage ### this df is produced when generating figure 2, must complete that first. 
MApca24 <- MApca24[, -c(6,7,9,10,12,13,15,16,18,19,21,22,24,25,27:30,32)]## remove SE and SD
View(MApca24)

#View(ABpca)


#### LISTS AND COLOR ASSIGNMENTS####



MA2ColorListTime <- c("2"="grey",
                      "4"="yellow",
                      "5"= "blue")



MA2ColorListPopulation<-c("C"="grey", "AT" = "red","LB" ="blue",
                          "MA"="purple")

#View(AB2ColorListPopulation)

MA2ColorListMDS=c("C"="grey", "AT" = "red","LB" ="blue",
                  "MA"="purple")

MA2ColorListTrait<-c(
  "LDf"="black",
  "LDf"="black",
  "Viab"="black",
  "SRM"="black",
  "SRF"="black",
  "TOTF"="black",
  "BS"="black")


#display.brewer.pal(8,"Dark2")

MA2TimeList=c("2","4","5")
MA2PopulationList=c("C","AT","LB","MA")
MA24TraitsList<-data.frame(ABBRV=c("LDm","LDf","Viab","SRM","SRF","Fec","Pig","LW","DW" ),LONG=c("Development Time (Males)","Development Time (Females)", "Viability","Starvation Resistance (Males)","Starvation Resistance (Females)","Fecundity","Pigmentation", "Lipid Weight", "Dry Weight"))
MA245TraitsList<-data.frame(ABBRV=c("LDm","LDf","Viab","SRM","SRF","Fec","Pig"),LONG=c("Development Time (Males)","Development Time (Females)", "Viability","Starvation Resistance (Males)","Starvation Resistance (Females)","Fecundity","Pigmentation"))

MA2SelectionLabels=c("TP2","TP4","TP5")
MA2PopulationLabels=c("Control","+AT","+LB","Microbial Addition" )
MA2TimeShape=c("2"=8,"4"=17,"5"=19)



#SETTING A CONSISTENT IN ALL FIGURES THEME

{
  {
    BASESIZE=8
    TITLESIZE=7
    FONTFAMILY="sans"
    TEXTCOLOR="black"
    AXISTEXTSIZE=7
    AXISTEXTCOLOR="grey40"
    AXISTITLESIZE=7
    LEGENDTEXTSIZE=6
    LABELSIZE=11
  }
  
  theme_evo <- theme_set(theme_bw(base_size=BASESIZE,
                                  base_family=FONTFAMILY))
  theme_evo <- theme_update(panel.background = element_blank(),
                            title = element_text(size=TITLESIZE, face= "bold"),
                            axis.text = element_text(size= AXISTEXTSIZE, color = AXISTEXTCOLOR, family = FONTFAMILY, face= "bold"),
                            axis.title = element_text(size= TITLESIZE, color = TEXTCOLOR, family = FONTFAMILY, face= "bold"),
                            legend.background = element_blank(),
                            legend.box.background = element_blank(),
                            legend.title=element_text(face= "bold",size = LEGENDTEXTSIZE),
                            legend.text=element_text(size = LEGENDTEXTSIZE),
                            legend.key =element_blank()
  )
  
}

# EMPTY PLOT
pBL<-ggplot()+theme_nothing()


# PHENOTYPIC PCA - FIGURE 3
{
  
  #### PREP####
  {
    # FUNTION TO FIND HULLS
    find_hull<-function(df) df[chull(df$PC1,df$PC2),]
    
    # READING THE DATA
    # View(ABpca2)
    
    
    # SORT THE DATA
    MApca24$Time<-factor(MApca24$timepoint,MA2TimeList)
    MApca24$Population=factor(MApca24$cage.treat,MA2PopulationList)
    MApca24$Population2=factor(MApca24$cage.treat2,MA2PopulationList)
    MApca24<-MApca24 %>% arrange(cage,Time,Population)
    
    MApca245$Time<-factor(MApca245$timepoint,MA2TimeList)
    MApca245$Population=factor(MApca245$cage.treat,MA2PopulationList)
    MApca245$Population2=factor(MApca245$cage.treat2,MA2PopulationList)
    MApca245<-MApca245 %>% arrange(cage,Time,Population)
    
    View(MApca245)
    #### CALCULATING THE PRINCIPAL COMPONENTS####
    MA2_C24<- MApca24 %>% subset(Population=="C")
    MA2_AT24<- MApca24 %>% subset(Population2=="AT")
    MA2_LB24<- MApca24 %>% subset(Population2=="LB")
    MA2_MA24<- MApca24 %>% subset(Population=="MA")
    
    MA2_C245<- MApca245 %>% subset(Population=="C")
    MA2_AT245<- MApca245 %>% subset(Population2=="AT")
    MA2_LB245<- MApca245 %>% subset(Population2=="LB")
    MA2_MA245<- MApca245 %>% subset(Population=="MA")
    
    MA24_2<- MApca24 %>% subset(Time=="2")
    MA24_4<- MApca24 %>% subset(Time=="4")
    
    MA245_2<- MApca245 %>% subset(Time=="2")
    MA245_4<- MApca245 %>% subset(Time=="4")
    MA245_5<- MApca245 %>% subset(Time=="5")
  

    MA2_C24.pca<-prcomp(MA2_C24[5:13],scale=TRUE)
    MA2_AT24.pca<-prcomp(MA2_AT24[5:13],scale=TRUE)
    MA2_LB24.pca<-prcomp(MA2_LB24[5:13],scale=TRUE)
    MA2_MA24.pca<-prcomp(MA2_MA24[5:13],scale=TRUE)
    
    MA2_C245.pca<-prcomp(MA2_C245[5:11],scale=TRUE)
    MA2_AT245.pca<-prcomp(MA2_AT245[5:11],scale=TRUE)
    MA2_LB245.pca<-prcomp(MA2_LB245[5:11],scale=TRUE)
    MA2_MA245.pca<-prcomp(MA2_MA245[5:11],scale=TRUE)
    
    MA24.pca<-prcomp(MApca24[5:13],scale=TRUE)
    MA245.pca<-prcomp(MApca245[5:11],scale=TRUE)
    
    # GETTING PCA COORDINATES TO A DATA FRAME
    ##All Vars
    MA2C.24.pca.df<-data.frame(cage=MA2_C24$cage,Population=MA2_C24$Population,Population2=MA2_C24$Population2,Time=MA2_C24$Time,MA2_C24.pca$x)
    MA2AT.24.pca.df<-data.frame(cage=MA2_AT24$cage,Population=MA2_AT24$Population,Population2=MA2_AT24$Population2,Time=MA2_AT24$Time,MA2_AT24.pca$x)
    MA2LB.24.pca.df<-data.frame(cage=MA2_LB24$cage,Population=MA2_LB24$Population,Population2=MA2_LB24$Population2,Time=MA2_LB24$Time,MA2_LB24.pca$x)
    MA2MA.24.pca.df<-data.frame(cage=MA2_MA24$cage,Population=MA2_MA24$Population,Population2=MA2_MA24$Population2,Time=MA2_MA24$Time,MA2_MA24.pca$x)
    
    MA2C.245.pca.df<-data.frame(cage=MA2_C245$cage,Population=MA2_C245$Population,Population2=MA2_C245$Population2,Time=MA2_C245$Time,MA2_C245.pca$x)
    MA2AT.245.pca.df<-data.frame(cage=MA2_AT245$cage,Population=MA2_AT245$Population,Population2=MA2_AT245$Population2,Time=MA2_AT245$Time,MA2_AT245.pca$x)
    MA2LB.245.pca.df<-data.frame(cage=MA2_LB245$cage,Population=MA2_LB245$Population,Population2=MA2_LB245$Population2,Time=MA2_LB245$Time,MA2_LB245.pca$x)
    MA2MA.245.pca.df<-data.frame(cage=MA2_MA245$cage,Population=MA2_MA245$Population,Population2=MA2_MA245$Population2,Time=MA2_MA245$Time,MA2_MA245.pca$x)
    
    MA2.24.pca.df<-data.frame(cage=MApca24$cage,Population=MApca24$Population,Population2=MApca24$Population2,Time=MApca24$Time,MA24.pca$x)
    MA2.245.pca.df<-data.frame(cage=MApca245$cage,Population=MApca245$Population,Population2=MApca245$Population2,Time=MApca245$Time,MA245.pca$x)
    
    View(MA2.24.pca.df)
    
    
    ####PCA PHENOTYPES COMBINED####
    MAPheno24.pca1=ggplot(MA2.24.pca.df, aes(x=PC1, y=PC2, shape = Time, color = Population2, group=Population)) +
      stat_ellipse()+
      geom_point(size=3)+
      #scale_colour_manual(name = "Assay Diet", 
                          #breaks = c("a", "b"),
                         # labels = c("Low Quality","High Quality"),
                          #values = c("grey","black"))+
     # scale_shape_manual(breaks=c("0", "2", "4"), 
                         #label=c("Founder","Timepoint 2 (Summer)", "Timepoint 4 (Fall)"), 
                        # values =c(21,19,17))+
      theme_cowplot()+
      theme(legend.position = "bottom", legend.box = "horizontal")
    
    
    MAPheno245.pca1=ggplot(MA2.245.pca.df, aes(x=PC1, y=PC2, shape = Time, color = Population2, group=Population)) +
      stat_ellipse()+
      geom_point(size=3)+
      #scale_colour_manual(name = "Treatment", 
                          #breaks = c("a", "b"),
                          #labels = c("Low Quality","High Quality"),
                          #values = c("grey","black"))+
      #scale_shape_manual(breaks=c("0", "2", "4"), 
                        # label=c("Founder","Timepoint 2 (Summer)", "Timepoint 4 (Fall)"), 
                         #values =c(21,19,17))+
      theme_cowplot()+
      theme(legend.position = "bottom", legend.box = "horizontal")
    
    MAPheno.PCA <- ggarrange(MAPheno24.pca1, MAPheno245.pca1, ncol=2, labels = c("A","B"), common.legend = FALSE, font.label = list(size = 20))
    
    MAPheno.PCA
    
  
    
    
    #### GETTING VECTORS TO A DATA FRAME####
    ##All Vars 
  
    MA2_C24.pca.vectors<-as.data.frame(MA2_C24.pca$rotation)
    MA2_C24.pca.vectors$ABBRV<-rownames(MA2_C24.pca.vectors)

    MA2_C245.pca.vectors<-as.data.frame(MA2_C245.pca$rotation)
    MA2_C245.pca.vectors$ABBRV<-rownames(MA2_C245.pca.vectors)
    
    MA2_AT24.pca.vectors<-as.data.frame(MA2_AT24.pca$rotation)
    MA2_AT24.pca.vectors$ABBRV<-rownames(MA2_AT24.pca.vectors)
    
    MA2_AT245.pca.vectors<-as.data.frame(MA2_AT245.pca$rotation)
    MA2_AT245.pca.vectors$ABBRV<-rownames(MA2_AT245.pca.vectors)
    
    MA2_LB24.pca.vectors<-as.data.frame(MA2_LB24.pca$rotation)
    MA2_LB24.pca.vectors$ABBRV<-rownames(MA2_LB24.pca.vectors)
    
    MA2_LB245.pca.vectors<-as.data.frame(MA2_LB245.pca$rotation)
    MA2_LB245.pca.vectors$ABBRV<-rownames(MA2_LB245.pca.vectors)
    
    MA2_MA24.pca.vectors<-as.data.frame(MA2_MA24.pca$rotation)
    MA2_MA24.pca.vectors$ABBRV<-rownames(MA2_MA24.pca.vectors)
    
    MA2_MA245.pca.vectors<-as.data.frame(MA2_MA245.pca$rotation)
    MA2_MA245.pca.vectors$ABBRV<-rownames(MA2_MA245.pca.vectors)
    
    
    # CALCULATING VAR EXPLAINED
    ##All Vars 
    MA2_C24.pca.var.explained<-MA2_C24.pca$sdev^2/sum(MA2_C24.pca$sdev^2)
    MA2_AT24.pcaa.var.explained<-MA2_AT24.pca$sdev^2/sum(MA2_AT24.pca$sdev^2)
    MA2_LB24.pcaa.var.explained<-MA2_LB24.pca$sdev^2/sum(MA2_LB24.pca$sdev^2)
    MA2_MA24.pcaa.var.explained<-MA2_MA24.pca$sdev^2/sum(MA2_MA24.pca$sdev^2)
    
    MA24.pca.var.explained<-MA24.pca$sdev^2/sum(MA24.pca$sdev^2)
    
    MA2_C245.pca.var.explained<-MA2_C245.pca$sdev^2/sum(MA2_C245.pca$sdev^2)
    MA2_AT245.pcaa.var.explained<-MA2_AT245.pca$sdev^2/sum(MA2_AT245.pca$sdev^2)
    MA2_LB245.pcaa.var.explained<-MA2_LB245.pca$sdev^2/sum(MA2_LB245.pca$sdev^2)
    MA2_MA245.pcaa.var.explained<-MA2_MA245.pca$sdev^2/sum(MA2_MA245.pca$sdev^2)
    
    MA245.pca.var.explained<-MA245.pca$sdev^2/sum(MA245.pca$sdev^2)
    
    ####calculating hulls####
    # SUBSETTING AND CALCULATING HULLS FOR EACH Population FOR PLOTTING
    MA2C.24.2.pca.df<-subset(MA2C.24.pca.df,Time=="2")
    MA2C.24.2.pca.hull<- MA2C.24.2.pca.df %>% ddply("Time",find_hull)
  
    MA2C.24.4.pca.df<-subset(MA2C.24.pca.df,Time=="4")
    MA2C.24.4.pca.hull<- MA2C.24.4.pca.df %>% ddply("Time",find_hull)
    
    MA2AT.24.2.pca.df<-subset(MA2AT.24.pca.df,Time=="2")
    MA2AT.24.2.pca.hull<- MA2AT.24.2.pca.df %>% ddply("Time",find_hull)
    
    MA2AT.24.4.pca.df<-subset(MA2AT.24.pca.df,Time=="4")
    MA2AT.24.4.pca.hull<- MA2AT.24.4.pca.df %>% ddply("Time",find_hull)
    
    MA2LB.24.2.pca.df<-subset(MA2LB.24.pca.df,Time=="2")
    MA2LB.24.2.pca.hull<- MA2LB.24.2.pca.df %>% ddply("Time",find_hull)
    
    MA2LB.24.4.pca.df<-subset(MA2LB.24.pca.df,Time=="4")
    MA2LB.24.4.pca.hull<- MA2LB.24.4.pca.df %>% ddply("Time",find_hull)
    
    MA2MA.24.2.pca.df<-subset(MA2MA.24.pca.df,Time=="2")
    MA2MA.24.2.pca.hull<- MA2MA.24.2.pca.df %>% ddply("Time",find_hull)
    
    MA2MA.24.4.pca.df<-subset(MA2MA.24.pca.df,Time=="4")
    MA2MA.24.4.pca.hull<- MA2MA.24.4.pca.df %>% ddply("Time",find_hull)
    
  
    
    #### HUlls using TIMWPOINT 2,4,5
    
    MA2C.245.2.pca.df<-subset(MA2C.245.pca.df,Time=="2")
    MA2C.245.2.pca.hull<- MA2C.245.2.pca.df %>% ddply("Time",find_hull)
    
    MA2C.245.4.pca.df<-subset(MA2C.245.pca.df,Time=="4")
    MA2C.245.4.pca.hull<- MA2C.245.4.pca.df %>% ddply("Time",find_hull)
    
    MA2C.245.5.pca.df<-subset(MA2C.245.pca.df,Time=="5")
    MA2C.245.5.pca.hull<- MA2C.245.5.pca.df %>% ddply("Time",find_hull)
    
    MA2AT.245.2.pca.df<-subset(MA2AT.245.pca.df,Time=="2")
    MA2AT.245.2.pca.hull<- MA2AT.245.2.pca.df %>% ddply("Time",find_hull)
    
    MA2AT.245.4.pca.df<-subset(MA2AT.245.pca.df,Time=="4")
    MA2AT.245.4.pca.hull<- MA2AT.245.4.pca.df %>% ddply("Time",find_hull)
    
    MA2AT.245.5.pca.df<-subset(MA2AT.245.pca.df,Time=="5")
    MA2AT.245.5.pca.hull<- MA2AT.245.5.pca.df %>% ddply("Time",find_hull)
    
    MA2LB.245.2.pca.df<-subset(MA2LB.245.pca.df,Time=="2")
    MA2LB.245.2.pca.hull<- MA2LB.245.2.pca.df %>% ddply("Time",find_hull)
    
    MA2LB.245.4.pca.df<-subset(MA2LB.245.pca.df,Time=="4")
    MA2LB.245.4.pca.hull<- MA2LB.245.4.pca.df %>% ddply("Time",find_hull)
    
    MA2LB.245.5.pca.df<-subset(MA2LB.245.pca.df,Time=="5")
    MA2LB.245.5.pca.hull<- MA2LB.245.5.pca.df %>% ddply("Time",find_hull)
   
    MA2MA.245.2.pca.df<-subset(MA2MA.245.pca.df,Time=="2")
    MA2MA.245.2.pca.hull<- MA2MA.245.2.pca.df %>% ddply("Time",find_hull)
    
    MA2MA.245.4.pca.df<-subset(MA2MA.245.pca.df,Time=="4")
    MA2MA.245.4.pca.hull<- MA2MA.245.4.pca.df %>% ddply("Time",find_hull)
    
    MA2MA.245.5.pca.df<-subset(MA2MA.245.pca.df,Time=="5")
    MA2MA.245.5.pca.hull<- MA2MA.245.5.pca.df %>% ddply("Time",find_hull)
    
    
  }
  
  #### PLOT Config####
  {
    # SIZES
    {
      F3_POINTSIZE=2
      F3_LINESIZE=1
      F3_FILLALPHA=0.3
      
      F3_XLLIM=-4.3
      F3_XULIM=3.15
      F3A_YLLIM=-3.0
      F3A_YULIM=4.45
      F3B_YLLIM=-3.45
      F3B_YULIM=4.0
      
      F3_LGDJUST=c(0,1)
      F3_LGDPOS=c(0,0) ###chnage if want legends on indivisual plots to c(0,0)
      F3_LGDSPCX=0.02
      F3_LGDSPCY=0.1
      F3_LGDKEYSIZE=0.6
      
      F3B_ARROWSIZE=.5
      F3_ARROWHEAD=0.2
      F3B_TEXTSIZE=2.5
      F3_VECTORSCALEX=1.2
      F3_VECTORSCALEY=1.2
      F3_VECTOR_ALPHA=0.4
      F3_CIRCLE_COLOR="grey90"
      
      F3_VECTORS_POS_L=0.0
      F3_VECTORS_POS_B=0.62
      F3_VECTORS_POS_T=1.02
      F3_VECTORS_POS_R=0.38
      
      F4_VECTORS_POS_L=1.0
      F4_VECTORS_POS_B=0.62
      F4_VECTORS_POS_T=1.02
      F4_VECTORS_POS_R=1.38
      
      
      F3_PDFHW=7.1
      F3_PDFHH=3.5
      F3_PDFVW=3.5
      F3_PDFVH=7.1
      
      # NEEDED FOR THE PCA VECTORS PLOT
      angle <- seq(-pi, pi, length = 50) 
      circle <- data.frame(x = sin(angle), y = cos(angle))
    }
    
    
    #### TIMEPOINT 2 & 4 ####

    
    MA2_24_MAC<-ggplot()+
      scale_color_manual(values = MA2ColorListPopulation) +
      scale_shape_manual(values = MA2TimeShape) +
      scale_fill_manual(values=MA2ColorListTime)+
      
      geom_polygon(data=MA2C.24.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.24.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2C.24.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.24.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2MA.24.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2MA.24.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2MA.24.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2MA.24.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_point(MA2C.24.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(MA2C.24.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(MA2MA.24.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      geom_point(MA2MA.24.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+

      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3A_YLLIM, F3A_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(MA24.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("PC2 -",percent(MA24.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("C VS. MA, TP 2 AND 4, All Traits") +
      theme(
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(.5, .5, .5,.5), 
        legend.text=element_text(size=10),
      )
    
    MA2_24_MAC
    
    MA2_24_ATLBC<-ggplot()+
      scale_color_manual(values = MA2ColorListPopulation) +
      scale_shape_manual(values = MA2TimeShape) +
      scale_fill_manual(values=MA2ColorListTime)+
      
      geom_polygon(data=MA2C.24.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.24.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2C.24.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.24.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2AT.24.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2AT.24.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2AT.24.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2AT.24.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2LB.24.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2LB.24.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2LB.24.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2LB.24.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_point(MA2C.24.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      geom_point(MA2C.24.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      geom_point(MA2AT.24.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+ 
      geom_point(MA2AT.24.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      geom_point(MA2LB.24.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+ 
      geom_point(MA2LB.24.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3A_YLLIM, F3A_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(MA24.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("PC2 -",percent(MA24.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("C VS. AT VS. LB, TP 2 AND 4, All Traits") +
      theme(
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(.5, .5, .5,.5), 
        legend.text=element_text(size=10),
      )
    
    MA2_24_ATLBC
    
  ####TIMEPOIINTS 2,4,5####
    MA2_245_MAC<-ggplot()+
      scale_color_manual(values = MA2ColorListPopulation) +
      scale_shape_manual(values = MA2TimeShape) +
      scale_fill_manual(values=MA2ColorListTime)+
      
      geom_polygon(data=MA2C.245.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.245.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2C.245.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.245.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2C.245.5.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.245.5.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2MA.245.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2MA.245.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2MA.245.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2MA.245.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2MA.245.5.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=MA2MA.245.5.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_point(MA2C.245.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(MA2C.245.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(MA2C.245.5.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(MA2MA.245.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      geom_point(MA2MA.245.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(MA2MA.245.5.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3A_YLLIM, F3A_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(MA245.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("PC2 -",percent(MA245.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("C VS. MA, TP 2, 4 AND 5 , All Traits") +
      theme(
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(.5, .5, .5,.5), 
        legend.text=element_text(size=10),
      )
    
    MA2_245_MAC
    
    MA2_245_ATLBC<-ggplot()+
      scale_color_manual(values = MA2ColorListPopulation) +
      scale_shape_manual(values = MA2TimeShape) +
      scale_fill_manual(values=MA2ColorListTime)+
      
      geom_polygon(data=MA2C.245.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.245.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2C.245.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.245.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2C.245.5.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2C.245.5.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2AT.245.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2AT.245.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2AT.245.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2AT.245.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2AT.245.5.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2AT.245.5.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2LB.245.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2LB.245.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2LB.245.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2LB.245.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=MA2LB.245.5.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population2)),size=F3_LINESIZE) +
      geom_polygon(data=MA2LB.245.5.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_point(MA2C.245.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      geom_point(MA2C.245.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      geom_point(MA2C.245.5.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      
      geom_point(MA2AT.245.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+ 
      geom_point(MA2AT.245.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      geom_point(MA2AT.245.5.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      
      
      geom_point(MA2LB.245.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+ 
      geom_point(MA2LB.245.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      geom_point(MA2LB.245.5.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population2),size=F3_POINTSIZE)+
      
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3A_YLLIM, F3A_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(MA245.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("PC2 -",percent(MA245.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("C VS. AT VS. LB, TP 2,4, AND 5 , All Traits") +
      theme(
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(.5, .5, .5,.5), 
        legend.text=element_text(size=10),
      )
    
    MA2_245_ATLBC
    
####COMBINE PCA PHENO PLOT 
    
    library(ggpubr)
    
    
    
    MAFig2PCA <- ggarrange(MA2_24_MAC,MA2_245_MAC,MA2_24_ATLBC,MA2_245_ATLBC,
                        labels = c("A", "B", "C", "D"),
                        ncol = 2, nrow = 2, common.legend = FALSE)
    MAFig2PCA
    
    save(Fig2PCA, file = "fig2.rdata")

