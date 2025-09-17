
pckg <- c('tidyverse','lme4','car','piecewiseSEM','lubridate','readxl','openxlsx','smatr','AICcmodavg','MuMIn','ggpubr','corrplot','arm','visreg','lmerTest','lmerTest','viridis','stringr','ggiraphExtra','ggsci','ggpubr','grid')

for(i in 1:length(pckg)) {
  if (!requireNamespace(pckg[i]))
    install.packages(pckg[i])
}


#package----
library(tidyverse)
library(lme4)
library(car)
library(piecewiseSEM)
library(lubridate)
library(readxl)
library(openxlsx)
library(smatr)
library(AICcmodavg)
library(MuMIn)#model selection
library(ggpubr)
library(corrplot)
library(arm)#standardize models
library(visreg)
library(lmerTest)
library(viridis)#included in tidyverse
library(stringr)
library(ggiraphExtra)
library(ggsci); library(ggpubr)
library(grid)# add panel spacing
library(metafor)
library(car)
library(ARTofR)
library(RColorBrewer)
library(scales)
source("stepba function.R")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                                  READ DATA                               ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#final plot data

DataFinal <- read.xlsx('Data/DataSecond.xlsx')%>%
  rename(Diversity = Richness)%>%
  filter(Diversity!=0)%>%
  mutate(Composition =paste(Diversity,Composition,sep='_'))%>%
  mutate(InflorescenceTime  = as.Date(InflorescenceTime, origin = "1899-12-30"))%>%
  mutate(FirstInflorescenceDay= yday(InflorescenceTime)-yday("2022-3-22"))%>%#yday：Julian Day of Year
  mutate(LogInflorescenceBiomass=log(InflorescenceBiomass),LogSeedBiomass=log(SeedBiomass),LogAboveBiomass=log(AboveBiomass),LogTotalAboveBiomass=log(TotalAboveBiomass),LogBelowBiomass=log(BelowBiomass),LogTotalBelowBiomass=log(TotalBelowBiomass),LogTotalBiomass=log(TotalBiomass),LogCormBiomass=log(CormBiomass+0.001),LogRhizomeBiomass=log(RhizomeBiomass),LogFine_rootBiomass=log(Fine_rootBiomass),LogPer_ramet_above=log(Per_ramet_above),LogSeedNumber=log(SeedNumber),LogInflorescenceNumber=log(InflorescenceNumber),LogCormNumber=log(CormNumber+0.1),AsexualBiomass=CormBiomass+RhizomeBiomass,LogAsexualBiomass=log(AsexualBiomass),LogRametNumber=log(RametNumber),InflorescenceSize=InflorescenceBiomass/InflorescenceNumber,LogInflorescenceSize=log(InflorescenceSize),SeedSize=SeedBiomass/SeedNumber,LogSeedSize=log(SeedSize),CormSize=CormBiomass/CormNumber)%>%
  mutate(CormSize = coalesce(CormSize, 0))%>%
  mutate(LogCormSize=log(CormSize+0.001))%>%
  mutate(LogCormNumber=log(CormNumber+0.1))%>%
  rename(AbovegroundVegetativeBiomass=AboveBiomass,LogAbovegroundVegetativeBiomass=LogAboveBiomass)%>%
  mutate(ReproductiveBiomass=InflorescenceBiomass+AsexualBiomass,LogReproductiveBiomass=log(ReproductiveBiomass))%>%
  mutate(VegetativeBiomass=AbovegroundVegetativeBiomass+Fine_rootBiomass,LogVegetativeBiomass=log(VegetativeBiomass))

summary(DataFinal)


#Phenological data
MonthlyData <- read.xlsx('Data/MonthlyData.xlsx')%>%
  rename(Diversity = Richness)%>%
  filter(Diversity!=0)%>%
  arrange(Plot,Composition, Repetition,Month)

MonthlyData<-MonthlyData%>%
    group_by(Plot) %>%
   mutate(Max_Inflorescence = cummax(Inflorescence),Max_Ramet = cummax(Ramet)) %>%
  ungroup() %>%
  dplyr::select(-Inflorescence,-Ramet) %>%
  rename(Inflorescence = Max_Inflorescence,Ramet=Max_Ramet)%>%
  mutate(Composition =paste(Diversity,Composition,sep='_'))%>%
   mutate(Repetition=as.factor(Repetition),Composition=as.factor(Composition),Ramet=as.numeric(Ramet),Inflorescence=as.numeric(Inflorescence))%>%#Diversity=as.factor(Diversity),
  mutate(LogRamet=log(Ramet),LogInflorescence=log(Inflorescence+1),LogHeight=log(MeanHeight),LogDiversity=log(Diversity),MeanHeight.scale=scale(MeanHeight),Diversity.scale=scale(Diversity),Month.scale=scale(Month))%>%
  mutate(Day = plyr::mapvalues(Month, c(4,5,6,7,8,9,10), c("2022-4-22","2022-5-22","2022-6-22","2022-7-22","2022-8-22","2022-9-22","2022-10-22")))%>%
  mutate(Day  = as.Date(Day, origin = "1899-12-30"))%>%
  mutate(Day= yday(Day)-yday("2022-3-22"))%>%#yday：Julian Day of Year
  mutate(Day.scale=scale(Day))

names(MonthlyData)
  
 unique(MonthlyData$Day)
summary(MonthlyData)
unique(MonthlyData$Month)

#soil data
Soil <- read.csv('Data/Soil.csv')%>%
  rename(Diversity = Richness)%>%
    mutate(Composition =paste(Diversity,Composition,sep='_'))
summary(Soil)


PlotTreatment <- read.xlsx('Data/PlotTreatment.xlsx')%>%
  rename(Diversity = Richness)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                   FIT ALLOMETRY OF Ramet/Inflorescence~HEIGHT                  ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#y=ax^b
#log(y)=log(a)+b⋅log(x)

str(MonthlyData)
##~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Ramet~Height  ----
##~~~~~~~~~~~~~~~~~~~~~~

RametLmer <- lmer(LogRamet~LogHeight*Diversity+(1|Composition),MonthlyData)
summary(RametLmer)
anova(RametLmer)


RametLmer.coefs <- fixef(RametLmer)
RametLmer.Slope <- RametLmer.coefs["LogHeight"]
RametLmer.Intercept <- exp(RametLmer.coefs["(Intercept)"])
cat("RametLmer.Slope =", RametLmer.Slope, "\n")
cat("RametLmer.Intercept =", RametLmer.Intercept, "\n")

#get coefficient
RametLmer.Summary <- as.data.frame(coef(summary(RametLmer)))
RametLmer.CI <- confint(RametLmer)[-c(1,2),]
RametLmer.Result <- as.data.frame(cbind(RametLmer.Summary,RametLmer.CI))
RametLmer.Result$Variable <- row.names(RametLmer.Result)
row.names(RametLmer.Result)
RametLmer.Result
RametLmer.Result <- RametLmer.Result %>%
  dplyr::select(Variable,Estimate,lower='2.5 %',upper='97.5 %',P='Pr(>|t|)')%>%
  filter(Variable!='(Intercept)')%>%
  mutate(Y='Ramet~Height')%>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0)))%>%
  mutate(Variable = plyr::mapvalues(Variable, c("LogHeight","Diversity","LogHeight:Diversity"), c("Ln(Height)","Diversity","Ln(Height) × Diversity")))%>%
  mutate(Colour=as.factor(Colour))
RametLmer.Result$Variable <- factor(RametLmer.Result$Variable,levels=c('Ln(Height) × Diversity','Diversity','Ln(Height)'))


#plot effect size
Plot.RametLmer.Result <- ggplot(data=RametLmer.Result,aes(Variable, Estimate, col=Colour)) +
  geom_errorbar(data=RametLmer.Result, mapping=aes(ymin=lower, ymax=upper), width=0, size=1) +
  geom_point(size=2, shape=19)+##0433ff
  geom_hline(aes(yintercept=0), colour="gray", linetype="dashed")+
  theme_bw()+
  scale_colour_manual(name = NULL, values = c('#199b26','#0433ff'))  + 
  scale_fill_manual(name = NULL, values = c('#199b26','#0433ff'))

Plot.RametLmer.Result<-Plot.RametLmer.Result+labs(x ="")+labs(y ="Effect size")+#
  theme(axis.title.x =element_text( size=9, colour="black"),
        axis.title.y=element_text(size=9, colour="black"))+
  theme(axis.text.x =element_text(size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"))+
  theme(panel.grid = element_blank())+
  theme(plot.margin=unit(x=c(0.2,0.4,0.2,0),units="cm"))+
  theme(legend.position="none")+
  theme(strip.text = element_text(size = 7))+
  labs(title = "Ln(Ramet)")+ 
  theme(plot.title = element_text(size = 9,hjust = 0.5, face = "bold"))

Plot.RametLmer.Result<-Plot.RametLmer.Result+coord_flip()#+
Plot.RametLmer.Result

common_legend_theme <- theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(2, 0, 2, 0, unit = "pt"), legend.key.width = unit(10, "pt"), legend.key.height = unit(10, "pt"), legend.text = element_text(size = 8), legend.title = element_text(size = 8))


#new data
new_LogHeight <- expand.grid(
  LogHeight=seq(min(MonthlyData$LogHeight),max(MonthlyData$LogHeight),length=1000),
  Diversity=c(1,2,4,8))%>%
  as.data.frame()

#predict Ramet
pred_RametLmer<-predict(RametLmer,newdata=new_LogHeight,re.form=~0)

ci_line.RametLmer<-bootMer(RametLmer,FUN=function(.) predict(.,newdata=new_LogHeight,
                                                               re.form=~0),nsim=1000)
ci_line.RametLmer<-apply(ci_line.RametLmer$t,2,function(x) x[order(x)][c(25,975)])
lb_RametLmer<-ci_line.RametLmer[1,]
ub_RametLmer<-ci_line.RametLmer[2,]

predict.RametLmer <-cbind(new_LogHeight,pred_RametLmer,lb_RametLmer,ub_RametLmer)%>%
  mutate(MeanHeight=exp(LogHeight),Ramet=exp(pred_RametLmer),lb_Ramet=exp(lb_RametLmer), ub_Ramet=exp(ub_RametLmer))




SpeciesCol%>%show_col

SpeciesCol <- brewer.pal(7, "BuPu")[c(3,4,5,6)]

SpeciesCol%>%show_col

fig.log.RametLmer <- ggplot(predict.RametLmer)+
  geom_point(
    data = MonthlyData,
    aes(x = LogHeight, y = LogRamet, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.1),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=LogHeight,y=pred_RametLmer,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=LogHeight,ymin=lb_RametLmer,ymax=ub_RametLmer,fill=as.factor(Diversity)),alpha=0.3)+
  xlab('Ln(Height (cm)) ')+
  ylab('Ln(Ramet)')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  common_legend_theme
fig.log.RametLmer


fig.RametLmer <- ggplot(predict.RametLmer)+
  geom_point(
    data = MonthlyData,
    aes(x = MeanHeight, y = Ramet, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.5),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=MeanHeight,y=Ramet,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=MeanHeight,ymin=lb_Ramet,ymax=ub_Ramet,fill=as.factor(Diversity)),alpha=0.2)+
  xlab('Height (cm)')+
  ylab('Ramet')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  common_legend_theme
fig.RametLmer




##~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Inflorescence~Height  ----
##~~~~~~~~~~~~~~~~~~~~~~

InflorescenceLmer<- lmer(LogInflorescence~LogHeight*Diversity+(1|Composition),MonthlyData)

summary(InflorescenceLmer)


InflorescenceLmer.coefs <- fixef(InflorescenceLmer)
InflorescenceLmer.Slope <- InflorescenceLmer.coefs["LogHeight"]
InflorescenceLmer.Intercept <- exp(InflorescenceLmer.coefs["(Intercept)"]) - 1
cat("InflorescenceLmer.Slope =", InflorescenceLmer.Slope, "\n")
cat("InflorescenceLmer.Intercept =", InflorescenceLmer.Intercept, "\n")

#get coefficient
InflorescenceLmer.Summary <- as.data.frame(coef(summary(InflorescenceLmer)))
InflorescenceLmer.CI <- confint(InflorescenceLmer)[-c(1,2),]
InflorescenceLmer.Result <- as.data.frame(cbind(InflorescenceLmer.Summary,InflorescenceLmer.CI))
InflorescenceLmer.Result$Variable <- row.names(InflorescenceLmer.Result)
row.names(InflorescenceLmer.Result)
InflorescenceLmer.Result
InflorescenceLmer.Result <- InflorescenceLmer.Result %>%
  dplyr::select(Variable,Estimate,lower='2.5 %',upper='97.5 %',P='Pr(>|t|)')%>%
  filter(Variable!='(Intercept)')%>%
  mutate(Y='Inflorescence~Height')%>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0)))%>%
  mutate(Variable = plyr::mapvalues(Variable, c("LogHeight","Diversity","LogHeight:Diversity"), c("Ln(Height)","Diversity","Ln(Height) × Diversity")))%>%
  mutate(Colour=as.factor(Colour))
InflorescenceLmer.Result$Variable <- factor(InflorescenceLmer.Result$Variable,levels=c('Ln(Height) × Diversity','Diversity','Ln(Height)'))


#plot effect size
Plot.InflorescenceLmer.Result <- ggplot(data=InflorescenceLmer.Result,aes(Variable, Estimate, col=Colour)) +
  geom_errorbar(data=InflorescenceLmer.Result, mapping=aes(ymin=lower, ymax=upper), width=0, size=1) +
  geom_point(size=2, shape=19)+
  geom_hline(aes(yintercept=0), colour="gray", linetype="dashed")+
  theme_bw()+
  scale_colour_manual(name = NULL, values = c('gray','#0433ff'))  + 
  scale_fill_manual(name = NULL, values = c('gray','#0433ff'))

Plot.InflorescenceLmer.Result<-Plot.InflorescenceLmer.Result+labs(x ="")+labs(y ="Effect size")+#
  theme(axis.title.x =element_text( size=9, colour="black"),
        axis.title.y=element_text(size=9, colour="black"))+
  theme(axis.text.x =element_text(size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"))+
  theme(panel.grid = element_blank())+
  theme(plot.margin=unit(x=c(0.2,0.4,0.2,0),units="cm"))+
  theme(legend.position="none")+
  theme(strip.text = element_text(size = 7))+
  labs(title = "Ln(Inflorescence)")+ 
  theme(plot.title = element_text(size = 9,hjust = 0.5, face = "bold"))

Plot.InflorescenceLmer.Result<-Plot.InflorescenceLmer.Result+coord_flip()
Plot.InflorescenceLmer.Result




#predict Inflorescence
pred_InflorescenceLmer<-predict(InflorescenceLmer,newdata=new_LogHeight,re.form=~0)

ci_line.InflorescenceLmer<-bootMer(InflorescenceLmer,FUN=function(.) predict(.,newdata=new_LogHeight,
                                                                                re.form=~0),nsim=1000)
ci_line.InflorescenceLmer<-apply(ci_line.InflorescenceLmer$t,2,function(x) x[order(x)][c(25,975)])
lb_InflorescenceLmer<-ci_line.InflorescenceLmer[1,]
ub_InflorescenceLmer<-ci_line.InflorescenceLmer[2,]

predict.InflorescenceLmer <-cbind(new_LogHeight,pred_InflorescenceLmer,lb_InflorescenceLmer,ub_InflorescenceLmer)%>%
  mutate(MeanHeight=exp(LogHeight),Inflorescence=exp(pred_InflorescenceLmer)-1,lb_Inflorescence=exp(lb_InflorescenceLmer)-1, ub_Inflorescence=exp(ub_InflorescenceLmer)-1)



fig.log.InflorescenceLmer <- ggplot(predict.InflorescenceLmer)+
  geom_point(
    data = MonthlyData,
    aes(x = LogHeight, y = LogInflorescence, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.1),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=LogHeight,y=pred_InflorescenceLmer,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=LogHeight,ymin=lb_InflorescenceLmer,ymax=ub_InflorescenceLmer,fill=as.factor(Diversity)),alpha=0.2)+
  xlab('Ln(Height (cm)) ')+
  ylab('Ln(Inflorescence)')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  common_legend_theme
fig.log.InflorescenceLmer

names(MonthlyData)
fig.InflorescenceLmer <- ggplot(predict.InflorescenceLmer)+
  geom_point(
    data = MonthlyData,
    aes(x = MeanHeight, y = Inflorescence, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.5),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=MeanHeight,y=Inflorescence,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=MeanHeight,ymin=lb_Inflorescence,ymax=ub_Inflorescence,fill=as.factor(Diversity)),alpha=0.2)+
  xlab('Height (cm)')+
  ylab('Inflorescence')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)+ 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  common_legend_theme
fig.InflorescenceLmer



##~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Inflorescence~Ramet  ----
##~~~~~~~~~~~~~~~~~~~~~~
InflorescenceLmer.Ramet<- lmer(LogInflorescence~LogRamet*Diversity+(1|Composition),MonthlyData)
summary(InflorescenceLmer.Ramet)


InflorescenceLmer.Ramet.coefs <- fixef(InflorescenceLmer.Ramet)
InflorescenceLmer.Ramet.Slope <- InflorescenceLmer.Ramet.coefs["LogRamet"]
InflorescenceLmer.Ramet.Intercept <- exp(InflorescenceLmer.Ramet.coefs["(Intercept)"]) - 1
cat("InflorescenceLmer.Ramet.Slope =", InflorescenceLmer.Ramet.Slope, "\n")
cat("InflorescenceLmer.Ramet.Intercept =", InflorescenceLmer.Ramet.Intercept, "\n")

#get coefficient
InflorescenceLmer.Ramet.Summary <- as.data.frame(coef(summary(InflorescenceLmer.Ramet)))

InflorescenceLmer.Ramet.CI <- confint(InflorescenceLmer.Ramet)[-c(1,2),]
InflorescenceLmer.Ramet.Result <- as.data.frame(cbind(InflorescenceLmer.Ramet.Summary,InflorescenceLmer.Ramet.CI))
InflorescenceLmer.Ramet.Result$Variable <- row.names(InflorescenceLmer.Ramet.Result)
row.names(InflorescenceLmer.Ramet.Result)
InflorescenceLmer.Ramet.Result
InflorescenceLmer.Ramet.Result <- InflorescenceLmer.Ramet.Result %>%
  dplyr::select(Variable,Estimate,lower='2.5 %',upper='97.5 %',P='Pr(>|t|)')%>%
  filter(Variable!='(Intercept)')%>%
  mutate(Y='Inflorescence~Ramet')%>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0)))%>%
  mutate(Variable = plyr::mapvalues(Variable, c("LogRamet","Diversity","LogRamet:Diversity"), c("Ln(Ramet)","Diversity","Ln(Ramet) × Diversity")))%>%
  mutate(Colour=as.factor(Colour))
InflorescenceLmer.Ramet.Result$Variable <- factor(InflorescenceLmer.Ramet.Result$Variable,levels=c('Ln(Ramet) × Diversity','Diversity','Ln(Ramet)'))


#plot effect size
Plot.InflorescenceLmer.Ramet.Result <- ggplot(data=InflorescenceLmer.Ramet.Result,aes(Variable, Estimate, col=Colour)) +
  geom_errorbar(data=InflorescenceLmer.Ramet.Result, mapping=aes(ymin=lower, ymax=upper), width=0, size=1) +
  geom_point(size=2, shape=19)+
  geom_hline(aes(yintercept=0), colour="gray", linetype="dashed")+
  theme_bw()+
  scale_colour_manual(name = NULL, values = c('#199b26','gray','#0433ff'))  + 
  scale_fill_manual(name = NULL, values = c('#199b26','gray','#0433ff'))

Plot.InflorescenceLmer.Ramet.Result<-Plot.InflorescenceLmer.Ramet.Result+labs(x ="")+labs(y ="Effect size")+
  theme(axis.title.x =element_text( size=9, colour="black"),
        axis.title.y=element_text(size=9, colour="black"))+
  theme(axis.text.x =element_text(size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"))+
  theme(panel.grid = element_blank())+
  theme(plot.margin=unit(x=c(0.2,0.4,0.2,0),units="cm"))+
  theme(legend.position="none")+
  theme(strip.text = element_text(size = 7))+
  labs(title = "Ln(Inflorescence)")+ 
  theme(plot.title = element_text(size = 9,hjust = 0.5, face = "bold"))

Plot.InflorescenceLmer.Ramet.Result<-Plot.InflorescenceLmer.Ramet.Result+coord_flip()#+
Plot.InflorescenceLmer.Ramet.Result



#new data
new_LogRamet<- expand.grid(
  LogRamet=seq(min(MonthlyData$LogRamet),max(MonthlyData$LogRamet),length=1000),
  Diversity=c(1,2,4,8))%>%
  as.data.frame()


#predict Inflorescence
pred_InflorescenceLmer.Ramet<-predict(InflorescenceLmer.Ramet,newdata=new_LogRamet,re.form=~0)

ci_line.InflorescenceLmer.Ramet<-bootMer(InflorescenceLmer.Ramet,FUN=function(.) predict(.,newdata=new_LogRamet,
                                                               re.form=~0),nsim=1000)
ci_line.InflorescenceLmer.Ramet<-apply(ci_line.InflorescenceLmer.Ramet$t,2,function(x) x[order(x)][c(25,975)])
lb_InflorescenceLmer.Ramet<-ci_line.InflorescenceLmer.Ramet[1,]
ub_InflorescenceLmer.Ramet<-ci_line.InflorescenceLmer.Ramet[2,]

predict.InflorescenceLmer.Ramet <-cbind(new_LogRamet,pred_InflorescenceLmer.Ramet,lb_InflorescenceLmer.Ramet,ub_InflorescenceLmer.Ramet)%>%
  mutate(Ramet=exp(LogRamet),Inflorescence=exp(pred_InflorescenceLmer.Ramet)-1,lb_Inflorescence=exp(lb_InflorescenceLmer.Ramet)-1, ub_Inflorescence=exp(ub_InflorescenceLmer.Ramet)-1)



fig.log.InflorescenceLmer.Ramet <- ggplot(predict.InflorescenceLmer.Ramet)+
  geom_point(
    data = MonthlyData,
    aes(x = LogRamet, y = LogInflorescence, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.1), 
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=LogRamet,y=pred_InflorescenceLmer.Ramet,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=LogRamet,ymin=lb_InflorescenceLmer.Ramet,ymax=ub_InflorescenceLmer.Ramet,fill=as.factor(Diversity)),alpha=0.2)+
  xlab('Ln(Ramet) ')+
  ylab('Ln(Inflorescence)')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  common_legend_theme

fig.log.InflorescenceLmer.Ramet


fig.InflorescenceLmer.Ramet <- ggplot(predict.InflorescenceLmer.Ramet)+
  geom_point(
    data = MonthlyData,
    aes(x = Ramet, y = Inflorescence, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.5),  
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=Ramet,y=Inflorescence,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=Ramet,ymin=lb_Inflorescence,ymax=ub_Inflorescence,fill=as.factor(Diversity)),alpha=0.2)+
  xlab('Ramet')+
  ylab('Inflorescence')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  common_legend_theme

fig.InflorescenceLmer.Ramet

##~~~~~~~~~~~~~~~~~~~~~~
##  ~ combine plot  ----
##~~~~~~~~~~~~~~~~~~~~~~
fig.allometry <- ggpubr::ggarrange(ggpubr::ggarrange(Plot.RametLmer.Result, fig.log.RametLmer + theme(legend.position="none"), fig.RametLmer + theme(legend.position="none"), Plot.InflorescenceLmer.Result, fig.log.InflorescenceLmer + theme(legend.position="none"), fig.InflorescenceLmer + theme(legend.position="none"), Plot.InflorescenceLmer.Ramet.Result, fig.log.InflorescenceLmer.Ramet + theme(legend.position="none"), fig.InflorescenceLmer.Ramet + theme(legend.position="none"), ncol=3, nrow=3, labels=c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'), font.label=list(size=10, color="black", face="bold"), widths=c(1.1,0.9,0.9), common.legend=FALSE), cowplot::get_legend(fig.log.RametLmer + guides(colour=guide_legend(nrow=1, ncol=4, byrow=TRUE), fill=guide_legend(nrow=1, ncol=4, byrow=TRUE)) + theme(legend.position="bottom", legend.box="horizontal", legend.margin=margin(2,0,2,0,"pt"), legend.key.width=unit(10,"pt"), legend.key.height=unit(10,"pt"), legend.text=element_text(size=8), legend.title=element_text(size=8))), ncol=1, heights=c(1,0.08))



fig.allometry <- ggpubr::ggarrange(Plot.RametLmer.Result,fig.log.RametLmer,fig.RametLmer,Plot.InflorescenceLmer.Result,fig.log.InflorescenceLmer,fig.InflorescenceLmer,Plot.InflorescenceLmer.Ramet.Result,fig.log.InflorescenceLmer.Ramet,fig.InflorescenceLmer.Ramet,ncol = 3,nrow=3,labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'),font.label = list(size = 10, color = "black", face = "bold"), widths = c(1.1, 0.9,0.9),common.legend = TRUE,legend = "bottom")
fig.allometry#Figure 4
#ggsave('Result/fig.allometry20250916.pdf',width = 18,height =17,units = c('cm'))




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##              THE INFLUENCE OF Diversity ON THE PEAK             ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Peak number of Inflorescences  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PeakInflorescence <- MonthlyData %>%
  group_by(Plot) %>%
  filter(Inflorescence == max(Inflorescence)) %>%# Peak number of Inflorescence
  slice_min(order_by = Day, n = 1)

summary_PeakInflorescence <- PeakInflorescence %>%
  mutate(Diversity<-as.character(Diversity))%>%
  group_by(Diversity) %>%
  summarise(
    meanDay = mean(Day),
    lowerDay = mean(Day) - qt(0.975, df = n() - 1) * sd(Day) / sqrt(n()),
    upperDay = mean(Day) + qt(0.975, df = n() - 1) * sd(Day) / sqrt(n()),
    meanInflorescence = mean(Inflorescence),
    lowerInflorescence = mean(Inflorescence) - qt(0.975, df = n() - 1) * sd(Inflorescence) / sqrt(n()),
    upperInflorescence= mean(Inflorescence) + qt(0.975, df = n() - 1) * sd(Inflorescence) / sqrt(n()),
    
  )%>%
  mutate(Diversity=as.numeric(Diversity))



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Peak Inflorescence Day  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PeakInflorescenceDayLmer <- lmer(Day~Diversity+(1|Composition),PeakInflorescence)
summary(PeakInflorescenceDayLmer)
MuMIn::r.squaredGLMM(PeakInflorescenceDayLmer)


#predict Peak flowering Day
pred_PeakInflorescenceDayLmer<-predict(PeakInflorescenceDayLmer,newdata=new_Diversity.PeakInflorescence,re.form=~0)

ci_line.PeakInflorescenceDayLmer<-bootMer(PeakInflorescenceDayLmer,FUN=function(.) predict(.,newdata=new_Diversity.PeakInflorescence,
                                                                                       re.form=~0),nsim=1000)
ci_line.PeakInflorescenceDayLmer<-apply(ci_line.PeakInflorescenceDayLmer$t,2,function(x) x[order(x)][c(25,975)])
lb_PeakInflorescenceDayLmer<-ci_line.PeakInflorescenceDayLmer[1,]
ub_PeakInflorescenceDayLmer<-ci_line.PeakInflorescenceDayLmer[2,]

predict.PeakInflorescenceDayLmer <-cbind(new_Diversity.PeakInflorescence,pred_PeakInflorescenceDayLmer,lb_PeakInflorescenceDayLmer,ub_PeakInflorescenceDayLmer)

summary(predict.PeakInflorescenceDayLmer)
summary(PeakInflorescence)

fig.PeakInflorescenceDayLmer<- ggplot()+
  #geom_point(aes(PeakInflorescence$Diversity,PeakInflorescence$Day),size=1,colour='#199b26',alpha=0.3)+#alpha = 1 / 4,
  #geom_jitter(aes(PeakInflorescence$Diversity,PeakInflorescence$Day),size=2,colour='#199b26',alpha=0.3,width = 0.2, height = 0)+#shape=1
  geom_point(data = summary_PeakInflorescence, aes(x = Diversity, y = meanDay), color = "#199b26", size = 2) +
  geom_errorbar(data = summary_PeakInflorescence, aes(x = Diversity, ymin = lowerDay, ymax = upperDay), width = 0.2, color = "#199b26") +
  geom_line(aes(x=predict.PeakInflorescenceDayLmer$Diversity,y=predict.PeakInflorescenceDayLmer$pred_PeakInflorescenceDayLmer),size=1,colour='#0433ff',linetype="dashed")+
  geom_ribbon(aes(x=predict.PeakInflorescenceDayLmer$Diversity,ymin=predict.PeakInflorescenceDayLmer$lb_PeakInflorescenceDayLmer,ymax=predict.PeakInflorescenceDayLmer$ub_PeakInflorescenceDayLmer),alpha=0.3,fill='#0433ff')+
  ylab('Peak inflorescence day')+
  xlab('Genotypic diversity')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  theme(legend.position="none")+
  scale_x_continuous(breaks =  c(1, 2, 4, 8))+
  #scale_y_continuous(limits = c(9, 10.3),breaks = c(9, 10))+
  annotate("text", x = 6.3, y = 199,label = "paste(italic(P), \" = 0.255\")", parse = TRUE,size = 2.5)+
  annotate("text", x = 6.3, y = 196,label = "paste(italic(R) [C] ^ 2 , \" = 0.171\")", parse = TRUE,size = 2.5)+
  annotate("text", x = 6.3, y = 193,label = "paste(italic(R) [M] ^ 2 , \" = 0.044\")", parse = TRUE,size = 2.5)
fig.PeakInflorescenceDayLmer


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ First Inflorescence day ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Date of first inflorescence appearanc

summary_FirstInflorescenceDay <- DataFinal %>%
  mutate(Diversity<-as.character(Diversity))%>%
  group_by(Diversity) %>%
  summarise(
    meanFirstInflorescenceDay = mean(FirstInflorescenceDay),
    lowerFirstInflorescenceDay = mean(FirstInflorescenceDay) - qt(0.975, df = n() - 1) * sd(FirstInflorescenceDay) / sqrt(n()),
    upperFirstInflorescenceDay = mean(FirstInflorescenceDay) + qt(0.975, df = n() - 1) * sd(FirstInflorescenceDay) / sqrt(n()),
  )%>%
  mutate(Diversity=as.numeric(Diversity))

DataFinal

FirstInflorescenceDayLmer <- lmer(FirstInflorescenceDay~Diversity+(1|Composition),DataFinal)
summary(FirstInflorescenceDayLmer)
anova(FirstInflorescenceDayLmer)
MuMIn::r.squaredGLMM(FirstInflorescenceDayLmer)

#new data
new_Diversity.first <- expand.grid(
  Diversity=seq(min(DataFinal$Diversity),max(DataFinal$Diversity),length=1000))%>%
  as.data.frame()
#predict Peak flowering Day
pred_FirstInflorescenceDayLmer<-predict(FirstInflorescenceDayLmer,newdata=new_Diversity.first,re.form=~0)

ci_line.FirstInflorescenceDayLmer<-bootMer(FirstInflorescenceDayLmer,FUN=function(.) predict(.,newdata=new_Diversity.first,
                                                                                     re.form=~0),nsim=1000)
ci_line.FirstInflorescenceDayLmer<-apply(ci_line.FirstInflorescenceDayLmer$t,2,function(x) x[order(x)][c(25,975)])
lb_FirstInflorescenceDayLmer<-ci_line.FirstInflorescenceDayLmer[1,]
ub_FirstInflorescenceDayLmer<-ci_line.FirstInflorescenceDayLmer[2,]

predict.FirstInflorescenceDayLmer <-cbind(new_Diversity.first,pred_FirstInflorescenceDayLmer,lb_FirstInflorescenceDayLmer,ub_FirstInflorescenceDayLmer)

summary(predict.FirstInflorescenceDayLmer)

fig.FirstInflorescenceDayLmer<- ggplot()+
  #geom_point(aes(PeakInflorescence$Diversity,PeakInflorescence$Day),size=1,colour='#199b26',alpha=0.3)+#alpha = 1 / 4,
  #geom_jitter(aes(DataFinal$Diversity,DataFinal$FirstInflorescenceDay),size=2,colour='#199b26',alpha=0.3,width = 0.1, height = 0)+#shape=1,
  geom_point(data = summary_FirstInflorescenceDay, aes(x = Diversity, y = meanFirstInflorescenceDay), color = "#199b26", size = 2) +
  geom_errorbar(data = summary_FirstInflorescenceDay, aes(x = Diversity, ymin = lowerFirstInflorescenceDay, ymax = upperFirstInflorescenceDay), width = 0.2, color = "#199b26") +
  geom_line(aes(x=predict.FirstInflorescenceDayLmer$Diversity,y=predict.FirstInflorescenceDayLmer$pred_FirstInflorescenceDayLmer),size=1,colour='#0433ff')+
  geom_ribbon(aes(x=predict.FirstInflorescenceDayLmer$Diversity,ymin=predict.FirstInflorescenceDayLmer$lb_FirstInflorescenceDayLmer,ymax=predict.FirstInflorescenceDayLmer$ub_FirstInflorescenceDayLmer),alpha=0.3,fill='#0433ff')+
  ylab('The first inflorescence day')+
  xlab('Genotypic diversity')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  theme(legend.position="none")+
  scale_x_continuous(breaks =  c(1, 2, 4, 8))+
  #scale_y_continuous(breaks = c(8,9, 10))+
  annotate("text", x = 6.3, y = 65,label = "paste(italic(P), \" = 0.002\")", parse = TRUE,size = 2.5)+
  annotate("text", x = 6.3, y = 63,label = "paste(italic(R) [C] ^ 2 , \" = 0.118\")", parse = TRUE,size = 2.5)+
  annotate("text", x = 6.3, y = 61,label = "paste(italic(R) [M] ^ 2 , \" = 0.118\")", parse = TRUE,size = 2.5)
fig.FirstInflorescenceDayLmer


##~~~~~~~~~~~~~~~~~~~~~
##  ~ Peak Height  ----
##~~~~~~~~~~~~~~~~~~~~~
PeakMeanHeight<- MonthlyData %>%
  group_by(Plot) %>%
  filter(MeanHeight == max(MeanHeight)) %>%
  slice_min(order_by = Day, n = 1)

summary_PeakMeanHeight <- PeakMeanHeight %>%
  mutate(Diversity<-as.character(Diversity))%>%
  group_by(Diversity) %>%
  summarise(
    meanDay = mean(Day),
    lowerDay = mean(Day) - qt(0.975, df = n() - 1) * sd(Day) / sqrt(n()),
    upperDay = mean(Day) + qt(0.975, df = n() - 1) * sd(Day) / sqrt(n()),
    meanMeanHeight = mean(MeanHeight),
    lowerMeanHeight = mean(MeanHeight) - qt(0.975, df = n() - 1) * sd(MeanHeight) / sqrt(n()),
    upperMeanHeight= mean(MeanHeight) + qt(0.975, df = n() - 1) * sd(MeanHeight) / sqrt(n()),
    
  )%>%
  mutate(Diversity=as.numeric(Diversity))

hist(PeakMeanHeight$MeanHeight)
PeakMeanHeightLmer <- lmer(MeanHeight~Diversity+(1|Composition),PeakMeanHeight)
summary(PeakMeanHeightLmer)
MuMIn::r.squaredGLMM(PeakMeanHeightLmer)

#new data
new_Diversity.PeakMeanHeight <- expand.grid(
  Diversity=seq(min(PeakMeanHeight$Diversity),max(PeakMeanHeight$Diversity),length=1000))%>%
  as.data.frame()

#predict Peak flowering Day
pred_PeakMeanHeightLmer<-predict(PeakMeanHeightLmer,newdata=new_Diversity.PeakMeanHeight,re.form=~0)

ci_line.PeakMeanHeightLmer<-bootMer(PeakMeanHeightLmer,FUN=function(.) predict(.,newdata=new_Diversity.PeakMeanHeight,
                                                                                     re.form=~0),nsim=1000)
ci_line.PeakMeanHeightLmer<-apply(ci_line.PeakMeanHeightLmer$t,2,function(x) x[order(x)][c(25,975)])
lb_PeakMeanHeightLmer<-ci_line.PeakMeanHeightLmer[1,]
ub_PeakMeanHeightLmer<-ci_line.PeakMeanHeightLmer[2,]

predict.PeakMeanHeightLmer <-cbind(new_Diversity.PeakMeanHeight,pred_PeakMeanHeightLmer,lb_PeakMeanHeightLmer,ub_PeakMeanHeightLmer)


summary(predict.PeakMeanHeightLmer)
summary(PeakMeanHeight)

fig.PeakMeanHeightLmer<- ggplot()+
  #geom_point(aes(PeakMeanHeight$Diversity,PeakMeanHeight$Day),size=1,colour='#199b26',alpha=0.3)+#alpha = 1 / 4,
  #geom_jitter(aes(PeakMeanHeight$Diversity,PeakMeanHeight$MeanHeight),size=2,colour='#199b26',alpha=0.3,width = 0.1, height = 0)+#shape=1,
  geom_point(data = summary_PeakMeanHeight, aes(x = Diversity, y = meanMeanHeight), color = "#199b26", size = 2) +
  geom_errorbar(data = summary_PeakMeanHeight, aes(x = Diversity, ymin = lowerMeanHeight, ymax = upperMeanHeight), width = 0.2, color = "#199b26") +
  
  geom_line(aes(x=predict.PeakMeanHeightLmer$Diversity,y=predict.PeakMeanHeightLmer$pred_PeakMeanHeightLmer),size=1,colour='#0433ff',linetype='dashed')+
  geom_ribbon(aes(x=predict.PeakMeanHeightLmer$Diversity,ymin=predict.PeakMeanHeightLmer$lb_PeakMeanHeightLmer,ymax=predict.PeakMeanHeightLmer$ub_PeakMeanHeightLmer),alpha=0.3,fill='#0433ff')+
  ylab('Peak height (cm)')+
  xlab('Genotypic diversity')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  theme(legend.position="none")+
  scale_x_continuous(breaks =  c(1, 2, 4, 8))+
  #scale_y_continuous(breaks = c(8,9, 10))+
  annotate("text", x = 2.2, y = 45.4,label = "paste(italic(P), \" = 0.201\")", parse = TRUE,size = 2.5)+
  annotate("text", x = 2.2, y = 44.5,label = "paste(italic(R) [C] ^ 2 , \" = 0.021\")", parse = TRUE,size = 2.5)+
  annotate("text", x = 2.2, y = 43.6,label = "paste(italic(R) [M] ^ 2 , \" = 0.021\")", parse = TRUE,size = 2.5)
fig.PeakMeanHeightLmer

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Peak Height Day  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PeakMeanHeightDayLmer <- lmer(Day~Diversity+(1|Composition),PeakMeanHeight)
summary(PeakMeanHeightDayLmer)
MuMIn::r.squaredGLMM(PeakMeanHeightDayLmer)


#predict Peak flowering Day
pred_PeakMeanHeightDayLmer<-predict(PeakMeanHeightDayLmer,newdata=new_Diversity.PeakMeanHeight,re.form=~0)

ci_line.PeakMeanHeightDayLmer<-bootMer(PeakMeanHeightDayLmer,FUN=function(.) predict(.,newdata=new_Diversity.PeakMeanHeight,
                                                                                               re.form=~0),nsim=1000)
ci_line.PeakMeanHeightDayLmer<-apply(ci_line.PeakMeanHeightDayLmer$t,2,function(x) x[order(x)][c(25,975)])
lb_PeakMeanHeightDayLmer<-ci_line.PeakMeanHeightDayLmer[1,]
ub_PeakMeanHeightDayLmer<-ci_line.PeakMeanHeightDayLmer[2,]

predict.PeakMeanHeightDayLmer <-cbind(new_Diversity.PeakMeanHeight,pred_PeakMeanHeightDayLmer,lb_PeakMeanHeightDayLmer,ub_PeakMeanHeightDayLmer)

summary(predict.PeakMeanHeightDayLmer)
summary(PeakMeanHeight)


##~~~~~~~~~~~~~~~~~~~~~
##  ~ Peak Ramet  ----
##~~~~~~~~~~~~~~~~~~~~~
PeakRamet <- MonthlyData %>%
  group_by(Plot) %>%
  filter(Ramet == max(Ramet)) %>%
  slice_min(order_by = Day, n = 1)
names(PeakRamet)

summary(PeakRamet)

summary_PeakRamet <- PeakRamet %>%
  mutate(Diversity<-as.character(Diversity))%>%
  group_by(Diversity) %>%
  summarise(
    meanDay = mean(Day),
    lowerDay = mean(Day) - qt(0.975, df = n() - 1) * sd(Day) / sqrt(n()),
    upperDay = mean(Day) + qt(0.975, df = n() - 1) * sd(Day) / sqrt(n()),
    meanRamet = mean(Ramet),
    lowerRamet = mean(Ramet) - qt(0.975, df = n() - 1) * sd(Ramet) / sqrt(n()),
    upperRamet= mean(Ramet) + qt(0.975, df = n() - 1) * sd(Ramet) / sqrt(n()),
    
  )%>%
  mutate(Diversity=as.numeric(Diversity))

hist(PeakRamet$Ramet)
hist(PeakRamet$LogRamet)
PeakRametLmer <- lmer(Ramet~Diversity+(1|Composition),PeakRamet)
summary(PeakRametLmer)
MuMIn::r.squaredGLMM(PeakRametLmer)

#new data
new_Diversity.PeakRamet <- expand.grid(
  Diversity=seq(min(PeakRamet$Diversity),max(PeakRamet$Diversity),length=1000))%>%
  as.data.frame()

#predict Peak flowering Day
pred_PeakRametLmer<-predict(PeakRametLmer,newdata=new_Diversity.PeakRamet,re.form=~0)

ci_line.PeakRametLmer<-bootMer(PeakRametLmer,FUN=function(.) predict(.,newdata=new_Diversity.PeakRamet,
                                                                               re.form=~0),nsim=1000)
ci_line.PeakRametLmer<-apply(ci_line.PeakRametLmer$t,2,function(x) x[order(x)][c(25,975)])
lb_PeakRametLmer<-ci_line.PeakRametLmer[1,]
ub_PeakRametLmer<-ci_line.PeakRametLmer[2,]

predict.PeakRametLmer <-cbind(new_Diversity.PeakRamet,pred_PeakRametLmer,lb_PeakRametLmer,ub_PeakRametLmer)


summary(predict.PeakRametLmer)
summary(PeakRamet)
summary(summary_PeakRamet)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Peak Ramet Day  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PeakRametDayLmer <- lmer(Day~Diversity+(1|Composition),PeakRamet)
summary(PeakRametDayLmer)
MuMIn::r.squaredGLMM(PeakRametDayLmer)


#predict Peak ramet Day
pred_PeakRametDayLmer<-predict(PeakRametDayLmer,newdata=new_Diversity.PeakRamet,re.form=~0)

ci_line.PeakRametDayLmer<-bootMer(PeakRametDayLmer,FUN=function(.) predict(.,newdata=new_Diversity.PeakRamet,
                                                                                         re.form=~0),nsim=1000)
ci_line.PeakRametDayLmer<-apply(ci_line.PeakRametDayLmer$t,2,function(x) x[order(x)][c(25,975)])
lb_PeakRametDayLmer<-ci_line.PeakRametDayLmer[1,]
ub_PeakRametDayLmer<-ci_line.PeakRametDayLmer[2,]

predict.PeakRametDayLmer <-cbind(new_Diversity.PeakRamet,pred_PeakRametDayLmer,lb_PeakRametDayLmer,ub_PeakRametDayLmer)

summary(predict.PeakRametDayLmer)
summary(PeakRamet)



##~~~~~~~~~~~~~~~~~~~~~~
##  ~ combine plot  ----
##~~~~~~~~~~~~~~~~~~~~~~


####plot together####

names(PeakMeanHeight)
fig.PeakDayLmer <- ggplot() +
  # Peak Mean Height
  geom_point(data = summary_PeakMeanHeight, aes(x = Diversity, y = meanDay, color = "Peak height day"), size = 2) +  
  geom_errorbar(data = summary_PeakMeanHeight, aes(x = Diversity, ymin = lowerDay, ymax = upperDay, color = "Peak height day"), width = 0.2) +
  geom_line(aes(x = predict.PeakMeanHeightDayLmer$Diversity, y = predict.PeakMeanHeightDayLmer$pred_PeakMeanHeightDayLmer, color = "Peak height day"), size = 1) +
  geom_ribbon(aes(x = predict.PeakMeanHeightDayLmer$Diversity, ymin = predict.PeakMeanHeightDayLmer$lb_PeakMeanHeightDayLmer, ymax = predict.PeakMeanHeightDayLmer$ub_PeakMeanHeightDayLmer, fill = "Peak height day"), alpha = 0.3) +
  
  # Peak Ramet 
  geom_point(data = summary_PeakRamet, aes(x = Diversity, y = meanDay, color = "Peak ramet day"), size = 2) + 
  geom_errorbar(data = summary_PeakRamet, aes(x = Diversity, ymin = lowerDay, ymax = upperDay, color = "Peak ramet day"), width = 0.2) +
  geom_line(aes(x = predict.PeakRametDayLmer$Diversity, y = predict.PeakRametDayLmer$pred_PeakRametDayLmer, color = "Peak ramet day"), size = 1) +
  geom_ribbon(aes(x = predict.PeakRametDayLmer$Diversity, ymin = predict.PeakRametDayLmer$lb_PeakRametDayLmer, ymax = predict.PeakRametDayLmer$ub_PeakRametDayLmer, fill = "Peak ramet day"), alpha = 0.3) +
  
  # Peak Inflorescence 
  geom_point(data = summary_PeakInflorescence, aes(x = Diversity, y = meanDay, color = "Peak inflorescence day"), size = 2) + 
  geom_errorbar(data = summary_PeakInflorescence, aes(x = Diversity, ymin = lowerDay, ymax = upperDay, color = "Peak inflorescence day"), width = 0.2) +
  geom_line(aes(x = predict.PeakInflorescenceDayLmer$Diversity, y = predict.PeakInflorescenceDayLmer$pred_PeakInflorescenceDayLmer, color = "Peak inflorescence day"), size = 1, linetype = "dashed") +  
  geom_ribbon(aes(x = predict.PeakInflorescenceDayLmer$Diversity, ymin = predict.PeakInflorescenceDayLmer$lb_PeakInflorescenceDayLmer, ymax = predict.PeakInflorescenceDayLmer$ub_PeakInflorescenceDayLmer, fill = "Peak inflorescence day"), alpha = 0.3) +
  
  ylab('Julian day of year') +
  xlab('Genotypic diversity') +
  theme_bw() +
  theme(axis.title = element_text(size = 9)) +
  theme(axis.text = element_text(size = 9, color = "black")) +
  theme(panel.grid = element_blank()) +
  
  scale_color_manual(name = "Ontogeny", values = c("Peak inflorescence day" = "#E69F00", "Peak height day" = "#009E73", "Peak ramet day" = "#0072B2")) +
  scale_fill_manual(name = "Ontogeny", values = c("Peak inflorescence day" = "#E69F00", "Peak height day" = "#009E73", "Peak ramet day" = "#0072B2")) +
  
  scale_x_continuous(breaks = c(1, 2, 4, 8)) +


  theme(legend.position = 'bottom', 
        legend.justification = "center",
        legend.title = element_text(size = 9), 
        legend.text = element_text(size = 9))  

fig.PeakDayLmer



#  PeakInflorescenceDayLmer
PeakInflorescence.Summary <- as.data.frame(coef(summary(PeakInflorescenceDayLmer)))
PeakInflorescence.CI <- confint(PeakInflorescenceDayLmer)[-c(1, 2), ] 
PeakInflorescence.Result <- as.data.frame(cbind(PeakInflorescence.Summary, PeakInflorescence.CI))
PeakInflorescence.Result$Variable <- row.names(PeakInflorescence.Result)
PeakInflorescence.Result <- PeakInflorescence.Result %>%
  dplyr::filter(Variable != "(Intercept)") %>% 
  dplyr::select(Variable, Estimate, lower = `2.5 %`, upper = `97.5 %`) %>%
  dplyr::mutate(Model = "Peak inflorescence day")

# PeakMeanHeightLmer 
PeakHeight.Summary <- as.data.frame(coef(summary(PeakMeanHeightDayLmer)))
PeakHeight.CI <- confint(PeakMeanHeightDayLmer)[-c(1, 2), ] 
PeakHeight.Result <- as.data.frame(cbind(PeakHeight.Summary, PeakHeight.CI))
PeakHeight.Result$Variable <- row.names(PeakHeight.Result)
PeakHeight.Result <- PeakHeight.Result %>%
  dplyr::filter(Variable != "(Intercept)") %>% 
  dplyr::select(Variable, Estimate, lower = `2.5 %`, upper = `97.5 %`) %>%
  dplyr::mutate(Model = "Peak height day")


PeakRamet.Summary <- as.data.frame(coef(summary(PeakRametDayLmer)))
PeakRamet.CI <- confint(PeakRametDayLmer)[-c(1, 2), ]  
PeakRamet.Result <- as.data.frame(cbind(PeakRamet.Summary, PeakRamet.CI))
PeakRamet.Result$Variable <- row.names(PeakRamet.Result)

PeakRamet.Result <- PeakRamet.Result %>%
  dplyr::filter(Variable != "(Intercept)") %>%  
  dplyr::select(Variable, Estimate, lower = `2.5 %`, upper = `97.5 %`) %>%
  dplyr::mutate(Model = "Peak ramet day")

AllModels.Results <- bind_rows(PeakInflorescence.Result, PeakHeight.Result, PeakRamet.Result) %>%
  dplyr::mutate(CI_range = (upper - lower) / 2) %>%
  dplyr::mutate(Label = paste0(round(Estimate, 2), " ± ", round(CI_range, 2)))

Plot.AllModels.EffectSize <- ggplot(data = AllModels.Results, aes(x = Model, y = Estimate, color = Model)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, size = 1) +
  geom_point(size = 3, shape = 19) +
  geom_hline(aes(yintercept = 0), colour = "gray", linetype = "dashed") + 
  theme_bw() +
  labs(x = "", y = "Effect size (genotypic diversity)") +
  #labs(title = "Effect of genotypic diversity on the days\nof peak height, ramet and inflorescence") + 
  theme(axis.title.x = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  scale_color_manual(name = "Legend", values = c("Peak inflorescence day" = "#E69F00", 
                                                    "Peak height day" = "#009E73", 
                                                    "Peak ramet day" = "#0072B2")) +
  coord_flip()+
  theme(legend.position = "none")

Plot.AllModels.EffectSize



####combine days figure####
fig.peakDay.effectSize <- ggpubr::ggarrange(fig.PeakDayLmer,Plot.AllModels.EffectSize,ncol = 2,nrow=1,labels=c('(a)', '(b)'),font.label = list(size = 10, color = "black", face = "bold"), widths = c(0.8, 1),common.legend = TRUE, legend = "bottom")
fig.peakDay.effectSize#Figure 3

#ggsave('Result/fig.peakDay.effectSize20250915.pdf', width =18, height =10, units = 'cm')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                      ALLOMETRY OF BIOMASS ALLOCATION                     ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

names(DataFinal)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ TotalAboveBiomass~TotalBelowBiomass  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TotalAbove.TotalBelow.Lmer <- lmer(LogTotalAboveBiomass~LogTotalBelowBiomass*Diversity+(1|Composition),DataFinal)
summary(TotalAbove.TotalBelow.Lmer)
anova(TotalAbove.TotalBelow.Lmer)


#get coefficient
TotalAbove.TotalBelow.Lmer.Summary <- as.data.frame(coef(summary(TotalAbove.TotalBelow.Lmer)))
TotalAbove.TotalBelow.Lmer.CI <- confint(TotalAbove.TotalBelow.Lmer)[-c(1,2),]
TotalAbove.TotalBelow.Lmer.Result <- as.data.frame(cbind(TotalAbove.TotalBelow.Lmer.Summary,TotalAbove.TotalBelow.Lmer.CI))
TotalAbove.TotalBelow.Lmer.Result$Variable <- row.names(TotalAbove.TotalBelow.Lmer.Result)
row.names(TotalAbove.TotalBelow.Lmer.Result)
TotalAbove.TotalBelow.Lmer.Result
TotalAbove.TotalBelow.Lmer.Result <- TotalAbove.TotalBelow.Lmer.Result %>%
  dplyr::select(Variable,Estimate,lower='2.5 %',upper='97.5 %',P='Pr(>|t|)')%>%
  filter(Variable!='(Intercept)')%>%
  mutate(Allometry ='TotalAboveBiomass~TotalBelowBiomass')%>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0)))%>%
  mutate(Variable = plyr::mapvalues(Variable, c("LogTotalBelowBiomass","Diversity","LogTotalBelowBiomass:Diversity"), c("Ln(Belowground biomass)","Diversity","Ln(Belowground biomass) × Diversity")))%>%
  mutate(Colour=as.factor(Colour))
TotalAbove.TotalBelow.Lmer.Result$Variable <- factor(TotalAbove.TotalBelow.Lmer.Result$Variable,levels=c('Ln(Belowground biomass) × Diversity','Diversity','Ln(Belowground biomass)'))


#plot effect size
Plot.TotalAbove.TotalBelow.Lmer.Result <- ggplot(data=TotalAbove.TotalBelow.Lmer.Result,aes(Variable, Estimate, col=Colour)) +
  geom_errorbar(data=TotalAbove.TotalBelow.Lmer.Result, mapping=aes(ymin=lower, ymax=upper), width=0, size=1) +
  geom_point(size=2, shape=19)+##0433ff
  geom_hline(aes(yintercept=0), colour="gray", linetype="dashed")+
  theme_bw()+
  scale_colour_manual(name = NULL, values = c('#199b26','#0433ff'))  + 
  scale_fill_manual(name = NULL, values = c('#199b26','#0433ff'))

Plot.TotalAbove.TotalBelow.Lmer.Result<-Plot.TotalAbove.TotalBelow.Lmer.Result+labs(x ="")+labs(y ="Effect size")+#
  theme(axis.title.x =element_text( size=9, colour="black"),
        axis.title.y=element_text(size=9, colour="black"))+
  theme(axis.text.x =element_text(size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"))+
  theme(panel.grid = element_blank())+
  theme(plot.margin=unit(x=c(0.2,0.4,0.2,0),units="cm"))+
  theme(legend.position="none")+
  theme(strip.text = element_text(size = 7))+
  labs(title = "Ln(Aboveground biomass)")+ 
  theme(plot.title = element_text(size = 9,hjust = 1, face = "bold"))

Plot.TotalAbove.TotalBelow.Lmer.Result<-Plot.TotalAbove.TotalBelow.Lmer.Result+coord_flip()#+
Plot.TotalAbove.TotalBelow.Lmer.Result


#new data
new_LogTotalBelowBiomass <- expand.grid(
  LogTotalBelowBiomass=seq(min(DataFinal$LogTotalBelowBiomass),max(DataFinal$LogTotalBelowBiomass),length=1000),
  Diversity=c(1,2,4,8))%>%
  as.data.frame()

#predict Ramet
pred_TotalAbove.TotalBelow.Lmer<-predict(TotalAbove.TotalBelow.Lmer,newdata=new_LogTotalBelowBiomass,re.form=~0)

ci_line.TotalAbove.TotalBelow.Lmer<-bootMer(TotalAbove.TotalBelow.Lmer,FUN=function(.) predict(.,newdata=new_LogTotalBelowBiomass,
                                                             re.form=~0),nsim=1000)
ci_line.TotalAbove.TotalBelow.Lmer<-apply(ci_line.TotalAbove.TotalBelow.Lmer$t,2,function(x) x[order(x)][c(25,975)])
lb_TotalAbove.TotalBelow.Lmer<-ci_line.TotalAbove.TotalBelow.Lmer[1,]
ub_TotalAbove.TotalBelow.Lmer<-ci_line.TotalAbove.TotalBelow.Lmer[2,]

predict.TotalAbove.TotalBelow.Lmer <-cbind(new_LogTotalBelowBiomass,pred_TotalAbove.TotalBelow.Lmer,lb_TotalAbove.TotalBelow.Lmer,ub_TotalAbove.TotalBelow.Lmer)%>%
  mutate(TotalBelowBiomass=exp(LogTotalBelowBiomass),TotalAboveBiomass=exp(pred_TotalAbove.TotalBelow.Lmer),lb_TotalAboveBiomass=exp(lb_TotalAbove.TotalBelow.Lmer), ub_TotalAboveBiomass=exp(ub_TotalAbove.TotalBelow.Lmer))




fig.log.TotalAbove.TotalBelow.Lmer <- ggplot(predict.TotalAbove.TotalBelow.Lmer)+
  geom_point(
    data = DataFinal,
    aes(x = LogTotalBelowBiomass, y = LogTotalAboveBiomass, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.1),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=LogTotalBelowBiomass,y=pred_TotalAbove.TotalBelow.Lmer,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=LogTotalBelowBiomass,ymin=lb_TotalAbove.TotalBelow.Lmer,ymax=ub_TotalAbove.TotalBelow.Lmer,fill=as.factor(Diversity)),alpha=0.3)+
  ylab('Ln(Aboveground biomass (g)) ')+
  xlab('Ln(Belowground biomass (g))')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  #scale_linetype_discrete(name = "Treatment")+
  theme(legend.position = c(0.8,0.25),legend.box = 'horizontal',legend.key.size = unit(0.5,'cm'),legend.box.spacing = unit(0.4,'cm'),legend.key.height = unit(0.3,'cm'),legend.box.margin = margin(0,0,0,0,'cm'),legend.box.just = 'bottom',legend.box.background = element_blank(),legend.spacing = unit(0,'cm'), legend.title = element_text(size=8),legend.text = element_text(size=8))
fig.log.TotalAbove.TotalBelow.Lmer



fig.TotalAbove.TotalBelow.Lmer <- ggplot(predict.TotalAbove.TotalBelow.Lmer)+
  geom_point(
    data = DataFinal,
    aes(x = TotalBelowBiomass, y = TotalAboveBiomass, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.5),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=TotalBelowBiomass,y=TotalAboveBiomass,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=TotalBelowBiomass,ymin=lb_TotalAboveBiomass,ymax=ub_TotalAboveBiomass,fill=as.factor(Diversity)),alpha=0.2)+
  ylab('Aboveground biomass (g)')+
  xlab('Belowground biomass (g)')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  theme(legend.position = c(0.8,0.25),legend.box = 'horizontal',legend.key.size = unit(0.5,'cm'),legend.box.spacing = unit(0.4,'cm'),legend.key.height = unit(0.3,'cm'),legend.box.margin = margin(0,0,0,0,'cm'),legend.box.just = 'bottom',legend.box.background = element_blank(),legend.spacing = unit(0,'cm'), legend.title = element_text(size=8),legend.text = element_text(size=8))

fig.TotalAbove.TotalBelow.Lmer


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ ReproductiveBiomass~VegetativeBiomass  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names(DataFinal)
Reproductive.Vegetative.Lmer <- lmer(LogReproductiveBiomass~LogVegetativeBiomass*Diversity+(1|Composition),DataFinal)
summary(Reproductive.Vegetative.Lmer)
anova(Reproductive.Vegetative.Lmer)


#get coefficient
Reproductive.Vegetative.Lmer.Summary <- as.data.frame(coef(summary(Reproductive.Vegetative.Lmer)))
Reproductive.Vegetative.Lmer.CI <- confint(Reproductive.Vegetative.Lmer)[-c(1,2),]
Reproductive.Vegetative.Lmer.Result <- as.data.frame(cbind(Reproductive.Vegetative.Lmer.Summary,Reproductive.Vegetative.Lmer.CI))
Reproductive.Vegetative.Lmer.Result$Variable <- row.names(Reproductive.Vegetative.Lmer.Result)
row.names(Reproductive.Vegetative.Lmer.Result)
Reproductive.Vegetative.Lmer.Result
Reproductive.Vegetative.Lmer.Result <- Reproductive.Vegetative.Lmer.Result %>%
  dplyr::select(Variable,Estimate,lower='2.5 %',upper='97.5 %',P='Pr(>|t|)')%>%
  filter(Variable!='(Intercept)')%>%
  mutate(Allometry ='ReproductiveBiomass~VegetativeBiomass')%>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0)))%>%
  mutate(Variable = plyr::mapvalues(Variable, c("LogVegetativeBiomass","Diversity","LogVegetativeBiomass:Diversity"), c("Ln(Vegetative biomass)","Diversity","Ln(==Vegetative biomass) × Diversity")))%>%
  mutate(Colour=as.factor(Colour))
Reproductive.Vegetative.Lmer.Result$Variable <- factor(Reproductive.Vegetative.Lmer.Result$Variable,levels=c('Ln(==Vegetative biomass) × Diversity','Diversity','Ln(Vegetative biomass)'))


#plot effect size
Plot.Reproductive.Vegetative.Lmer.Result <- ggplot(data=Reproductive.Vegetative.Lmer.Result,aes(Variable, Estimate, col=Colour)) +
  geom_errorbar(data=Reproductive.Vegetative.Lmer.Result, mapping=aes(ymin=lower, ymax=upper), width=0, size=1) +
  geom_point(size=2, shape=19)+##0433ff
  geom_hline(aes(yintercept=0), colour="gray", linetype="dashed")+
  theme_bw()+
  scale_colour_manual(name = NULL, values = c('gray','#0433ff'))  + 
  scale_fill_manual(name = NULL, values = c('gray','#0433ff'))

Plot.Reproductive.Vegetative.Lmer.Result<-Plot.Reproductive.Vegetative.Lmer.Result+labs(x ="")+labs(y ="Effect size")+#
  theme(axis.title.x =element_text( size=9, colour="black"),
        axis.title.y=element_text(size=9, colour="black"))+
  theme(axis.text.x =element_text(size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"))+
  theme(panel.grid = element_blank())+
  theme(plot.margin=unit(x=c(0.2,0.4,0.2,0),units="cm"))+
  theme(legend.position="none")+
  theme(strip.text = element_text(size = 7))+
  labs(title = "Ln(Reproductive biomass)")+ 
  theme(plot.title = element_text(size = 9,hjust = 1, face = "bold"))

Plot.Reproductive.Vegetative.Lmer.Result<-Plot.Reproductive.Vegetative.Lmer.Result+coord_flip()#+
Plot.Reproductive.Vegetative.Lmer.Result


#new data
new_LogVegetativeBiomass <- expand.grid(
  LogVegetativeBiomass=seq(min(DataFinal$LogVegetativeBiomass),max(DataFinal$LogVegetativeBiomass),length=1000),
  Diversity=c(1,2,4,8))%>%
  as.data.frame()

#predict Ramet
pred_Reproductive.Vegetative.Lmer<-predict(Reproductive.Vegetative.Lmer,newdata=new_LogVegetativeBiomass,re.form=~0)

ci_line.Reproductive.Vegetative.Lmer<-bootMer(Reproductive.Vegetative.Lmer,FUN=function(.) predict(.,newdata=new_LogVegetativeBiomass,
                                                                                               re.form=~0),nsim=1000)
ci_line.Reproductive.Vegetative.Lmer<-apply(ci_line.Reproductive.Vegetative.Lmer$t,2,function(x) x[order(x)][c(25,975)])
lb_Reproductive.Vegetative.Lmer<-ci_line.Reproductive.Vegetative.Lmer[1,]
ub_Reproductive.Vegetative.Lmer<-ci_line.Reproductive.Vegetative.Lmer[2,]

predict.Reproductive.Vegetative.Lmer <-cbind(new_LogVegetativeBiomass,pred_Reproductive.Vegetative.Lmer,lb_Reproductive.Vegetative.Lmer,ub_Reproductive.Vegetative.Lmer)%>%
  mutate(VegetativeBiomass=exp(LogVegetativeBiomass),ReproductiveBiomass=exp(pred_Reproductive.Vegetative.Lmer),lb_ReproductiveBiomass=exp(lb_Reproductive.Vegetative.Lmer), ub_ReproductiveBiomass=exp(ub_Reproductive.Vegetative.Lmer))



fig.log.Reproductive.Vegetative.Lmer <- ggplot(predict.Reproductive.Vegetative.Lmer)+
  geom_point(
    data = DataFinal,
    aes(x = LogVegetativeBiomass, y = LogReproductiveBiomass, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.1),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=LogVegetativeBiomass,y=pred_Reproductive.Vegetative.Lmer,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=LogVegetativeBiomass,ymin=lb_Reproductive.Vegetative.Lmer,ymax=ub_Reproductive.Vegetative.Lmer,fill=as.factor(Diversity)),alpha=0.3)+
  ylab('Ln(Reproductive biomass (g)) ')+
  xlab('Ln(Vegetative biomass (g))')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  #scale_linetype_discrete(name = "Treatment")+
  theme(legend.position = c(0.8,0.25),legend.box = 'horizontal',legend.key.size = unit(0.5,'cm'),legend.box.spacing = unit(0.4,'cm'),legend.key.height = unit(0.3,'cm'),legend.box.margin = margin(0,0,0,0,'cm'),legend.box.just = 'bottom',legend.box.background = element_blank(),legend.spacing = unit(0,'cm'), legend.title = element_text(size=8),legend.text = element_text(size=8))
fig.log.Reproductive.Vegetative.Lmer



fig.Reproductive.Vegetative.Lmer <- ggplot(predict.Reproductive.Vegetative.Lmer)+
  geom_point(
    data = DataFinal,
    aes(x = VegetativeBiomass, y = ReproductiveBiomass, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.5),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=VegetativeBiomass,y=ReproductiveBiomass,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=VegetativeBiomass,ymin=lb_ReproductiveBiomass,ymax=ub_ReproductiveBiomass,fill=as.factor(Diversity)),alpha=0.2)+
  ylab('Reproductive biomass (g)')+
  xlab('Vegetative biomass (g)')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  theme(legend.position = c(0.8,0.25),legend.box = 'horizontal',legend.key.size = unit(0.5,'cm'),legend.box.spacing = unit(0.4,'cm'),legend.key.height = unit(0.3,'cm'),legend.box.margin = margin(0,0,0,0,'cm'),legend.box.just = 'bottom',legend.box.background = element_blank(),legend.spacing = unit(0,'cm'), legend.title = element_text(size=8),legend.text = element_text(size=8))
fig.Reproductive.Vegetative.Lmer



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ InflorescenceBiomass~AsexualBiomass  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names(DataFinal)
Inflorescence.Asexual.Lmer <- lmer(LogInflorescenceBiomass~LogAsexualBiomass*Diversity+(1|Composition),DataFinal)
summary(Inflorescence.Asexual.Lmer)
anova(Inflorescence.Asexual.Lmer)


#get coefficient
Inflorescence.Asexual.Lmer.Summary <- as.data.frame(coef(summary(Inflorescence.Asexual.Lmer)))
Inflorescence.Asexual.Lmer.CI <- confint(Inflorescence.Asexual.Lmer)[-c(1,2),]
Inflorescence.Asexual.Lmer.Result <- as.data.frame(cbind(Inflorescence.Asexual.Lmer.Summary,Inflorescence.Asexual.Lmer.CI))
Inflorescence.Asexual.Lmer.Result$Variable <- row.names(Inflorescence.Asexual.Lmer.Result)
row.names(Inflorescence.Asexual.Lmer.Result)
Inflorescence.Asexual.Lmer.Result
Inflorescence.Asexual.Lmer.Result <- Inflorescence.Asexual.Lmer.Result %>%
  dplyr::select(Variable,Estimate,lower='2.5 %',upper='97.5 %',P='Pr(>|t|)')%>%
  filter(Variable!='(Intercept)')%>%
  mutate(Allometry ='InflorescenceBiomass~AsexualBiomass')%>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0)))%>%
  mutate(Variable = plyr::mapvalues(Variable, c("LogAsexualBiomass","Diversity","LogAsexualBiomass:Diversity"), c("Ln(Asexual reproduction)","Diversity","Ln(=Asexual reproduction) × Diversity")))%>%
  mutate(Colour=as.factor(Colour))
Inflorescence.Asexual.Lmer.Result$Variable <- factor(Inflorescence.Asexual.Lmer.Result$Variable,levels=c('Ln(=Asexual reproduction) × Diversity','Diversity','Ln(Asexual reproduction)'))


#plot effect size
Plot.Inflorescence.Asexual.Lmer.Result <- ggplot(data=Inflorescence.Asexual.Lmer.Result,aes(Variable, Estimate, col=Colour)) +
  geom_errorbar(data=Inflorescence.Asexual.Lmer.Result, mapping=aes(ymin=lower, ymax=upper), width=0, size=1) +
  geom_point(size=2, shape=19)+##0433ff
  geom_hline(aes(yintercept=0), colour="gray", linetype="dashed")+
  theme_bw()+
  scale_colour_manual(name = NULL, values = c('gray','#0433ff'))  + 
  scale_fill_manual(name = NULL, values = c('gray','#0433ff'))

Plot.Inflorescence.Asexual.Lmer.Result<-Plot.Inflorescence.Asexual.Lmer.Result+labs(x ="")+labs(y ="Effect size")+#
  theme(axis.title.x =element_text( size=9, colour="black"),
        axis.title.y=element_text(size=9, colour="black"))+
  theme(axis.text.x =element_text(size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"))+
  theme(panel.grid = element_blank())+
  theme(plot.margin=unit(x=c(0.2,0.4,0.2,0),units="cm"))+
  theme(legend.position="none")+
  theme(strip.text = element_text(size = 7))+
  labs(title = "Ln(Inflorescence biomass)")+ 
  theme(plot.title = element_text(size = 9,hjust = 1, face = "bold"))

Plot.Inflorescence.Asexual.Lmer.Result<-Plot.Inflorescence.Asexual.Lmer.Result+coord_flip()#+

Plot.Inflorescence.Asexual.Lmer.Result


#new data
new_LogAsexualBiomass <- expand.grid(
  LogAsexualBiomass=seq(min(DataFinal$LogAsexualBiomass),max(DataFinal$LogAsexualBiomass),length=1000),
  Diversity=c(1,2,4,8))%>%
  as.data.frame()

#predict Ramet
pred_Inflorescence.Asexual.Lmer<-predict(Inflorescence.Asexual.Lmer,newdata=new_LogAsexualBiomass,re.form=~0)

ci_line.Inflorescence.Asexual.Lmer<-bootMer(Inflorescence.Asexual.Lmer,FUN=function(.) predict(.,newdata=new_LogAsexualBiomass,
                                                                                                   re.form=~0),nsim=1000)
ci_line.Inflorescence.Asexual.Lmer<-apply(ci_line.Inflorescence.Asexual.Lmer$t,2,function(x) x[order(x)][c(25,975)])
lb_Inflorescence.Asexual.Lmer<-ci_line.Inflorescence.Asexual.Lmer[1,]
ub_Inflorescence.Asexual.Lmer<-ci_line.Inflorescence.Asexual.Lmer[2,]

predict.Inflorescence.Asexual.Lmer <-cbind(new_LogAsexualBiomass,pred_Inflorescence.Asexual.Lmer,lb_Inflorescence.Asexual.Lmer,ub_Inflorescence.Asexual.Lmer)%>%
  mutate(AsexualBiomass=exp(LogAsexualBiomass),InflorescenceBiomass=exp(pred_Inflorescence.Asexual.Lmer),lb_InflorescenceBiomass=exp(lb_Inflorescence.Asexual.Lmer), ub_InflorescenceBiomass=exp(ub_Inflorescence.Asexual.Lmer))




fig.log.Inflorescence.Asexual.Lmer <- ggplot(predict.Inflorescence.Asexual.Lmer)+
  geom_point(
    data = DataFinal,
    aes(x = LogAsexualBiomass, y = LogInflorescenceBiomass, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.1),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=LogAsexualBiomass,y=pred_Inflorescence.Asexual.Lmer,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=LogAsexualBiomass,ymin=lb_Inflorescence.Asexual.Lmer,ymax=ub_Inflorescence.Asexual.Lmer,fill=as.factor(Diversity)),alpha=0.3)+
  ylab('Ln(Inflorescence biomass (g)) ')+
  xlab('Ln(Asexual reproduction (g))')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  #scale_linetype_discrete(name = "Treatment")+
  theme(legend.position = c(0.8,0.77),legend.box = 'horizontal',legend.key.size = unit(0.5,'cm'),legend.box.spacing = unit(0.4,'cm'),legend.key.height = unit(0.3,'cm'),legend.box.margin = margin(0,0,0,0,'cm'),legend.box.just = 'bottom',legend.box.background = element_blank(),legend.spacing = unit(0,'cm'), legend.title = element_text(size=8),legend.text = element_text(size=8))
fig.log.Inflorescence.Asexual.Lmer



fig.Inflorescence.Asexual.Lmer <- ggplot(predict.Inflorescence.Asexual.Lmer)+
  geom_point(
    data = DataFinal,
    aes(x = AsexualBiomass, y = InflorescenceBiomass, col = as.factor(Diversity)),
    alpha = 0.3, size = 0.5,
    position = position_jitter(width = 0, height = 0.5),
    inherit.aes = FALSE
  ) +
  geom_line(aes(x=AsexualBiomass,y=InflorescenceBiomass,col=as.factor(Diversity)),size=1)+
  geom_ribbon(aes(x=AsexualBiomass,ymin=lb_InflorescenceBiomass,ymax=ub_InflorescenceBiomass,fill=as.factor(Diversity)),alpha=0.2)+
  ylab('Inflorescence biomass (g)')+
  xlab('Asexual reproduction (g)')+
  theme_bw()+
  theme(axis.title = element_text(size = 9))+
  theme(axis.text = element_text(size =9, color = "black"))+
  theme(panel.grid = element_blank())+
  scale_colour_manual(name = 'Diversity', values = SpeciesCol)  + 
  scale_fill_manual(name = 'Diversity', values = SpeciesCol)+
  theme(legend.background = element_blank())+
  #scale_linetype_discrete(name = "Treatment")+
  theme(legend.position = c(0.8,0.77),legend.box = 'horizontal',legend.key.size = unit(0.5,'cm'),legend.box.spacing = unit(0.4,'cm'),legend.key.height = unit(0.3,'cm'),legend.box.margin = margin(0,0,0,0,'cm'),legend.box.just = 'bottom',legend.box.background = element_blank(),legend.spacing = unit(0,'cm'), legend.title = element_text(size=8),legend.text = element_text(size=8))
fig.Inflorescence.Asexual.Lmer




##~~~~~~~~~~~~~~~~~~~~~~
##  ~ combine plot  ----
##~~~~~~~~~~~~~~~~~~~~~~

fig.allometry.TotalBiomass <- ggpubr::ggarrange(ggpubr::ggarrange(Plot.TotalAbove.TotalBelow.Lmer.Result, fig.log.TotalAbove.TotalBelow.Lmer + theme(legend.position="none"), fig.TotalAbove.TotalBelow.Lmer + theme(legend.position="none"), Plot.Reproductive.Vegetative.Lmer.Result, fig.log.Reproductive.Vegetative.Lmer + theme(legend.position="none"), fig.Reproductive.Vegetative.Lmer + theme(legend.position="none"), Plot.Inflorescence.Asexual.Lmer.Result, fig.log.Inflorescence.Asexual.Lmer + theme(legend.position="none"), fig.Inflorescence.Asexual.Lmer + theme(legend.position="none"), ncol=3, nrow=3, labels=c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'), font.label=list(size=10, color="black", face="bold"), widths=c(1.6,0.9,0.9), common.legend=FALSE), cowplot::get_legend(fig.log.TotalAbove.TotalBelow.Lmer + guides(colour=guide_legend(nrow=1, ncol=4, byrow=TRUE), fill=guide_legend(nrow=1, ncol=4, byrow=TRUE)) + theme(legend.position="bottom", legend.box="horizontal", legend.margin=margin(2,0,2,0,"pt"), legend.key.width=unit(10,"pt"), legend.key.height=unit(10,"pt"), legend.text=element_text(size=8), legend.title=element_text(size=8))), ncol=1, heights=c(1,0.08))

#Figure S5

#ggsave('Result/fig.allometry.TotalBiomass20250916.pdf',width = 18,height =18,units = c('cm'))



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                              SMA OF EACH PLOT                            ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ sma: Ramet~MeanHeight of each plot  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(MonthlyData)
names(MonthlyData)
sma.Ramet <- sma(Ramet~MeanHeight*Plot,data=MonthlyData,log='xy')

summary(sma.Ramet)
smaResult.Ramet <- coef(sma.Ramet)
smaResult.Ramet$Plot <- row.names(smaResult.Ramet)

smaResult.Ramet <-  smaResult.Ramet%>%merge(PlotTreatment)
str(smaResult.Ramet)
smaResult.Ramet <- smaResult.Ramet%>%
  rename(slope.Ramet='slope',elevation.Ramet='elevation')
  
cor(smaResult.Ramet$slope,smaResult.Ramet$elevation)
slope.Ramet.lm <- lm(slope.Ramet~Diversity,smaResult.Ramet)
summary(slope.Ramet.lm)
hist(smaResult.Ramet$slope.Ramet)

elevation.Ramet.lm <- lm(elevation.Ramet~Diversity,smaResult.Ramet)
summary(elevation.Ramet.lm)
hist(smaResult.Ramet$elevation.Ramet)


pval.Ramet <- as.data.frame(sma.Ramet$pval)%>%
  gather(key = Plot,value = pval)%>%
  mutate(Plot = str_sub(Plot, 2, -1))%>%
  mutate(Plot = str_replace_all(Plot, "\\.", "-"))

r2.Ramet<- as.data.frame(sma.Ramet$r2)%>%
  gather(key = Plot,value = r2)%>%
  mutate(Plot = str_sub(Plot, 2, -1))%>%
  mutate(Plot = str_replace_all(Plot, "\\.", "-"))

number.Ramet<- as.data.frame(sma.Ramet$n)%>%
  gather(key = Plot,value = number)%>%
  mutate(Plot = str_sub(Plot, 2, -1))%>%
  mutate(Plot = str_replace_all(Plot, "\\.", "-"))

sma.Ramet.stat <- merge(pval.Ramet,smaResult.Ramet,by="Plot")%>%
  merge(number.Ramet,by="Plot")%>%
  merge(r2.Ramet,by="Plot")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ sma: Inflorescence~MeanHeight of each plot  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sma.Inflorescence <- sma(Inflorescence+1~MeanHeight*Plot,data=MonthlyData,log='xy')#GR
summary(sma.Inflorescence)
smaResult.Inflorescence <- coef(sma.Inflorescence)
smaResult.Inflorescence$Plot <- row.names(smaResult.Inflorescence)

smaResult.Inflorescence <-  smaResult.Inflorescence%>%merge(PlotTreatment)
str(smaResult.Inflorescence)
smaResult.Inflorescence <- smaResult.Inflorescence%>%
  rename(slope.Inflorescence='slope',elevation.Inflorescence='elevation')

cor(smaResult.Inflorescence$slope,smaResult.Inflorescence$elevation)
slope.Inflorescence.lm <- lm(slope.Inflorescence~Diversity,smaResult.Inflorescence)
summary(slope.Inflorescence.lm)
hist(smaResult.Inflorescence$slope.Inflorescence)

elevation.Inflorescence.lm <- lm(elevation.Inflorescence~Diversity,smaResult.Inflorescence)
summary(elevation.Inflorescence.lm)
hist(smaResult.Inflorescence$elevation.Inflorescence)


pval.Inflorescence <- as.data.frame(sma.Inflorescence$pval)%>%
  gather(key = Plot,value = pval)%>%
  mutate(Plot = str_sub(Plot, 2, -1))%>%
  mutate(Plot = str_replace_all(Plot, "\\.", "-"))

r2.Inflorescence<- as.data.frame(sma.Inflorescence$r2)%>%
  gather(key = Plot,value = r2)%>%
  mutate(Plot = str_sub(Plot, 2, -1))%>%
  mutate(Plot = str_replace_all(Plot, "\\.", "-"))

number.Inflorescence<- as.data.frame(sma.Inflorescence$n)%>%
  gather(key = Plot,value = number)%>%
  mutate(Plot = str_sub(Plot, 2, -1))%>%
  mutate(Plot = str_replace_all(Plot, "\\.", "-"))

sma.Inflorescence.stat <- merge(pval.Inflorescence,smaResult.Inflorescence,by="Plot")%>%
  merge(number.Inflorescence,by="Plot")%>%
  merge(r2.Inflorescence,by="Plot")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ sma: Inflorescence~Ramet of each plot  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sma.Inflorescence.Ramet <- sma(Inflorescence+1~Ramet*Plot,data=MonthlyData,log='xy')#GR
summary(sma.Inflorescence.Ramet)
smaResult.Inflorescence.Ramet <- coef(sma.Inflorescence.Ramet)
smaResult.Inflorescence.Ramet$Plot <- row.names(smaResult.Inflorescence.Ramet)

smaResult.Inflorescence.Ramet <-  smaResult.Inflorescence.Ramet%>%merge(PlotTreatment)
str(smaResult.Inflorescence.Ramet)
smaResult.Inflorescence.Ramet <- smaResult.Inflorescence.Ramet%>%
  rename(slope.Inflorescence.Ramet='slope',elevation.Inflorescence.Ramet='elevation')

cor(smaResult.Inflorescence.Ramet$slope.Inflorescence.Ramet,smaResult.Inflorescence.Ramet$elevation.Inflorescence.Ramet)
slope.Inflorescence.Ramet.lm <- lm(slope.Inflorescence.Ramet~Diversity,smaResult.Inflorescence.Ramet)
summary(slope.Inflorescence.Ramet.lm)
hist(smaResult.Inflorescence.Ramet$slope.Inflorescence.Ramet)

elevation.Inflorescence.Ramet.lm <- lm(elevation.Inflorescence.Ramet~Diversity,smaResult.Inflorescence.Ramet)
summary(elevation.Inflorescence.Ramet.lm)
hist(smaResult.Inflorescence.Ramet$elevation.Inflorescence.Ramet)


pval.Inflorescence.Ramet <- as.data.frame(sma.Inflorescence.Ramet$pval)%>%
  gather(key = Plot,value = pval)%>%
  mutate(Plot = str_sub(Plot, 2, -1))%>%
  mutate(Plot = str_replace_all(Plot, "\\.", "-"))

r2.Inflorescence.Ramet<- as.data.frame(sma.Inflorescence.Ramet$r2)%>%
  gather(key = Plot,value = r2)%>%
  mutate(Plot = str_sub(Plot, 2, -1))%>%
  mutate(Plot = str_replace_all(Plot, "\\.", "-"))

number.Inflorescence.Ramet<- as.data.frame(sma.Inflorescence.Ramet$n)%>%
  gather(key = Plot,value = number)%>%
  mutate(Plot = str_sub(Plot, 2, -1))%>%
  mutate(Plot = str_replace_all(Plot, "\\.", "-"))

sma.Inflorescence.Ramet.stat <- merge(pval.Inflorescence.Ramet,smaResult.Inflorescence.Ramet,by="Plot")%>%
  merge(number.Inflorescence.Ramet,by="Plot")%>%
  merge(r2.Inflorescence.Ramet,by="Plot")



smaResult <- smaResult.Inflorescence%>%
  merge(smaResult.Ramet)%>%
  merge(smaResult.Inflorescence.Ramet)




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                    THE CORRELATION AMONG SOIL                  ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


names(Soil)


cor.SoilData<- Soil%>%
  dplyr::select(-Plot,-Diversity,-Composition,-Repetition)
names(cor.SoilData)
colnames(cor.SoilData) <- c("Soil organic carbon" ,"Total nitrogen", "Ammonium","Nitrate","Available phosphorus" , "Available potassium", "Microbial biomass carbon" , "Bacterial Chao1" ,"Bacterial diversity","Bacterial PCoA1" ,"Fungal Chao1" ,"Fungal diversity","Fungal PCoA1")


str(cor.SoilData)
M <- cor(cor.SoilData)


res1 <- cor.mtest(M, conf.level = 0.95)
res2 <- cor.mtest(M, conf.level = 0.99)

corrplot::corrplot(M,type = "full",p.mat = res1$p,method='pie', order = "hclust",addrect  = 4,hclust.method="ward.D2",
                   rect.col = "red",tl.srt = 30, cl.cex =1.3,
                   tl.cex =1.2,tl.col='black',tl.pos = "lt",cl.pos = "b")
corrplot::corrplot(M,add = TRUE, type = "lower", method = "number",number.font = 2, order = "hclust",
                   tl.cex = 0.8,
                   addrect = 2,hclust.method="ward.D2",
                   diag = F,tl.pos = "n", cl.pos = "n",cl.cex =1)
#save by hand ,Soil.CorPlot.pdf
#Figure S1




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                              COMBINE PLOT DATA                           ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PlotData <- DataFinal%>%
  merge(smaResult)%>%
  mutate(BiomassRation=TotalAboveBiomass/TotalBelowBiomass)%>%
  merge(Soil)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                     THE INFLUENCE OF SOIL ON ALLOMETRY                   ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

names(Soil)


PlotDataScale <- PlotData%>%
  mutate_at(c('Diversity','SOC','TN','NH4','NO3','AP','AK','MBC','B_chao','B_shannon','B_PC1','F_chao','F_shannon','F_PC1','slope.Ramet','slope.Inflorescence','slope.Inflorescence.Ramet','RametNumber'), ~(scale(.) %>% as.vector))

####slope.Ramet####

detach("package:lmerTest", unload=TRUE) 
slope.Ramet.soil.lmer <-lmer(slope.Ramet~Diversity + SOC + TN + NH4 + NO3 + AP + AK + MBC + B_chao + B_shannon + B_PC1 + F_chao + F_shannon + F_PC1+(1|Composition),PlotDataScale, na.action = "na.fail")
summary(slope.Ramet.soil.lmer)
stepba(slope.Ramet.soil.lmer)
library(lmerTest)
slope.Ramet.soil.final <-lmer(slope.Ramet~NO3 + AP + F_shannon+(1|Composition),PlotDataScale, na.action = "na.fail")


summary(slope.Ramet.soil.final)


slope.Ramet.soil.Summary <- as.data.frame(coef(summary(slope.Ramet.soil.final)))
slope.Ramet.soil.CI <- confint(slope.Ramet.soil.final)[-c(1,2),]
slope.Ramet.soil <- as.data.frame(cbind(slope.Ramet.soil.Summary,slope.Ramet.soil.CI))
slope.Ramet.soil$Variable <- row.names(slope.Ramet.soil)

row.names(slope.Ramet.soil)
slope.Ramet.soil
slope.Ramet.soil <- slope.Ramet.soil %>%
  dplyr::select(Variable,Estimate,lower='2.5 %',upper='97.5 %',P='Pr(>|t|)')%>%
  filter(Variable!='(Intercept)')%>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0)))%>%
  mutate(Title='The scaling exponent')
slope.Ramet.soil



####intercept.Ramet####

detach("package:lmerTest", unload=TRUE) 
elevation.Ramet.soil.lmer <-lmer(elevation.Ramet~Diversity + SOC + TN + NH4 + NO3 + AP + AK + MBC + B_chao + B_shannon + B_PC1 + F_chao + F_shannon + F_PC1+(1|Composition),PlotDataScale, na.action = "na.fail")
summary(elevation.Ramet.soil.lmer)
stepba(elevation.Ramet.soil.lmer)
library(lmerTest)
elevation.Ramet.soil.final <-lmer(elevation.Ramet~NO3 + AP + F_shannon+(1|Composition),PlotDataScale, na.action = "na.fail")

summary(elevation.Ramet.soil.final)



elevation.Ramet.soil.Summary <- as.data.frame(coef(summary(elevation.Ramet.soil.final)))
elevation.Ramet.soil.CI <- confint(elevation.Ramet.soil.final)[-c(1,2),]
elevation.Ramet.soil <- as.data.frame(cbind(elevation.Ramet.soil.Summary,elevation.Ramet.soil.CI))
elevation.Ramet.soil$Variable <- row.names(elevation.Ramet.soil)

row.names(elevation.Ramet.soil)
elevation.Ramet.soil
elevation.Ramet.soil <- elevation.Ramet.soil %>%
  dplyr::select(Variable,Estimate,lower='2.5 %',upper='97.5 %',P='Pr(>|t|)')%>%
  filter(Variable!='(Intercept)')%>%
  mutate(Colour = ifelse(lower > 0 & upper > 0, 1,
                         ifelse(lower < 0 & upper < 0, -1, 0)))%>%
  mutate(Title='The intercept')
elevation.Ramet.soil



####Plot####
unique(Allometry.Soil$Variable)
#combine data
Allometry.Soil <-  slope.Ramet.soil%>%
  rbind(elevation.Ramet.soil)%>%
  #rbind(slope.Inflorescence.Ramet.soil)%>%
  mutate(Variable = plyr::mapvalues(Variable, c("NO3","AP","F_shannon" ), c("Nitrate","Available phosphorus","Fungal diversity" )))

Allometry.Soil$Variable <- factor(Allometry.Soil$Variable ,levels = c("Nitrate","Available phosphorus","Fungal diversity" ))  

Allometry.Soil$Colour <- factor(Allometry.Soil$Colour,levels=c('-1','1'))
Allometry.Soil$Title <- factor(Allometry.Soil$Title,levels=c('The scaling exponent','The intercept'))

#plot
colors2 = c('#0433ff','gray','#199b26')


Plot.Allometry.Soil <- ggplot(data=Allometry.Soil,aes(Variable, Estimate, col=Colour)) +
  facet_wrap(~ Title,  nrow = 1,scales = "free_x") + #scales = "free_y",
  geom_errorbar(data=Allometry.Soil, mapping=aes(ymin=lower, ymax=upper), width=0, size=1) +
  geom_point(size=2, shape=19)+##0433ff
  geom_hline(aes(yintercept=0), colour="gray", linetype="dashed")+
  theme_bw()+
  scale_colour_manual(name = NULL, values = c('#199b26','#0433ff'))  + 
  scale_fill_manual(name = NULL, values = c('#199b26','#0433ff'))

Plot.Allometry.Soil<-Plot.Allometry.Soil+labs(x ="")+labs(y ="Effect size")+
  theme(axis.title.x =element_text( size=9, colour="black"),
        axis.title.y=element_text(size=9, colour="black"))+
  theme(axis.text.x =element_text(size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"))+
  theme(panel.grid = element_blank())+
  theme(plot.margin=unit(x=c(0.2,0.4,0.2,0),units="cm"))+
  theme(legend.position="none")+
  theme(strip.text = element_text(size = 8))+
  theme(panel.spacing = unit(0.5, "cm")) 

Plot.Allometry.Soil<-Plot.Allometry.Soil+coord_flip()#+


Plot.Allometry.Soil#Figure S3

#ggsave('Result/Plot.Allometry.Soil20250113.pdf', width =10, height =5, units = 'cm')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                                    SEM                                   ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#prepare SEM DATA
names(PlotData)
library(lmerTest)
detach("package:lmerTest", unload=TRUE) 

PlotData$Composition<- as.factor(PlotData$Composition)
summary(PlotData)
str(PlotData)
names(PlotData)
summary(Soil)

sem1<-psem(
  lmer(LogInflorescenceNumber~LogRametNumber+(1|Composition),  data=PlotData),
  lmer(LogCormNumber~LogRametNumber+(1|Composition), data=PlotData),
  lmer(LogRametNumber~slope.Ramet+elevation.Ramet +(1|Composition), data=PlotData),
  lmer(slope.Ramet~F_shannon+AP+NO3+(1|Composition), data=PlotData),
  lmer(elevation.Ramet~F_shannon+AP+NO3+(1|Composition), data=PlotData),
  lmer(F_shannon~Diversity+(1|Composition), data=PlotData),
  lmer(AP~Diversity+(1|Composition), data=PlotData),
  lmer(NO3~Diversity+(1|Composition), data=PlotData),
  slope.Ramet%~~%elevation.Ramet,
  data=PlotData
)
summary(sem1)


sem2<-psem(
  lm(LogInflorescenceNumber~LogRametNumber,  data=PlotData),
  lmer(LogCormNumber~LogRametNumber+(1|Composition), data=PlotData),
  lmer(LogRametNumber~slope.Ramet+elevation.Ramet +(1|Composition), data=PlotData),
  lmer(slope.Ramet~F_shannon+AP+NO3+(1|Composition), data=PlotData),
  lmer(elevation.Ramet~F_shannon+AP+NO3+(1|Composition), data=PlotData),
  lm(F_shannon~Diversity, data=PlotData),
  lmer(AP~Diversity+(1|Composition), data=PlotData),
  lmer(NO3~Diversity+(1|Composition), data=PlotData),
  slope.Ramet%~~%elevation.Ramet,
  data=PlotData
)
summary(sem2)


sem3<-psem(
  lm(LogInflorescenceNumber~LogRametNumber,  data=PlotData),
  lmer(LogCormNumber~LogRametNumber+(1|Composition), data=PlotData),
  lmer(LogRametNumber~slope.Ramet+elevation.Ramet +Diversity+(1|Composition), data=PlotData),
  lmer(slope.Ramet~F_shannon+AP+NO3+(1|Composition), data=PlotData),
  lmer(elevation.Ramet~F_shannon+AP+NO3+(1|Composition), data=PlotData),
  lm(F_shannon~Diversity, data=PlotData),
  lmer(AP~Diversity+(1|Composition), data=PlotData),
  lmer(NO3~Diversity+(1|Composition), data=PlotData),
  slope.Ramet%~~%elevation.Ramet,
  data=PlotData
)
summary(sem3)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                            FIG: CHANGE WITH DAYS                         ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mean
daily_means <- MonthlyData %>%
  group_by(Day, Diversity) %>%
  summarise(
    MeanHeight = mean(MeanHeight, na.rm = TRUE),
    Ramet = mean(Ramet, na.rm = TRUE),
    Inflorescence = mean(Inflorescence, na.rm = TRUE)
  )

str(MonthlyData)

MonthlyData$Diversity <- factor(MonthlyData$Diversity)
daily_means$Diversity <- factor(daily_means$Diversity)

# Mean Height vs Day
fig.MeanHeight <- ggplot(MonthlyData, aes(x = Day, y = MeanHeight, color = Diversity)) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_point(data = daily_means, aes(x = Day, y = MeanHeight), size = 1) + 
  geom_line(data = daily_means, aes(x = Day, y = MeanHeight)) + 
  ylab('Mean height (cm)') +
  xlab('Julian day of year') +
  scale_colour_manual(name = 'Diversity', values = SpeciesCol) +  
  theme_bw() +
  theme(axis.title = element_text(size = 9)) +
  theme(axis.text = element_text(size = 9, color = "black")) +
  theme(panel.grid = element_blank()) +
  theme(legend.position = c(0.8,0.25),legend.box = 'horizontal',legend.key.size = unit(0.5,'cm'),legend.box.spacing = unit(0.4,'cm'),legend.key.height = unit(0.3,'cm'),legend.box.margin = margin(0,0,0,0,'cm'),legend.box.just = 'bottom',legend.box.background = element_blank(),legend.spacing = unit(0,'cm'), legend.title = element_text(size=8),legend.text = element_text(size=8))
  

# Ramet vs Day
fig.Ramet <- ggplot(MonthlyData, aes(x = Day, y = Ramet, color = Diversity)) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_point(data = daily_means, aes(x = Day, y = Ramet), size = 1) + 
  geom_line(data = daily_means, aes(x = Day, y = Ramet)) + 
  ylab('Ramet number') +
  xlab('Julian day of year') +
  scale_colour_manual(name = 'Diversity', values = SpeciesCol) +  
  theme_bw() +
  theme(axis.title = element_text(size = 9)) +
  theme(axis.text = element_text(size = 9, color = "black")) +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none")

# Inflorescence vs Day
fig.Inflorescence <- ggplot(MonthlyData, aes(x = Day, y = Inflorescence, color = Diversity)) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_point(data = daily_means, aes(x = Day, y = Inflorescence), size = 1) + 
  geom_line(data = daily_means, aes(x = Day, y = Inflorescence)) + 
  ylab('Inflorescence number') +
  xlab('Julian day of year') +
  scale_colour_manual(name = 'Diversity', values = SpeciesCol) +  
  theme_bw() +
  theme(axis.title = element_text(size = 9)) +
  theme(axis.text = element_text(size = 9, color = "black")) +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none")


fig.day <- ggpubr::ggarrange(
  fig.MeanHeight, fig.Ramet, fig.Inflorescence,
  ncol = 3, nrow = 1, labels = c('(a)', '(b)', '(c)'),
  font.label = list(size = 10, color = "black", face = "bold"),
  common.legend = T
)

fig.day
#ggsave('Result/fig.day.pdf',width = 15,height =6,units = c('cm'))#Figure S4
