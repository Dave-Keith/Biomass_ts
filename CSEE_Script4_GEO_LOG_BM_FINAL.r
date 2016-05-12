#####CSEE Script 4 Biomass geo log 
#########################################################################################################
###CSEE####
####see Biomass_script april 27.R
#########LOG GEO#################LOG GEO##################LOG GEO#################
#adjusted GEOMEAN

rm(list=ls(all=T))

library(matlab)
library(lme4)
library(arm)
detach("package:nlme")

load("d:/Dropbox/ASD/Analyses/F_SSB.RData")
#load("d:/Dropbox/ASD/Analyses/F_SSB.RData")
geo.mean <- function(x,n) exp(mean(log(x),na.rm=T))
#help(write.csv)
setwd("d:/Dropbox/My_Papers/Biomass_and_Age/Results")
unique.stocks <- as.character(unique(ASD_FR$Stock.ID))
num.stocks <- length(unique.stocks)
GEO <- zeros(num.stocks,6)

#trouble="COD4VsW" (i=18) and "HAD4X5Y" (i=41)

for(i in 1:num.stocks)  
{
  if(i != 69 && i != 72 && i != 125 && i != 19 && i !=18 && i != 9 && i != 29 && i !=41)
  {
    name <- unique.stocks[i]
    order<-as.character(ASD_FR$Order[ASD_FR$Stock.ID == name ][1])
    st <- ASD_FR[ASD_FR$Stock.ID == name,grep("BM",names(ASD_FR))] 
    check <- unlist(c(st))
    if(length(na.omit(check)) > 0)
    {
      # Need to get rid of rows full of NA's...
      kay <- 0
      #k=18
      for(k in 1:length(st[,1]))
      {
        if(length(na.omit(as.numeric(st[k,]))) == 0) kay <- c(kay,k)
        
      }
      # kay is the 
      if(length(kay) > 1) st <- st[-kay,]  
      
      if(min(st,na.rm=T) == 0) 
      {
        mins <- 0.1*min(st[st>0],na.rm=T)
        j=17
        for(j in 1:length(st[1,]))
        {
          if(length(which(na.omit(st[,j])==0)) >0) st[(st[,j]==0),j] <- mins
        }
      }
      ages <- which(colSums(st,na.rm=T) >0)-1 
      quan <- quantile(ages) 
      LQ <- round(quan[1]):round(quan[2])
      UQ <- round(quan[4]):round(quan[5]) 
      years <- 5
      ts_len <- length(st[,1]) 
      BM_LQ_start <- st[1:years,LQ+1]
      BM_UQ_start <-  st[1:years,UQ+1]
      BM_LQ_now <- st[(ts_len-years+1):ts_len,LQ+1]
      BM_UQ_now <- st[(ts_len-years+1):ts_len,UQ+1]
      BM_LQ_start_tot <- rowSums(BM_LQ_start)
      BM_UQ_start_tot <- rowSums(BM_UQ_start)
      BM_LQ_now_tot <- rowSums(BM_LQ_now)
      BM_UQ_now_tot <- rowSums(BM_UQ_now)
      BM_geo_LQS <- log(geo.mean(BM_LQ_start_tot,years)) #geomean of youngest over first 5 years
      BM_geo_UQS <- log(geo.mean(BM_UQ_start_tot,years)) #geomean of oldest over first 5 years
      BM_geo_LQN <- log(geo.mean(BM_LQ_now_tot,years)) #geomean of youngest over last 5 years (now)
      BM_geo_UQN <- log(geo.mean(BM_UQ_now_tot,years)) #geomean of oldest over last 5 years (now)
      BM_dec_LQ <- (BM_geo_LQN - BM_geo_LQS) 
      BM_dec_UQ <- (BM_geo_UQN - BM_geo_UQS) 
      percent_difference_old<-((BM_geo_UQN-BM_geo_UQS)/BM_geo_UQS)*100
      percent_difference_young<-((BM_geo_LQN-BM_geo_LQS)/BM_geo_LQS)*100 #geometric mean
      GEO[i,]<-c(name, order, BM_dec_UQ, percent_difference_old, BM_dec_LQ, percent_difference_young)
      
    }
  }
}



bad.stocks <- c(which(GEO[,1] == 0), attr(na.omit(GEO[,3]),"na.action"))
Geo <- GEO[-bad.stocks,]



write.csv(Geo, file="CSEE_Script4_GEO_LOG_BM_FINAL_NEW.csv")
# save(GEO,file="GEO.rdata")

Geo.2 <- rbind(Geo[,1:4],Geo[,c(1,2,5,6)])
Geo.2 <- as.data.frame(Geo.2,stringsAsFactors=F)
Geo.2$Age <- c(rep("Old",length(Geo.2[,1])/2),rep("Young",length(Geo.2[,1])/2))
colnames(Geo.2) <- c("Stock","Order","DBM","PDBM","Age")
str(Geo.2)
Geo.2$DBM <- as.numeric(Geo.2$DBM)
Geo.2$PDBM <- as.numeric(Geo.2$PDBM)

qqnorm(Geo.2$DBM)
qqline(Geo.2$DBM)

#mod <- lm(PDBM ~ Age*Order,Geo.2)
#mod <- lm(DBM ~ Age*Order,Geo.2)
#best.mod <- step(mod)
#summary(mod)
#summary(best.mod)
#plot(mod)
Geo.2$Age <- as.factor(Geo.2$Age)
Geo.2$Order <- as.factor(Geo.2$Order)

DBM<-Geo.2$DBM
Age<-Geo.2$Age
Order<-Geo.2$Order
# lets peak at the data...
boxplot(DBM~Age)
boxplot(DBM~Order)
lm.mod <- lm(DBM~Age*Order,Geo.2)
summary(lm.mod)
# This picks out the best model based on AIC model selection...
best.mod <- step(lm.mod)
# So we see the trend is negative, as expected, but doesn't differ across orders, or between Age classes 
summary(best.mod)
# Note a couple points with huge leverage, are these errors?
par(mfrow=c(2,2))
plot(best.mod)
par(mfrow=c(1,1))
#Certainly seems to be nothing, those p-values are huge, I think that's why you get the 0's in the lmer, nothing
# seems to be going on!!
anova(lm.mod)
mod.bm.log <- lmer(DBM ~ Age + (1+Age|Order)) 
summary(mod.bm.log)
anova(mod.bm.log)

ranef(mod.bm.log)
fixef(mod.bm.log)
tot<-(ranef(mod.bm.log))$Order + fixef(mod.bm.log)
tot

  
oldflat<-Geo.2[Geo.2$Order == "Pleuronectiformes" & Geo.2$Age == "Old",]
mean(oldflat$DBM)
youngflat<-Geo.2[Geo.2$Order == "Pleuronectiformes" & Geo.2$Age == "Young",]
mean(youngflat$DBM)
oldgad<-Geo.2[Geo.2$Order == "Gadiformes" & Geo.2$Age == "Old",]
mean(oldgad$DBM)
younggad<-Geo.2[Geo.2$Order == "Gadiformes" & Geo.2$Age == "Young",]
mean(younggad$DBM)
oldherr<-Geo.2[Geo.2$Order == "Clupeiformes" & Geo.2$Age == "Old",]
mean(oldherr$DBM)
youngherr<-Geo.2[Geo.2$Order == "Clupeiformes" & Geo.2$Age == "Young",]
mean(youngherr$DBM)
oldperc<-Geo.2[Geo.2$Order == "Perciformes" & Geo.2$Age == "Old",]
mean(oldperc$DBM)
youngperc<-Geo.2[Geo.2$Order == "Perciformes" & Geo.2$Age == "Young",]
mean(youngperc$DBM)
oldscor<-Geo.2[Geo.2$Order == "Scorpaeniformes" & Geo.2$Age == "Old",]
mean(oldscor$DBM)
youngscor<-Geo.2[Geo.2$Order == "Scorpaeniformes" & Geo.2$Age == "Young",]
mean(youngscor$DBM)

######## percentage diff in geomean using logged biomass of individual ########
######%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#######


mod.bm.log.p <- lmer(PDBM ~ Age + (1+Age|Order),Geo.2) 
summary(mod.bm.log.p)
anova(mod.bm.log.p)
ranef(mod.bm.log.p)
fixef(mod.bm.log.p)
tot.p<-(ranef(mod.bm.log.p))$Order + fixef(mod.bm.log)
tot.p

qqnorm(Geo.2$PDBM)
qqline(Geo.2$PDBM)

oldflat<-Geo.2[Geo.2$Order == "Pleuronectiformes" & Geo.2$Age == "Old",]
mean(oldflat$PDBM)
youngflat<-Geo.2[Geo.2$Order == "Pleuronectiformes" & Geo.2$Age == "Young",]
mean(youngflat$PDBM)
oldgad<-Geo.2[Geo.2$Order == "Gadiformes" & Geo.2$Age == "Old",]
mean(oldgad$PDBM)
younggad<-Geo.2[Geo.2$Order == "Gadiformes" & Geo.2$Age == "Young",]
mean(younggad$PDBM)
oldherr<-Geo.2[Geo.2$Order == "Clupeiformes" & Geo.2$Age == "Old",]
mean(oldherr$PDBM)
youngherr<-Geo.2[Geo.2$Order == "Clupeiformes" & Geo.2$Age == "Young",]
mean(youngherr$PDBM)
oldperc<-Geo.2[Geo.2$Order == "Perciformes" & Geo.2$Age == "Old",]
mean(oldperc$PDBM)
youngperc<-Geo.2[Geo.2$Order == "Perciformes" & Geo.2$Age == "Young",]
mean(youngperc$PDBM)
oldscor<-Geo.2[Geo.2$Order == "Scorpaeniformes" & Geo.2$Age == "Old",]
mean(oldscor$PDBM)
youngscor<-Geo.2[Geo.2$Order == "Scorpaeniformes" & Geo.2$Age == "Young",]
mean(youngscor$PDBM)




#REMOVE SCORPAENIFORMES
#remove osmeriformes and perciformes and scorpaeniformes
colnames(Geo.2) <- c("Stock","Order","DBM","PDBM","Age")
str(Geo.2)
Geo.2$Order
rm <- c(which(Geo.2$Order == "Osmeriformes"),which(Geo.2$Order == "Scorpaeniformes")) #eerst samenvoegen
NWE<- Geo.2[-rm,] #dan eruit halen
colnames(NWE) <- c("Stock","Order","DBM","PDBM","Age")
str(NWE)
NWE$Order
NWE$Age <- as.factor(NWE$Age)
NWE$Order <- as.factor(NWE$Order)

names(NWE)
qqnorm(NWE$DBM)
qqline(NWE$DBM)
str(NWE)
names(NWE)


mod1 <- lmer(DBM ~ Age + (1+ Age|Order),NWE) #My current favourite one
summary(mod1)
anova(mod1)
ranef(mod1)
fixef(mod1)
tot<-(ranef(mod1)$Order + fixef(mod1))
tot



