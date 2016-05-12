#########################################################################################################
###
###             Here I'm processing the data to get it how I want it for the analyses and figures...
rm(list=ls(all=T))
#ASD <- read.table("/home/keithdm/Documents/PhD/R_programs/Age_structured_database/ASD_Final.csv",
#                  header=T,sep =",")#,na.strings="-9999")
load("C:/Users/Emilie/Dropbox/ASD/Analyses/Processed_DB.RData")
#load("D:/Dropbox/ASD/Analyses/Processed_DB.RData")
# This ran the below loop to tidy up the database into something usable!
#setwd("D:/Dropbox/ASD/Analyses")
#load("F_ssb.RData")
# I'll get fancy to show you everything R, first lets make a function that calculates Geometric mean
geo.mean <- function(x,n) prod(x)^(1/n) #so this is taking the product of the sums of the years and then ^1/n
# Try this with i = 5, just rund the stuff between the {}, don't run the whole for loop yet..
i=5
# This will do some Biomass calcs for each stock...
#for(i in 1:num.stocks)
#{
   # First pick the stock we want, look at unique.stocks, and stocks if you need to understand these better
    name <- unique.stocks[i]
    # Now pick the Biomass or abundance for this stock
    st <- BM.calc[stocks == name,] 
#BM.calc is the Num.0:Num.31 for all the stocks, so including those for which there is no Num available
#BM.calc[stocks == name,] is the biomass in numbers for stock i=5
#stock i=5 is AMPL5YZ, which has st from NUM.1:Num.11
    # Now if there is no data we skip everything using this if command, if you run this after runing the above 2 lines with i = 1 you'll see nothing happens
#    if(length(na.omit(st)) > 0)
#{
    # Now pick the youngest age classes, going for lowest quantile, while ignoring NA's
    # This tells you which age classes have data, the minus one corrects for year 0...
    ages <- which(colSums(st,na.rm=T) >0)-1 #so ages for i=5 results in ages 1-11   
#colSums(st,na.rm=T) shows the column sums for each Num.x category
    # This calculates the quantiles for the ages for which we have data
    quan <- quantile(ages) 
    # quan[2] gives 3.5 here, which rounds of to 4 in the next step
# THis picks out the lower and upper quantiles
    LQ <- round(quan[1]):round(quan[2])
    UQ <- round(quan[4]):round(quan[5]) #so this gives the ages 8 9 10 11 
    # Now we can go pick out the data that belongs to these quantiles
    # First lets define the number of years we want to average across at 
    # the start and end of the time series.  Figure out how long
    # the time series is, by looking at the length of one column of the BM data
    years <- 5
    ts_len <- length(st[,1]) #for i=5 this is 28, so how many years we have data for 
    # Now pick out the first five years for the upper and lower quantiles.
    # 1:years picks rows 1:5 which are the first 5 years of data
    # LQ+1 picks out the colums we want, remember +1 because of 0 age class..
    BM_LQ_start <- st[1:years,LQ+1]
    BM_UQ_start <-  st[1:years,UQ+1]
    # To pick the last 5 years, go from the length of the time series -5(+1) to 
    # the end of the time series which ts_len captures...
    BM_LQ_now <- st[(ts_len-years+1):ts_len,LQ+1]
    BM_UQ_now <- st[(ts_len-years+1):ts_len,UQ+1]
    # Now calculate the mean for each year...
    BM_LQ_start_tot <- rowSums(BM_LQ_start)
    BM_UQ_start_tot <- rowSums(BM_UQ_start)
    BM_LQ_now_tot <- rowSums(BM_LQ_now)
    BM_UQ_now_tot <- rowSums(BM_UQ_now)
    # And this is how you calculate a geometric mean using the function I made up above, x is the time series, and n is the number of years
    # over which we have data...
    BM_geo_LQS <- geo.mean(BM_LQ_start_tot,years) #geomean of youngest over first 5 years
    BM_geo_UQS <- geo.mean(BM_UQ_start_tot,years) #geomean of oldest over first 5 years
    BM_geo_LQN <- geo.mean(BM_LQ_now_tot,years) #geomean of youngest over last 5 years (now)
    BM_geo_UQN <- geo.mean(BM_UQ_now_tot,years) #geomean of oldest over last 5 years (now)
    # Now quickly calculate the total decline, negative if decline, positive is growth...
    BM_dec_LQ <- (BM_geo_LQN - BM_geo_LQS) 
    BM_dec_UQ <- (BM_geo_UQN - BM_geo_UQS) 
#}
percent_difference_old<-((BM_geo_UQN-BM_geo_UQS)/BM_geo_UQS)*100
percent_difference_young<-((BM_geo_LQN-BM_geo_LQS)/BM_geo_LQS)*100 #geometric mean
#results for in graph  ***************************^^^^^^^^^^^^^^^^^^^^^^^^^********
name
#old
BM_dec_UQ
percent_difference_old
#young
BM_dec_LQ
percent_difference_young

#!# BM_young_now_then<-c(BM_geo_LQN,BM_geo_LQS)
#!# BM_old_now_then<-c(BM_geo_UQN,BM_geo_UQS)
#!# max.young<-as.numeric(max(BM_young_now_then))
#!# max.old<-max(BM_old_now_then)
#!# BM_dec_LQ_plot<-round((BM_dec_LQ), digits=3)
#!# BM_dec_UQ_plot<-round((BM_dec_UQ), digits=3)       
# Simple Horizontal Bar Plot with Added Labels 
#!# par(mfrow=c(2,1))
#!# barplot(BM_young_now_then, main="Change in geometric mean abundance of young fish", horiz=TRUE,names.arg=c("now", "then"), col=c("lightblue","grey"), xlim=c(0,max.young))
#!# mtext(side=3,paste(BM_dec_LQ_plot,name), cex=1.2)
#!# barplot(BM_old_now_then, main="Change in geometric mean abundance of old fish", horiz=TRUE,names.arg=c("now", "then"), col=c("darkblue","grey"), xlim=c(0,max.old))
#!# mtext(side=3,paste(BM_dec_UQ_plot,name), cex=1.2)

#############################     overall trend of log of total abundance    #######################
#this is taking the sum of each year of all ages and plotting that over time
BM<-(st[1:length(ages)+1])
BM_tot<-rowSums(BM)
x<-1:length(BM_tot)
#!# plot(BM_tot~x, type="l",ylab="Total abundance x 1000", xlab="Year", main="Trend in total abundance over time")
#!# mtext(side=3,paste(name), cex=1.2)
#!# abline(lm(BM_tot~x))
lm.tot.mod<-lm(lm(BM_tot~x))
summary(lm.tot.mod)
#now in log
BM<-(st[1:length(ages)+1])
BM_tot<-rowSums(BM)
BM_tot_log<-log(BM_tot)
x<-1:length(BM_tot_log)
#!# plot(BM_tot_log~x, type="l",ylab="log total abundance young x 1000", xlab="Year", main="Trend in total abundance (log) over time")
#!# mtext(side=3,paste(name), cex=1.2)
#!# abline(lm(BM_tot_log~x))
lm.tot.log.mod<-lm(lm(BM_tot_log~x))
summary(lm.tot.log.mod)

#change in trend over time for oldest quartile fish
BM_UQ <-  st[UQ+1]
AV_BM_UQ <-rowMeans(BM_UQ, na.rm = TRUE)
y=AV_BM_UQ
x = seq(1,length(AV_BM_UQ))
AV_BM_UQ = cbind(x,y)
AV_BM_UQ=as.matrix(AV_BM_UQ)
y=as.matrix(AV_BM_UQ)
x=as.matrix(x)
lm.mod.old<-lm(log(AV_BM_UQ[,2])~AV_BM_UQ[,1])
summary(lm.mod.old)
#!# plot(log(AV_BM_UQ[,2])~AV_BM_UQ[,1], type="l",ylab="Mean abundance (log) old x 1000", xlab="Year", main="Trend of old fish abundance (log) over time")
#!# mtext(side=3,paste(name), cex=1.2)
#!# abline(lm(log(AV_BM_UQ[,2])~AV_BM_UQ[,1]), col="darkblue")

#change in trend over time for youngest quartile fish
BM_LQ <- st[1:years]
AV_BM_LQ <-rowMeans(BM_LQ, na.rm = TRUE)
y=AV_BM_LQ
x = seq(1,length(AV_BM_LQ))
AV_BM_LQ = cbind(x,y)
AV_BM_LQ=as.matrix(AV_BM_LQ)
y=as.matrix(AV_BM_LQ)
x=as.matrix(x)
lm.mod.young<-lm(log(AV_BM_LQ[,2])~AV_BM_LQ[,1])
summary(lm.mod.young)

#!# plot(log(AV_BM_LQ[,2])~AV_BM_LQ[,1], type="l",ylab="mean abundance (log) of young x 1000", xlab="Year", main="Trend of young fish abundance (log) over time")
#!# mtext(side=3,paste(name), cex=1.2)
#!# mod <- lm(log(AV_BM_LQ[,2])~AV_BM_LQ[,1])
#!# lines(coef(mod)[1] +coef(mod)[2] *AV_BM_LQ[,1], col="lightblue")
#plot(AV_BM_LQ[,2]~AV_BM_LQ[,1], type="l",ylab="Mean of nr of young x 1000", xlab="Year", main="Young fish abundance over time")
#abline(lm(AV_BM_LQ[,2]~AV_BM_LQ[,1]),col="lightblue")


############
    # Can see both decline, about 7400 tonnes less in the Youngest age classes, and 8700 fewer tonnes in the oldest now
#    } # End the if loop 
 #      } # End the for loop
#save(ASD_FR,file="F_ssb.RData")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#Current biomass of young fish (1st quantile) and old fish (4th quantile)  compared to the maximum of 
#young and old fish over the assessed years, using a running average
#average of last three years compared to running average of three years
#or I do highest mean of the youngest ages over the last 5 years compared to alltime max
#years<-5
#BM_LQ_now <- st[(ts_len-years+1):ts_len,LQ+1] #this is the biomass of youngest q. for last 5 years
#AV_BM_LQ_now<-rowMeans(BM_LQ_now, na.rm = T)
#Max.now<-max(AV_BM_LQ_now)
#compare that to biomass across all the years for the youngest quantile: 
#BM_LQ <- st[1:years,ages+1] wrong, this gives me all ages for first 5 years and I need all years for first ages
#BM_LQ <- st[1:years]
#AV_BM_LQ <-rowMeans(BM_LQ, na.rm = TRUE)
#Max.tot<-max(AV_BM_LQ)
#percentofmax<-Max.now/Max.tot #current max is x% of what it used to be
####################################################################################################################
# last 5 year running average of youngest quantile versus maximum 5 year running average of youngest quantile
#to show which ones I pick
BM_tot_matrix<-as.matrix(BM_tot)
years=5
BM_LQ <- st[1:years]
AV_BM_LQ <-rowMeans(BM_LQ, na.rm = TRUE)
i=1
RA <- NA
for(i in 1:(ts_len-4))
{
  RA[i]= (AV_BM_LQ[i]+AV_BM_LQ[i+1]+AV_BM_LQ[i+2]+AV_BM_LQ[i+3]+AV_BM_LQ[i+4])/5
}
peryoung <- RA[ts_len-4]/max(RA,na.rm=T)
peryoung #percent of maximum 5 year running average
#plot in barplot
RA_now<-RA[ts_len-4]
RA_max<-max(RA,na.rm=T)
RA_now_max_young<-c(RA_now, RA_max)
# Simple Horizontal Bar Plot with Added Labels 

barplot(RA_now_max_young, main="Running average of young fish: current vs. max", horiz=TRUE,names.arg=c("now", "max"), col=c("lightblue","grey"), xlim=c(0,max(RA_now_max_young)))
mtext(side=3,paste(name), cex=1.2)
#######################  old   #####################################
# last 5 year running average of oldest quartile versus maximum 5 year running average of oldest quantile
BM_UQ <-  st[UQ+1]
AV_BM_UQ <-rowMeans(BM_UQ, na.rm = TRUE)
i=1
RA <- NA
for(i in 1:(ts_len-4))
{
  RA[i]= (AV_BM_UQ[i]+AV_BM_UQ[i+1]+AV_BM_UQ[i+2]+AV_BM_UQ[i+3]+AV_BM_UQ[i+4])/5
}
perold <- RA[ts_len-4]/max(RA,na.rm=T) #RA[ts_len-4] is the running average of the last five years of the oldest quantile
perold #percent of maximum 5 year running average

RA_now<-RA[ts_len-4]
RA_max<-max(RA,na.rm=T)
RA_now_max_old<-c(RA_now, RA_max)
# Simple Horizontal Bar Plot with Added Labels 
barplot(RA_now_max_old, main="Running average of old fish: current vs. max", horiz=TRUE,names.arg=c("now", "max"), col=c("darkblue","grey"), xlim=c(0,max(RA_now_max_old)))
mtext(side=3,paste(name), cex=1.2)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#



#plot(residuals(lm.mod))
lm.mod<-lm(log(AV_BM_LQ[,2])~AV_BM_LQ[,1])
#summary(lm.mod)
#plot(lm.mod)
res<-as.matrix(residuals(lm.mod))
acf.res.mod<-acf(res)

cors <-as.numeric(unlist(acf.res.mod[1])[1])
#pacf(res)

#lm.mod2<-lm(AV_BM_LQ[,2]~AV_BM_LQ[,1])
#summary(lm.mod2)
#because of dependence
library(nlme)
gls.mod<-gls(log(AV_BM_LQ[,2])~AV_BM_LQ[,1], correlation=corAR1(cors))
#plot(residuals(gls.mod))
acf(gls.mod$residuals)
#summary(gls.mod)
#names(gls.mod)
#plot(gls.mod)
#print(gls.mod)



#plot(lm.mod)
#plot(residuals(lm.mod))
#res<-as.matrix(residuals(lm.mod))
#acf.res.mod<-acf(res)
#cors <-as.numeric(unlist(acf.res.mod[1])[1])
#pacf(res)

#plot(AV_BM_UQ[,2]~AV_BM_UQ[,1], type="l",ylab="Mean abundance old x 1000", xlab="Year", main="Mean of old fish abundance over time")
#abline(lm(AV_BM_UQ[,2]~AV_BM_UQ[,1]), col="darkblue")
#lm.mod2<-lm(AV_BM_UQ[,2]~AV_BM_UQ[,1])
#summary(lm.mod2)

#gls.mod<-gls(log(AV_BM_UQ[,2])~AV_BM_UQ[,1], correlation=corAR1(cors))
#resid(gls.mod)
#acf(gls.mod$residuals)
#summary(gls.mod)


####################
#As a rule, negative skewness indicates that the mean of the data values is 
#less than the median, and the data distribution is left-skewed; positive skewness 
#would indicates that the mean of the data values is larger 
#than the median, and the data distribution is right-skewed. 
library(moments)                # load the moments package 
#skewness()               # apply the skewness function


ts_len <- length(st[,1]) #for i=5 this is 28, so how many years we have data for 
BM<-(st[1:length(ages)+1])
row.names(BM)<-1:ts_len


#colors()

pdf("agedist6.pdf") 
for(i in 1:ts_len)
{ 
  BM_year<-as.matrix(BM[i,])
  ages<-colnames(BM_year)
  barplot(BM_year, xlab="Age Category", col="seagreen3")  
  mtext(side=3,paste("Skew = ",round(skewness(as.numeric(BM_year)),digits=2),"          Year = ",row.names(BM_year)), cex=1.2)
mtext(side=4,paste(name), cex=1.2)
  }
dev.off()


#BM_year<-as.matrix(BM[6,])
#skewness(BM_year[1,])

skew <- NA
for(i in 1:ts_len)
{ 
  BM_year<-as.matrix(BM[i,])
  skew[i]<-(skewness(BM_year[1,]))
}


year<-1:ts_len
plot(year[skew<0],skew[skew<0], type="p", pch=19, col="darkblue",xlab="year",ylab="skew",
     ylim=c(min(skew),max(skew)), xlim=c(min(year), max(year)))
mtext(side=3,paste(name), cex=1.2)
points(year[skew>=0],skew[skew>=0], pch=19, col="lightblue")
abline(0,0)


#regression from scratch
#n=length(AV_BM_LQ[,1])
#Z=cbind(rep(1,n),AV_BM_LQ[,1])
#Z=as.matrix(Z)
#n=length(Z[,1])
#k=length(Z[1,])
#p=k+1

#compute Z'Z
#ZtZ=t(Z)%*%Z
#ZtZi=solve(ZtZ)
#Zty=t(Z)%*%x
#betahat=ZtZi%*%Zty
#xhat=Z%*%betahat

#plot the results and the data
#plot(AV_BM_LQ[,1],AV_BM_LQ[,2], type="b")
#lines(AV_BM_LQ[,1], xhat, col="blue")

#residual analysis
#e=x-xhat

#residual plots
#par(mfrow=c(2,2))
#plot(xhat,e,type="b")
#plot(Z[,2],e,type="b")
#qqnorm(e)
#e.acf=acf(e)


#significance testing
#MSRes<-(sum((e)^2))/(n-p) #compute mean square error of the residuals
#test significance of the slope
#beta1hat=betahat[2]
#sebeta1hat=sqrt(MSRes*ZtZi[2,2])
#t=beta1hat/sebeta1hat
#pvalols=1-pt(t,n-p)
#pvalols

#assessing and modelling the residual autocorrelation
#V=correlation matrix of the residuals
#determine lag 1 autocorrelation
#rho=e.acf$acf[2]
#or rho=cor(e[2:n], e[1:n-1])

#make corresponding correlation matrix for the errors
#rhovec=matrix(NA,1,n)
#for (i in 1:n) 
#  {
#  rhovec[i]=rho^(i-1)
#  }
#rhovec=as.numeric(rhovec) 
#V=toeplitz(rhovec)
#a.eig<-eigen(V)
#a.eig
#make the square root matrix required for transformations
#K<- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*%solve(a.eig$vectors)
#Ki=solve(K) #inverse of K
#Ki
#Generalized LS regression
#transform the data
#x=Ki%*%x
#Z=Ki%*%Z
#x
#Z
#complete the GLS regression
#hint using these transformed variables, you can directly apply the ols code above
#plot(x)
#compute Z'Z
#ZtZ=t(Z)%*%Z
#ZtZi=solve(ZtZ)
#Zty=t(Z)%*%x
#betahat=ZtZi%*%Zty
#xhat=Z%*%betahat

#plot the results and the data
#plot the results and the data
#plot(AV_BM_LQ[,1],AV_BM_LQ[,2], type="b")
#lines(AV_BM_LQ[,1], xhat, col="blue")

#residual analysis
#e=x-xhat

#residual plots
#par(mfrow=c(2,2))
#plot(xhat,e,type="b")
#plot(Z[,2],e,type="b")
#qqnorm(e)
#e.acf=acf(e)
#plot(GT[,1],GT[,2], type="b")
#lines(GT[,1], xhat, col="blue")

#residual analysis
#e=x-xhat
#residual plots
#par(mfrow=c(2,2))
#plot(xhat,e,type="b")
#plot(Z[,2],e,type="b")
#qqnorm(e)
#e.acf=acf(e)
#residuals versus Z looks funny
#so log axis of Z
#Z
#log.Z=log(Z)
#residual plots with logged Z
#par(mfrow=c(2,2))
#plot(xhat,e,type="b")
#plot(log.Z[,2],e,type="b")
#qqnorm(e)
#e.acf=acf(e)
#still looks funny, so cut off the axis beyond which there are outliers
#residual plots with logged Z and getting rid of outliers (1st and 39th value of Z)
#par(mfrow=c(2,2))
#plot(xhat,e,type="b")
#plot(Z[,2],xlim=c(1268,1334),e,type="b")
#qqnorm(e)
#e.acf=acf(e)
#Z

#significance testing
#MSRes<-(sum((e)^2))/(n-p) #compute mean square error of the residuals
#test significance of the slope
#beta1hat=betahat[2]
#sebeta1hat=sqrt(MSRes*ZtZi[2,2])
#t=beta1hat/sebeta1hat
#pvalgls=1-pt(t,n-p)
#t
#pvalgls
#show matrix that is used for transformation of the data (to take out autocorrelation)
#V
#fit<-lm()
#plot(t,fit$resid, type="o", ylab="detrended global air temperature after transformation of data using GLS")


#meaningful division, based on reproductive value
#how to calculate reproductive value?

