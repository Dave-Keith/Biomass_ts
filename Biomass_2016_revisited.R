### This script is used to develop the analysis looking at the trends in some measure of abundance for each age class of a
###  a population, Questions I want to answer
### 1:  What is the pattern in change of biomass/aboundance over time, does it differ between old and young
### 2:  What is the current status of old vs. young component of these populations current (relate to max/min historical levels)
### 3:  There might be a cool question about catch in here (though that could be another paper), what percentage of our catch 
###     is coming from young/old components of the population and how has that changed over time??


#############  Section 1 Data Processing#############  Section 1 Data Processing#############  Section 1 Data Processing

rm(list=ls(all=T))
library(matlab)
library(lme4)
library(arm)
library(gamm4)
library(lattice)
library(akima)
library(mapplots)
# For linux
#direct = "/media/sf_data/r/"
yr = as.numeric(format(Sys.time(), "%Y")) -1

# re-load rprofile if you have one...
source(".Rprofile")
#Also load in a couple of R files that will help with plotting, these are from the TESA GAMM course by Zuur that occured in Jan 2016
# See DK to get these if you need them.  I placed the scripts in the same folder as my .Rprofile so they would be easily accessable.
# Put these aren't mine to make public!!
direct = "d:/r/"
source(paste(direct,"DIYBiplot.R",sep=""))
source(paste(direct,"HighstatLibV9.R",sep=""))
#setwd("d:/Dropbox/My_Papers/Biomass_and_age/")
setwd("d:/Github/Current_papers/Biomass_ts")

ASD_MGMT<-read.csv("ASD_MGMT.csv", header=T,stringsAsFactors = F)
# Currently the mid-point of most time series is 1989...
mean.year <- round(mean(ASD_MGMT$Year))
ASD_MGMT$Year_cen <- ASD_MGMT$Year-mean.year
geo.mean <- function(x,n) prod(x)^(1/n) 

unique.stocks <- as.character(unique(ASD_MGMT$Stock.ID))
num.stocks <- length(unique.stocks)


# What are we going to pick as our variable to look at?
# Options are any of the age specific variables... "SSB","BM","Num","Catch","FM",
var <- "BM"

db <- NULL
ages <- NULL
age.quan <- NULL
young.ages <- NULL
old.ages <- NULL
i=18
# so we want to run a loop to grab the biomass estiamtes for each stock
for(i in 1:num.stocks)  
{

  # Several of these stocks have a male and female entry, so we'll exclude the males
  # Get the stock information....
  name <- unique.stocks[i]
  if(is.element(name,c("ALPLAICBSAIm","ARFLOUNDBSAIm","GHALBSAIm","NRSOLEEBSAIm","YSOLEBSAIm"))==F)
  {
    # Grab the names of the stock descriptors to keep and then subset the data by stock and to pull the appropriate
    # variable defined by var above.
    dat.names <- c("Stock.ID","Management","Area","Order","Family","Genus","Species","LME","Year","Year_cen")
    dat <- subset(ASD_MGMT, Stock.ID == name ,select = c(dat.names,names(ASD_MGMT[,grep(var,names(ASD_MGMT))])))
    # Determine which (if any) columns and rows have data    
    cols <- which(colSums(dat[,grep(var,names(dat))],na.rm=T) > 0)
    rows <- which(rowSums(subset(dat,select = -c(which(names(dat) %in% dat.names))),na.rm=T) == 0)
    # Remove any columns with no data
    if(length(cols) > 0)
    {
      db[[name]] <- dat[,c(1:length(dat.names),cols + length(dat.names))]
    } # end if(length(cols) > 0)
    # Remove any years with no data
    if(length(rows) > 0)
    {
      db[[name]] <- db[[name]][-rows,]
    } # end if(length(rows) > 0)
    
    # If we have data...
    if(is.null(db[[name]]) == F)
    {
      # OK, this is a bizarro way to get the ages but it does the trick! Following line grabs the variable names...
      ages[[name]]  <- as.numeric(substr(names(db[[name]][(grep("[[:digit:]]",
                                                                names(db[[name]])))]),start=(nchar(var)+2),stop=(nchar(var)+3)))
      var.names <- names(db[[name]][(grep("[[:digit:]]",names(db[[name]])))])
      
      # Now we break the ages into old and young quantiles, this splits the data into the oldest 25% and youngest 25% of ages
      # Note that we don't think of + groups and that we round so that often it's not 25% (e.g. sometimes 20%, others may be 30%)
      age.quan[[name]] <- quantile(ages[[name]]) 
      #First the young
      young.ages[[name]] <- round(age.quan[[name]][1]):round(age.quan[[name]][2])
      var.names.young <- var.names[1:length(young.ages[[name]])]
      # Now for the old
      old.ages[[name]] <- round(age.quan[[name]][4]):round(age.quan[[name]][5]) 
      var.names.old <- var.names[(length(ages[[name]]) - length(old.ages[[name]] ) +1) : length(ages[[name]])]
      
      # A few differnt ways to define our response variable, stan might be the most statisically 
      # pleasing one, I think ratio is the easiest transformation to wrap my head around and to
      # use when comparing models but statistical properties are wonderful
      # The offset might work nicely in terms of the model too and avoid transformation concerns...
      db[[name]]$total <- rowSums(subset(db[[name]],select = -c(which(names(db[[name]]) %in% dat.names))),na.rm=T)
      # Replace 0's with half minimum non-zero value in time series just so log works 
      if(any(db[[name]]$total == 0))db[[name]]$total[db[[name]]$total == 0] <- 0.5*min(db[[name]]$total[db[[name]]$total > 0])
      db[[name]]$total_ratio   <- db[[name]]$total/max(db[[name]]$total,na.rm=T)
      db[[name]]$total_stan   <- scale(log(db[[name]]$total))
      db[[name]]$total_offset   <- max(db[[name]]$total,na.rm=T)
      
      db[[name]]$old   <- rowSums(subset(db[[name]],select = c(which(names(db[[name]]) %in% var.names.old))),na.rm=T)
      # Replace 0's with half minimum value in time series just so log works 
      if(any(db[[name]]$old == 0))db[[name]]$old[db[[name]]$old == 0] <- 0.5*min(db[[name]]$old[db[[name]]$old > 0])
      db[[name]]$old_ratio   <- db[[name]]$old/max(db[[name]]$old,na.rm=T)
      db[[name]]$old_stan   <- scale(log(db[[name]]$old))
      db[[name]]$old_offset   <- max(db[[name]]$old,na.rm=T)
      
      db[[name]]$young <- rowSums(subset(db[[name]],select = c(which(names(db[[name]]) %in% var.names.young))),na.rm=T)
      # Replace 0's with half minimum value in time series just so log works 
      if(any(db[[name]]$young == 0))db[[name]]$young[db[[name]]$young == 0] <- 0.5*min(db[[name]]$young[db[[name]]$young > 0])
      db[[name]]$young_ratio   <- db[[name]]$young/max(db[[name]]$young,na.rm=T)
      db[[name]]$young_offset   <- max(db[[name]]$young,na.rm=T)
      db[[name]]$young_stan   <- scale(log(db[[name]]$young))
      # Now we can remove the individual age data as that's not especially necessary to keep
      db[[name]] <- subset(db[[name]],select = -c(which(names(db[[name]]) %in% var.names)))
    } # end if(length(cols) > 0 && length(rows) > 0)
  } # end if(is.element(name,c("ALPLAICBSAIm","ARFLOUNDBSAIm","GHALBSAIm","NRSOLEEBSAIm","YSOLEBSAIm"))==F)
} # end for(i in 1:num.stocks)  
# Make sure the columns all the same length
range(lapply(db,ncol))

db <- do.call("rbind",db)

#############  End Section 1 Data Processing#############  End Section 1 Data Processing#############  End Section 1 Data Processing


#############  Section 2 Data Exploration#############  Section 2 Data Exploration#############  Section 2 Data Exploration
#############  Section 2 Data Exploration#############  Section 2 Data Exploration#############  Section 2 Data Exploration
#############  Section 2 Data Exploration#############  Section 2 Data Exploration#############  Section 2 Data Exploration
### OK, now we have the data, time to start walking through the models properly
### First lets take a look at our data....

# Let's look and see what we have over time
# Well even in 1946 we have 3 populations which is nice, by 1950 we have 5, 1955 there are 10 and at the
# peak we have 1000, drops back to 10 by 2010 but hopefully we can update that!!  We have 100 in the late 90's and  over
# 20 stocks from 1961 onwards, hit 50 in 1976 and are at 80 in 1982
aggregate(total_stan~Year,db,FUN=length)

range(aggregate(total_stan~Stock.ID,db,FUN=length)$V1) # Everything has at least 10 years
# On average the time series are 33 years long, nice!
median(aggregate(total_stan~Stock.ID,db,FUN=length)$V1)

# How does this break down by other covariates?  ICES and NOAA dominate the analysis, maybe 10% of data are other locations
aggregate(total_stan~Management,db,FUN=length)

# First we need to set up the plotting area, I hate lattice so I'm doing this my way...
# Set the number of rows
nr <- ceiling(sqrt(length(unique(db$Stock.ID))))
# Set the number of columns, the funky little command is used to add one to the nc if nr is a perfect square
# As we are printing nchains + 1
ifelse(sqrt(length(unique(db$Stock.ID))) %% 1==0,  nc <- ceiling(sqrt(length(unique(db$Stock.ID))))+1, 
       nc <- ceiling(sqrt(length(unique(db$Stock.ID)))))
windows(15,15)
ly <- layout(matrix(c(1:(nr*nc)), nr, nc, byrow = T))
layout.show(ly)
par(mar=c(0,0,2,0))
for(i in 1:length(unique(db$Stock.ID))) 
{
  with(subset(db,Stock.ID==unique(db$Stock.ID)[i]),plot(old_stan~Year,pch=16,cex=0.01,xaxt="n",yaxt="n",bty="U"))
  title(unique(db$Stock.ID)[i],cex.main=0.7)
} # end for(i in 1:length(unique(mw$year))) 

# Now look at the young folks...
windows(15,15)
ly <- layout(matrix(c(1:(nr*nc)), nr, nc, byrow = T))
layout.show(ly)
par(mar=c(0,0,2,0))
for(i in 1:length(unique(db$Stock.ID))) 
{
  with(subset(db,Stock.ID==unique(db$Stock.ID)[i]),plot(young_stan~Year,pch=16,cex=0.01,xaxt="n",yaxt="n",bty="U"))
  title(unique(db$Stock.ID)[i],cex.main=0.7)
} # end for(i in 1:length(unique(mw$year))) 

# Now look at the total
windows(15,15)
ly <- layout(matrix(c(1:(nr*nc)), nr, nc, byrow = T))
layout.show(ly)
par(mar=c(0,0,2,0))
for(i in 1:length(unique(db$Stock.ID))) 
{
  with(subset(db,Stock.ID==unique(db$Stock.ID)[i]),plot(total_stan~Year,pch=16,cex=0.01,xaxt="n",yaxt="n",bty="U"))
  title(unique(db$Stock.ID)[i],cex.main=0.7)
} # end for(i in 1:length(unique(mw$year))) 


# Now the histograms of the shell heights, because the bin size is 10's for these the y-axis * 10 gives us
# the proportion within each bin.  The last time we saw a big spike in any one bin was really 2007 (about 30% were 100-110)
# Otherwise is it oddly smooth proportion of the samples in each size bin.
windows(20,15)
ly <- layout(matrix(c(1:(nr*nc)), nr, nc, byrow = T))
par(mar=c(2,1,1,2))
xmax <- ceiling(max(mw$sh)/10)*10) # This rounds the xmax to the tens place (e.g. 172 becomes 180).
for(i in 1:length(unique(mw$year))) 
{
  if(is.element(i,seq(nr,nr*nc,by=nr))==F)   
  {
    with(subset(mw,year==sort(unique(mw$year))[i]),hist(sh,xlim=c(50,xmax),cex=0.5,ylim=c(0,0.05),main="",yaxt="n",xaxt="n",freq=F))
    with(subset(mw,year==sort(unique(mw$year))[i]),text(0.9*xmax,0.8*0.05,unique(year),cex=1.5))
    axis(2,labels=F)
    axis(1,cex.axis=0.8)
  } # end if(i %in% seq(nr,nr*nc,by=nr))   
  
  if(i %in% seq(nr,nr*nc,by=nr))   
  {
    with(subset(mw,year==sort(unique(mw$year))[i]),hist(sh,xlim=c(50,xmax),cex=0.5,ylim=c(0,0.05),main="",yaxt="n",xaxt="n",freq=F))
    with(subset(mw,year==sort(unique(mw$year))[i]),text(0.9*xmax,0.8*0.05,unique(year),cex=1.5))
    axis(4)
    axis(1,cex.axis=0.8)
  } # end if(i %in% seq(nr,nr*nc,by=nr))   
  if(i == length(unique(mw$year))) axis(4)
} # end for(i in 1:length(unique(mw$year))) 

# These all seems reasonable, certainly some difference in what we are seeing depending on the population...


# Now a dotplot of the data looking for outliers and such...
windows(11,8.5)
Vars <- c("Year_cen","total","total_ratio","total_stan","total_offset","old"  ,        "old_ratio"   ,
          "old_stan"  ,   "old_offset"  , "young"     ,   "young_ratio" , "young_offset", "young_stan"  )
Mydotplot(db[,Vars])
# You can see the problem with using the raw data, the actual values are so askew for some stocks
# The ratio is o.k. but having that artifical ceiling at 1 is a pain, the standardized data does look nice
# but of course interpretation can be a bit funny.

# Now let's look for co-linearity in our data, of course there should be all sorts of it.
# Can see that the young Biomass is very much correlated with total biomass, but so is old
# so that's fine.  Not much here since we don't have any real numeric covaraiates...
windows(11,8.5)
Vars <- c("Year_cen","total","total_ratio","total_stan","total_offset","old"  ,        "old_ratio"   ,
          "old_stan"  ,   "old_offset"  , "young"     ,   "young_ratio" , "young_offset", "young_stan"  )
pairs(db[,Vars],lower.panel=panel.cor)


aggregate(young_stan~Year,db,FUN=mean) # How does the mean change each year
windows(11,8.5)
plot(aggregate(young_stan~Year,db,FUN=mean),main="Young Standaridzed mean") # Can see a pretty steady decline in this!!
windows(11,8.5)
plot(aggregate(old_stan~Year,db,FUN=mean),main="Old Standaridzed mean") # See that rebound in 2000, curious
windows(11,8.5)
plot(aggregate(total_stan~Year,db,FUN=mean),main="Total Standaridzed mean") # more cyclic then the others, again some 2000 rebound evidence.
windows(11,8.5)
plot(aggregate(young_ratio~Year,db,FUN=mean),main="Young Ratio mean") # Can see a pretty steady decline in this!!
windows(11,8.5)
plot(aggregate(old_ratio~Year,db,FUN=mean),main="Old Ratio mean") # collapse as we move into the 70's (given new stocks coming
# online this collapse is a bit surprising to me!!)
windows(11,8.5)
plot(aggregate(total_ratio~Year,db,FUN=mean),main="Total Ratio mean") #


# Looking at the ratio data
aggregate(old_ratio~Year,db,FUN=sd) # How does the variance change each year
aggregate(young_ratio~Year,db,FUN=sd) # How does the variance change each year
aggregate(total_ratio~Year,db,FUN=sd) # How does the variance change each year
windows(11,8.5)
par(mfrow=c(3,1))
plot(aggregate(old_ratio~Year,db,FUN=sd),main="Old ratio Variance") # Can see early in the time series we might have a problem 
# since the stocks are all starting at maximum in late 1940's, but it evens out very quickly.  There is a steady decline
# which is likely a function of the high ratios disappearing.
plot(aggregate(young_ratio~Year,db,FUN=sd),main="Young ratio Variance") # For the young the data actually looks solid variance wize.
plot(aggregate(total_ratio~Year,db,FUN=sd),main="Total ratio Variance") # For the total see it is low early as well, good by the 60's

# Now the standarized data
aggregate(old_stan~Year,db,FUN=sd) # How does the variance change each year
aggregate(young_stan~Year,db,FUN=sd) # How does the variance change each year
aggregate(total_stan~Year,db,FUN=sd) # How does the variance change each year
windows(11,8.5)
par(mfrow=c(3,1))
plot(aggregate(old_stan~Year,db,FUN=sd),main="Old stan Variance") # Slightly surprising pattern, but again after 1960 very steady
plot(aggregate(young_stan~Year,db,FUN=sd),main="Young stan Variance") # Very low early then steady since 60's
plot(aggregate(total_stan~Year,db,FUN=sd),main="Total stan Variance") # Steady since the 1960's

# Now let's look at the data against some of our covariates...

# First up management alone
aggregate(old_ratio~Management,db,FUN=mean)
aggregate(young_ratio~Management,db,FUN=mean)
aggregate(total_ratio~Management,db,FUN=mean)
# These are interesting somewhat surprisingly, can see DFO/NOAA stocks and NAFO stocks are clearly doing worse
# DFO may also be a bit low.  ICES and NOAA look to be doing better big picture.  Everyone's old are depressed more
# than 60%, not sure any are much different but DFO does look a bit lower
windows(11,8.5)
par(mfrow=c(3,1))
boxplot(old_ratio~Management,db,main="old")
boxplot(young_ratio~Management,db,main="young")
boxplot(total_ratio~Management,db,main="total")

# Stan isn't interesting as they all tend to 0 (because they are standardized by stock, need to remove year to get interesting)

# Can see the database is mostly Herring and Cod, but only Scorpaeniformes is really poorly represented.
table(subset(db,select=c("Management","Order")))

# Here's boxplot of the ratios by year and management, can see DFO again looks bad!!
# DK note:  Notice also that NOAA in the early years has data that looks very strange!!
windows(15,8.5)
par(mfrow=c(3,1))
boxplot(old_ratio~interaction(Year,Management),db)
boxplot(young_ratio~interaction(Year,Management),db)
boxplot(total_ratio~interaction(Year,Management),db)

# Here's the standarized data...
windows(15,8.5)
par(mfrow=c(3,1))
boxplot(old_stan~interaction(Year,Management),db)
boxplot(young_stan~interaction(Year,Management),db)
boxplot(total_stan~interaction(Year,Management),db)

# What about taxonomy, really only enough probably to look at herring vs. cod...


# First up order alone, looks like herrings and cods been hit the hardest, 
# The old ages again showing this, if we look at total most things averaging around 50% or greater.
aggregate(old_ratio~Order,db,FUN=mean)
aggregate(young_ratio~Order,db,FUN=mean)
aggregate(total_ratio~Order,db,FUN=mean)
# These are interesting somewhat surprisingly, can see DFO/NOAA stocks and NAFO stocks are clearly doing worse
# DFO may also be a bit low.  ICES and NOAA look to be doing better big picture.  Everyone's old are depressed more
# than 60%, not sure any are much different but DFO does look a bit lower
windows(11,8.5)
par(mfrow=c(3,1))
boxplot(old_ratio~Order,db,main="old")
boxplot(young_ratio~Order,db,main="young")
boxplot(total_ratio~Order,db,main="total")

# Stan isn't interesting as they all tend to 0 (because they are standardized by stock, need to remove year to get interesting)

# Here's boxplot of the ratios by year and Order, can see Gads look worst, Clups were bad but seemingly on comback trail.
# There is an upswing in old individuals for Clups, Gads, and Pleuro's that you don't see for young or total biomass.
# DK note:  Notice also that Perciformes in the early years has data that looks very strange, this must be the NOAA pattern too...
windows(15,8.5)
par(mfrow=c(3,1))
boxplot(old_ratio~interaction(Year,Order),db)
boxplot(young_ratio~interaction(Year,Order),db)
boxplot(total_ratio~interaction(Year,Order),db)

# Here's the standarized data...
windows(15,8.5)
par(mfrow=c(3,1))
boxplot(old_stan~interaction(Year,Order),db)
boxplot(young_stan~interaction(Year,Order),db)
boxplot(total_stan~interaction(Year,Order),db)

####################  END SECTION 1 ####################  END SECTION 1 .####################  END SECTION 1 ##########################
####################  END SECTION 1 ####################  END SECTION 1 .####################  END SECTION 1 ##########################


# Well let's start at the start, how has abundance of old farts changed over time...
mod.1.old.stan <- lm(old_stan~Year,data=db) # note their used to be a mod.1 when we had outliers in the data, they've been taken care of so mod.1 is gone..
mod.1.old.rat <- lm(old_ratio~Year,data=db) # note their used to be a mod.1 when we had outliers in the data, they've been taken care of so mod.1 is gone..

# What's it look like, you can see the models are "significant' but they hardly explain anything...
summary(mod.1.old.stan) 
summary(mod.1.old.rat) 

# How are the residuals, honestly they ain't bad!!
windows(11,8.5)
par(mfrow=c(2,2))
plot(mod.1.old.stan) 

# Now these are terrible!!
windows(11,8.5)
par(mfrow=c(2,2))
plot(mod.1.old.rat) 


# The fit of the model to the data...
windows(11,8.5)
plot(old_ratio~Year,data=db,pch=16)

E1.old.rat <- resid(mod.1.old.rat, type = "pearson")
Disp.rat <- sum(E1.old.rat^2) / mod.1.old.rat$df.res
Disp.rat

E1.old.stan <- resid(mod.1.old.stan, type = "pearson")
Disp.stan <- sum(E1.old.stan^2) / mod.1.old.stan$df.res
Disp.stan

F1.old.stan <- fitted(mod.1.old.stan)
F1.old.rat <- fitted(mod.1.old.rat)

# You can kinda see the residuals are too high below fitted values of -0.2 and above around 0.2 (not enough low values up there)
windows(11,8.5)
par(mfrow=c(1,1))
plot(x = F1.old.stan, 
     y = E1.old.stan,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h=0,col="blue",lty=2)

# What is the residual trend?
gam.E1.old.stan <- gam(E1.old.stan~te(F1.old.stan))
summary(gam.E1.old.stan) # Only 2% of residual variance explained, but significant trend..
windows(11,8.5)
plot(gam.E1.old.stan) # Clearly a trend in residuals with E being too high at both low and high fitted values...

# What is the residual trend?
gam.E1.old.rat <- gam(E1.old.rat~te(F1.old.rat))
summary(gam.E1.old.rat) # Only 3% of residual variance explained
windows(11,8.5)
plot(gam.E1.old.rat) #Clearly residuals are too high at both low and high fitted values

#Next how do the residuals look versus Year , I think it's related to the above, not enough low residuals early (when populations)
# were doing good, and not enough low late, which is suggesting a rebound for the old individuals here I think...
windows(11,8.5)
plot(E1.old.stan ~db$Year)
abline(h=0,col="blue",lty=2)

# Let's do a gam on these residuals...
center_year <- db$Year_cen
gam.E1.old.stan.year <- gam(E1.old.stan~te(center_year))
summary(gam.E1.old.stan.year) # Only 2% of residual variance explained, but significant trend..
windows(11,8.5)
plot(gam.E1.old.stan.year) # Same trend more or less as we saw above...

# What do the residuals look like versus Managment...

windows(11,8.5)
par(mfrow=c(2,1))
boxplot(E1.old.stan~db$Management, main="Old Stan") # nothing exciting here and not sure there could be...
boxplot(E1.old.rat~db$Management, main="Old Ratio") # NAFO and DFO's look low when we look this way...

# How about against taxonomy...

windows(11,8.5)
par(mfrow=c(2,1))
boxplot(E1.old.stan~db$Order, main="Old Stan") # nothing exciting here and not sure there could be
boxplot(E1.old.rat~db$Order, main="Old Ratio") # Clupeiforms and Gadiforms both looking low...


# What's our buddy autocorrelation saying about the residuals...
# Looks like we have at least an AR1, possibly an AR2...
windows(11,8.5)
par(mfrow=c(2,2))
acf(E1.old.stan)
pacf(E1.old.stan)
acf(E1.old.rat)
pacf(E1.old.rat)


##########  So really what we have here are counts which change over time, what we want is a glm, but given the 
##########  major differences in counts between populations we need to include an offset, so a poisson glm with offset seems reasonable alternative..
##########  Using the quasipoisson because the data are not integers and is likely super over-dispersed anyway...
mod.2.old.poisson <- glm(old~Year,data=db,offset=log(old_offset),family=quasipoisson) 

summary(mod.2.old.poisson)


# For the model checking I need to look at the model we developed with Njal!!!
E2.old.poisson <- resid(mod.2.old.poisson, type = "deviance")
Disp.stan <- sum(E2.old.poisson^2) / mod.2.old.poisson$df.res
Disp.stan

# I think that's actually the fitted values of interest... I think...
F2.old.poisson <- fitted(mod.2.old.poisson)/db$old_offset

# How are the residuals, gross!!
windows(11,8.5)
par(mfrow=c(2,2))
plot(mod.2.old.poisson) 


# You can kinda see the residuals are too high below fitted values of -0.2 and above around 0.2 (not enough low values up there)
windows(11,8.5)
par(mfrow=c(1,1))
plot(x = F2.old.poisson, 
     y = E2.old.poisson,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h=0,col="blue",lty=2)

































# Some initial poking around I'll revisit later.
plot(db[[name]]$total~db[[name]]$Year,log="y",type="o",col="black",
     pch=15,ylim=c(min(db[[name]]$old[db[[name]]$old>0]),max(db[[name]]$total)))
lines(db[[name]]$old~db[[name]]$Year,type="o",col="blue",pch=16)
lines(db[[name]]$young~db[[name]]$Year,type="o",col="orange",pch=17)

tst.off <- glm(young~(Year_cen) ,data=db,family=poisson)
plot(tst.off)
summary(tst.off)
summary(tst.off$lme)
tst <- gamm4(young_stan~s(Year_cen), random = ~ (1 | Stock.ID),data=db)
summary(tst$mer)
summary(tst$gam)
plot(tst$gam)
plot(tst$mer) # Nice residuals using stan!!

tst2 <- gamm4(old_stan~s(Year_cen), random = ~ (1 | Stock.ID),data=db)
summary(tst2$mer)
summary(tst2$gam) 

# So overall this suggests a decline occured for the first part of the time series followed by a recovery since the late 70's
# the recovery seems to have slowed more recently
plot(tst2$gam)
plot(tst2$mer) # Nice residuals using stan again!!

# Might be a bit of auto-correlation there, maybe ;-)
acf(tst2$gam$residuals)
# save(GEO,file="GEO.rdata")
CSEE.ALL <- rbind(CSEE[,1:6],CSEE[,c(1,2,3,4,7,8)]) #divide young and old
CSEE.ALL <- as.data.frame(CSEE.ALL,stringsAsFactors=F)
CSEE.ALL$Age <- c(rep("Old",length(CSEE.ALL[,1])/2),rep("Young",length(CSEE.ALL[,1])/2))
colnames(CSEE.ALL) <- c("name", "Order", "Species", "MGMT", "DBM", "PDBM", "Age")
str(CSEE.ALL)
CSEE.ALL$DBM <- as.numeric(CSEE.ALL$DBM)
CSEE.ALL$PDBM <- as.numeric(CSEE.ALL$PDBM)
qqnorm(CSEE.ALL$DBM)
qqline(CSEE.ALL$DBM)

# summary(as.factor(CSEE.ALL$Order))

oldflat<-CSEE.ALL[CSEE.ALL$Order == "Pleuronectiformes" & CSEE.ALL$Age == "Old",]
sum(oldflat$DBM)
#total sum of biomass is -11.79725
#exp(11.79)=131926.5, so we lost that in tonnes of biomass?
mean(oldflat$DBM)
median(oldflat$DBM)

summary
youngflat<-CSEE.ALL[CSEE.ALL$Order == "Pleuronectiformes" & CSEE.ALL$Age == "Young",]
sum(youngflat$DBM)
summary(youngflat)
mean(youngflat$DBM)
median(youngflat$DBM)
oldgad<-CSEE.ALL[CSEE.ALL$Order == "Gadiformes" & CSEE.ALL$Age == "Old",]
sum(oldgad$DBM)
mean(oldgad$DBM)
median(oldgad$DBM)
younggad<-CSEE.ALL[CSEE.ALL$Order == "Gadiformes" & CSEE.ALL$Age == "Young",]
mean(younggad$DBM)
median(younggad$DBM)

#1054950834 tons?!?! DAVe is this correct?
younggad$name
#exp(20.77676)=1054,950,834 tonnes
oldherr<-CSEE.ALL[CSEE.ALL$Order == "Clupeiformes" & CSEE.ALL$Age == "Old",]
mean(oldherr$DBM)
median(oldherr$DBM)
sum(oldherr$DBM)
#exp(3.510065)=33.454
youngherr<-CSEE.ALL[CSEE.ALL$Order == "Clupeiformes" & CSEE.ALL$Age == "Young",]
mean(youngherr$DBM)
median(youngherr$DBM)
oldperc<-CSEE.ALL[CSEE.ALL$Order == "Perciformes" & CSEE.ALL$Age == "Old",]
mean(oldperc$DBM)
median(oldperc$DBM)
youngperc<-CSEE.ALL[CSEE.ALL$Order == "Perciformes" & CSEE.ALL$Age == "Young",]
mean(youngperc$DBM)
median(youngperc$DBM)
oldscor<-CSEE.ALL[CSEE.ALL$Order == "Scorpaeniformes" & CSEE.ALL$Age == "Old",]
mean(oldscor$DBM)
youngscor<-CSEE.ALL[CSEE.ALL$Order == "Scorpaeniformes" & CSEE.ALL$Age == "Young",]
mean(youngscor$DBM)
youngosm<-CSEE.ALL[CSEE.ALL$Order == "Osmeriformes" & CSEE.ALL$Age == "Young",]
mean(youngosm$DBM)

oldflat<-CSEE.ALL[CSEE.ALL$Order == "Pleuronectiformes" & CSEE.ALL$Age == "Old",]
mean(oldflat$PDBM)
median(oldflat$PDBM)
youngflat<-CSEE.ALL[CSEE.ALL$Order == "Pleuronectiformes" & CSEE.ALL$Age == "Young",]
mean(youngflat$PDBM)
median(youngflat$PDBM)
oldgad<-CSEE.ALL[CSEE.ALL$Order == "Gadiformes" & CSEE.ALL$Age == "Old",]
mean(oldgad$PDBM)
median(oldgad$PDBM)
younggad<-CSEE.ALL[CSEE.ALL$Order == "Gadiformes" & CSEE.ALL$Age == "Young",]
mean(younggad$PDBM)
median(younggad$PDBM)
oldherr<-CSEE.ALL[CSEE.ALL$Order == "Clupeiformes" & CSEE.ALL$Age == "Old",]
mean(oldherr$PDBM)
median(oldherr$PDBM)
youngherr<-CSEE.ALL[CSEE.ALL$Order == "Clupeiformes" & CSEE.ALL$Age == "Young",]
mean(youngherr$PDBM)
median(youngherr$PDBM)
oldperc<-CSEE.ALL[CSEE.ALL$Order == "Perciformes" & CSEE.ALL$Age == "Old",]
mean(oldperc$PDBM)
median(oldperc$PDBM)
youngperc<-CSEE.ALL[CSEE.ALL$Order == "Perciformes" & CSEE.ALL$Age == "Young",]
mean(youngperc$PDBM)
median(youngperc$PDBM)


# so i will remove the scorpaeniformes, bc it may change the model
#REMOVE SCORPAENIFORMES
CSEE.ALL$Order
rm <- c(which(CSEE.ALL$Order == "Scorpaeniformes")) #eerst samenvoegen
CSEE.ALL<- CSEE.ALL[-rm,] #dan eruit halen
CSEE.ALL$Order
CSEE.ALL$MGMT<- as.factor(CSEE.ALL$MGMT)
CSEE.ALL$Age <- as.factor(CSEE.ALL$Age)
CSEE.ALL$Order <- as.factor(CSEE.ALL$Order)


# remove NAFO, ICCAT #
rm <- c(which(CSEE.ALL$MGMT == "NAFO"))
CSEE.ALL<- CSEE.ALL[-rm,]
CSEE.ALL$MGMT
summary(CSEE.ALL)
write.csv(CSEE.ALL, file="CSEE_ALL_BM.csv")

# models #
# lets peak at the data.."name", "Order", "Species", "MGMT", "DBM", "PDBM", "Age")
rm(list=ls(all=T))
CSEE.ALL<-read.csv("CSEE_ALL_BM.csv", header=T)
name<-CSEE.ALL$name
Order<-CSEE.ALL$Order
Species<-CSEE.ALL$Species
MGMT<-CSEE.ALL$MGMT
DBM<-CSEE.ALL$DBM
PDBM<-CSEE.ALL$PDBM
Age<-CSEE.ALL$Age

summary(CSEE.ALL$MGMT)
# GRAPH 1
summary(lm(DBM~Age,CSEE.ALL))
summary.aov(lm(DBM~Age,CSEE.ALL))
boxplot(DBM~Age, main="Difference in Log of Biomass", col=(c("blue3", "lightblue")))

#Graph 2
summary.aov(lm(DBM~MGMT,CSEE.ALL))
summary(lm(DBM~MGMT,CSEE.ALL))
boxplot(DBM~MGMT, main="Difference in Log of Biomass", col=(c("coral","seagreen", "aliceblue","cadetblue")))

#Graph 3
summary(lm(DBM~Order, CSEE.ALL))
summary.aov(lm(DBM~Order, CSEE.ALL))
boxplot(DBM~Order, main="Difference in Log of Biomass", col=(c("salmon", "lightpink", "hotpink", "palevioletred4")))

#Graph 4
summary(lm(DBM~Order*Age, CSEE.ALL))
summary.aov(lm(DBM~Order*Age, CSEE.ALL))
boxplot(DBM~Age*Order, main="Log Biomass in tonnes",  col=(c("blue3", "lightblue")))

#########################################
#                                     ###
#       fancy plot DBM~....         ###
#                                     ###
#########################################

summary(lm(DBM~Age-1)) #gives estimate for young and old
exp(-0.337)-1    #28% we lost of old biomass
exp(-0.377)-1  #we lost 30.9% of our young biomass

#I want Graph 1 in a fancy plot incl 95% CI
dbmage<-(lm(DBM~Age-1))
coef(dbmage)
ci_dbmage<-confint(dbmage,level=.95)

lci_dbmage<-ci_dbmage[,1]
uci_dbmage<-ci_dbmage[,2]
mn_dbmage<- coef(summary.lm(dbmage))[,1]

elci_dbmage<-(exp(lci_dbmage)-1)*100
euci_dbmage<-(exp(uci_dbmage)-1)*100
emn_dbmage<-(exp(mn_dbmage)-1)*100

par(mar=c(5,5,5,13),cex=1)
#old
plot(emn_dbmage[c(1)],1.5, ylim=c(1,length(emn_dbmage)+1),
     bty="U",pch=19,xlim=c(min(elci_dbmage),max(euci_dbmage)),ylab="",yaxt="n",xlab="", cex=1)
#young
points(emn_dbmage[c(2)],2.5,pch=19, font=2)
#points(mn_BM_I_peak[c(5:8)],6:9,pch=19)
# Now I add the CI's
segments(elci_dbmage[c(1)],1.5,euci_dbmage[c(1)],1.5)
segments(elci_dbmage[c(2)],2.5,euci_dbmage[c(2)],2.5)
# Throw a line in to divide them
abline(h=c(2),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
axis(2,at=c(1.5,2.5),labels=c("Old","Young"),las=1,tck=0, cex.axis=1.2)
mtext(side=1,"Loss of Biomass in %",line=3, cex=1.2)
#axis(4,at=c(1,2,3,4,5,6,7,8),labels=c("Clupeiformes","Gadiformes","Perciformes","Pleuronectiformes",
  #                                    "Clupeiformes","Gadiformes","Perciformes","Pleuronectiformes"),las=1)

#222222222222222222222222222222222222222222222222222222222222222
###GRAPH 2 fancy
summary.aov(lm(DBM~MGMT,CSEE.ALL))
summary(lm(DBM~MGMT,CSEE.ALL))
boxplot(DBM~MGMT, main="Difference in Log of Biomass", col=(c("coral","seagreen", "aliceblue","cadetblue")))
summary(MGMT)
dm<-(lm(DBM~MGMT-1))
coef(dm)
ci_dm<-confint(dm,level=.95)

lci_dm<-ci_dm[,1]
uci_dm<-ci_dm[,2]
mn_dm<- coef(summary.lm(dm))[,1]

elci_dm<-(exp(lci_dm)-1)*100
euci_dm<-(exp(uci_dm)-1)*100
emn_dm<-(exp(mn_dm)-1)*100

par(mar=c(5,5,5,13),cex=1)
#mgmt
plot(emn_dm[c(1:4)],1.5:4.5, ylim=c(1,length(emn_dm)+1),
     bty="U",pch=19,xlim=c(min(elci_dm),max(euci_dm)),ylab="",yaxt="n",xlab="", cex=1)
#points(mn_BM_I_peak[c(5:8)],6:9,pch=19)
# Now I add the CI's
segments(elci_dm[c(1:4)],1.5:4.5,euci_dm[c(1:4)],1.5:4.5)
# Throw a line in to divide them
abline(h=c(2,3,4),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
mtext(side=1,"Change in Biomass (%)",line=3, cex=1.2)
axis(4,at=c(1.5,2.5,3.5,4.5),labels=c("DFO (8)","DFO+NOAA (6)","ICES (104)","NOAA (64)"),las=1, cex.axis=1.2)


#GRAPH 3
#DBM~ORDER
dbmage<-(lm(DBM~Order-1))
coef(dbmage)
ci_dbmage<-confint(dbmage,level=.95)

lci_dbmage<-ci_dbmage[,1]
uci_dbmage<-ci_dbmage[,2]
mn_dbmage<- coef(summary.lm(dbmage))[,1]

elci_dbmage<-(exp(lci_dbmage)-1)*100
euci_dbmage<-(exp(uci_dbmage)-1)*100
emn_dbmage<-(exp(mn_dbmage)-1)*100

par(mar=c(5,5,5,13),cex=1)
#order
plot(emn_dbmage[c(1:4)],1.5:4.5, ylim=c(1,length(emn_dbmage)+1),
     bty="U",pch=19,xlim=c(min(elci_dbmage),max(euci_dbmage)),ylab="",yaxt="n",xlab="", cex=1)
#points(mn_BM_I_peak[c(5:8)],6:9,pch=19)
# Now I add the CI's
segments(elci_dbmage[c(1:4)],1.5:4.5,euci_dbmage[c(1:4)],1.5:4.5)
# Throw a line in to divide them
abline(h=c(2,3,4),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
mtext(side=1,"Change in Biomass (%)",line=3, cex=1.2)
axis(4,at=c(1.5,2.5,3.5,4.5),labels=c("Clupeiformes","Gadiformes","Perciformes","Pleuronectiformes"),las=1, cex.axis=1.2)


#I want graph 4 in a fancy plot, incl 95% CI
dbmageor<-(lm(DBM~Age*Order-Age-Order-1))
coef(dbmageor)
ci_dbmageor<-confint(dbmageor,level=.95)

lci_dbmageor<-ci_dbmageor[,1]
uci_dbmageor<-ci_dbmageor[,2]
mn_dbmageor<- coef(summary.lm(dbmageor))[,1]

elci_dbmageor<-(exp(lci_dbmageor)-1)*100
euci_dbmageor<-(exp(uci_dbmageor)-1)*100
emn_dbmageor<-(exp(mn_dbmageor)-1)*100


# And now I do the plot..
par(mar=c(5,5,5,13),cex=1)
# So first I pull out of the mn_BM_I_peak / mn_dbmageor object the first 4 taxonomic groups
#(you'll have to customize yourself) using the c(1:4) from the DI category
#old
OLD<-c(1,3,5,7)
plot(emn_dbmageor[c(1,3,5,7)],1:4, ylim=c(1,length(emn_dbmageor)+1),
    bty="U",pch=19,xlim=c(min(elci_dbmageor),max(euci_dbmageor)),ylab="",yaxt="n",xlab="", cex=1)

#plot(mn_BM_I_peak[c(1:4)],1:4,ylim=c(1,length(mn_BM_I_peak)+2),
 #    bty="U",pch=19,xlim=c(min(LCI_BMI_peak),max(UCI_BMI_peak)),ylab="",yaxt="n",xlab="") 

# Then I get the ones that are NDD / age y
#young
YOUNG<-c(2,4,6,8)
points(emn_dbmageor[c(2,4,6,8)],6:9,pch=19, cex=1)
#points(mn_BM_I_peak[c(5:8)],6:9,pch=19)

# And finally the PDD / age old
#points(mn_BM_I_peak[c(9:10)],11:12,pch=19)
# Now I add the CI's

segments(elci_dbmageor[c(1,3,5,7)],1:4,euci_dbmageor[c(1,3,5,7)],1:4)
segments(elci_dbmageor[c(2,4,6,8)],6:9,euci_dbmageor[c(2,4,6,8)],6:9)

# Throw a line in to divide them
abline(h=c(5),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
axis(2,at=c(2.5,7.5,11.5),labels=c("Old","Young",""),las=1,tck=0, cex.axis=1.2)
mtext(side=1,"Loss of Biomass in %",line=3, cex=1.2)
axis(4,at=c(1,2,3,4,6,7,8,9),labels=c("Clupeiformes","Gadiformes","Perciformes","Pleuronectiformes",
                                            "Clupeiformes","Gadiformes","Perciformes","Pleuronectiformes"),las=1, cex.axis=1.2)
#segments(LCI_BMI_peak[c(1:4)],1:4,UCI_BMI_peak[c(1:4)],1:4)
#segments(LCI_BMI_peak[c(5:8)],6:9,UCI_BMI_peak[c(5:8)],6:9)
#segments(LCI_BMI_peak[c(9:10)],11:12,UCI_BMI_peak[c(9:10)],11:12)
as.list(emn_dbmageor)

#summary.aov(lm(DBM~Species*Age,CSEE.ALL))
#summary(lm(DBM~Species*Age,CSEE.ALL))
#summary.aov(lm(DBM~Age*MGMT,CSEE.ALL))

#boxplot(DBM~Order*MGMT, main="Difference in Log of Biomass",) 
        #col=(c("coral","salmon","seagreen","lightpink", "aliceblue","hotpink","cadetblue","palevioletred4")))
#percentage loss of log biomass
boxplot(PDBM~Species)
summary(lm(PDBM~Age*Order))
boxplot(PDBM~Age*Order)
boxplot(PDBM~Order)
boxplot(PDBM~Age)
#boxplot(PDBM~MGMT,  main="Log of Biomass", col=(c("coral","seagreen", "aliceblue","cadetblue")))

summary(lm(DBM~Age,CSEE.ALL))
lm.mod.1 <- lm(DBM~Age*Order,CSEE.ALL)
summary.aov(lm.mod.1)
plot(DBM~Age*Order)
summary(lm(DBM~Age*Order,CSEE.ALL))
best.mod <- step(lm.mod.1)
anova(lm.mod.1)

lm.mod.2 <-lm(DBM~Age*MGMT,CSEE.ALL)
summary(lm.mod.2)
plot(DBM~MGMT)
?pairwise.t.test
best.mod <- step(lm.mod.1)
anova(lm.mod.2)

lm.mod.3 <-lm(DBM~Age*MGMT*Order,CSEE.ALL)
summary(lm.mod.3)
best.mod <- step(lm.mod.3)
plot(lm.mod.3)
anova(lm.mod.3)

lm.mod.4 <-lm(DBM~MGMT, CSEE.ALL)
summary(lm.mod.4)

best.mod <- step(lm.mod.4)
plot(lm.mod.4)
anova(lm.mod.4)

lm.mod.5 <-lm(DBM~Age,CSEE.ALL)
summary(lm.mod.5)
best.mod <- step(lm.mod.1)
anova(lm.mod.2)

lm.mod<-lm(DBM~Age*MGMT*Order)
summary(lm.mod)
#???pairwise.t.test(DBM,MGMT)

#THE ONLY SIGNIFICANT MODEL!!!!!
lm.mod.6<-lm(DBM~Age*Species)
summary(lm.mod.6)
summary.aov(lm.mod.6)

lm.mod.7<-lm(DBM~Age+Species+MGMT)
summary(lm.mod.7)
summary.aov(lm.mod.7)

best<-step(lm.mod.7, direction="both")
best<-step(lm.mod.7, direction="forward")
best<-step(lm.mod.7, direction="backward")

lm.mod.8<-lm(DBM~Age*Species)
summary(lm.mod.8)
best<-step(lm.mod.8, direction="both")
best<-step(lm.mod.8, direction="forward")
best<-step(lm.mod.8, direction="backward")
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

#mod.bm.log <- lmer(DBM ~ Age + (1+Age|Order),CSEE.ALL) 
#summary(mod.bm.log)
#anova(mod.bm.log)
#ranef(mod.bm.log)
#fixef(mod.bm.log)
#tot<-(ranef(mod.bm.log))$Order + fixef(mod.bm.log)
#tot
#qqnorm(CSEE.ALL$PDBM)
#qqline(CSEE.ALL$PDBM)

#lets try everything
#transform data
ALLBMLOG <- data.frame(DBM,Order, Species, MGMT, Age)
m.full<-lm(DBM ~ ., ALLBMLOG) #the "." means "everything else"
summary.aov(m.full)
summary(m.full)
#does other direction matter?no tiny bit only
ALLBMLOGDIR<-data.frame(DBM,Age,MGMT,Order,Species)
mfulldir<-lm(DBM~.,ALLBMLOGDIR)
summary.aov(mfulldir)
summary(mfulldir)

step(mfulldir,direction="both")
step(m.full,direction="both")
step(m.full, direction="forward")
step(m.full, direction="backward")

mod.best<-lm(DBM~Age)
mod.best
summary(mod.best)


#SIMPLE CORRELATIONS TO SEE WHATs going on, needs to be numeric

Order<-as.numeric(CSEE.ALL$Order)
Species<-as.numeric(CSEE.ALL$Species)
MGMT<-as.numeric(CSEE.ALL$MGMT)
DBM<-CSEE.ALL$DBM
PDBM<-CSEE.ALL$PDBM
Age<-as.numeric(CSEE.ALL$Age)

cor(DBM,Age)
cor(DBM,MGMT)
cor(DBM,Species)
cor(DBM,Order)

Order<-as.character(CSEE.ALL$Order)
Species<-as.character(CSEE.ALL$Species)
MGMT<-as.character(CSEE.ALL$MGMT)
DBM<-CSEE.ALL$DBM
PDBM<-CSEE.ALL$PDBM
Age<-as.character(CSEE.ALL$Age)



###################################################
###                                             ###
###     BIOMASS LOGGED ON SPECIES LEVEL         ###
###                   SSBL                      ###
###                                             ###
###################################################

cod<-subset(CSEE.ALL, CSEE.ALL$Species == "morhua")
had<-subset(CSEE.ALL, CSEE.ALL$Species=="aeglefinus")
her<-subset(CSEE.ALL, CSEE.ALL$Species=="harengus")
sole<-subset(CSEE.ALL, CSEE.ALL$Species=="vulgaris")
plai<-subset(CSEE.ALL, CSEE.ALL$Species=="platessa")

cod<-as.data.frame(cod)
had<-as.data.frame(had)
her<-as.data.frame(her)
sole<-as.data.frame(sole)
plai<-as.data.frame(plai)

SPP <- rbind(cod,had,her,sole,plai)
write.csv(SPP, file="SPPBMLOG.csv")
summary(SPP)

#gives 13 morhua, 12 aegle, 11 hareng, 7 vulg, 5 plat, total of 48 stocks
#load new dataframe containing only those species

rm(list=ls(all=T))
SPLIM<-read.csv("SPPBMLOG.csv", header=T)
summary(SPLIM)


DBM<-SPLIM$DBM
PDNUM<-SPLIM$PDBM
Order<-SPLIM$Order
name<-SPLIM$name
Species<-SPLIM$Species
MGMT<-SPLIM$MGMT
Age<-SPLIM$Age
summary(lm(DBM~Species))
summary.aov(lm(DBM~Species))


#GRAPH 5
boxplot(DBM~Species, main="Difference in Log of Biomass for Selected Species",
        col=(c("lightpink", "salmon", "lightpink","palevioletred4", "palevioletred4")))
summary(Species)
summary(lm(DBM~Age*Species))

dbmage<-(lm(DBM~Species-1))
coef(dbmage)
ci_dbmage<-confint(dbmage,level=.95)

lci_dbmage<-ci_dbmage[,1]
uci_dbmage<-ci_dbmage[,2]
mn_dbmage<- coef(summary.lm(dbmage))[,1]

elci_dbmage<-(exp(lci_dbmage)-1)*100
euci_dbmage<-(exp(uci_dbmage)-1)*100
emn_dbmage<-(exp(mn_dbmage)-1)*100


par(mar=c(5,5,5,13),cex=1)
#old
plot(emn_dbmage[c(1:5)],1.5:5.5, ylim=c(1,length(emn_dbmage)+1),
     bty="U",pch=19,xlim=c(min(elci_dbmage),max(euci_dbmage)),ylab="",yaxt="n",xlab="", cex=1)
#points(mn_BM_I_peak[c(5:8)],6:9,pch=19)
# Now I add the CI's
segments(elci_dbmage[c(1:5)],1.5:5.5,euci_dbmage[c(1:5)],1.5:5.5)
# Throw a line in to divide them
abline(h=c(2,3,4,5),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
mtext(side=1,"Loss of Biomass in %",line=3, cex=1.2)
axis(4,at=c(1.5:5.5),labels=c("aeglefinus (12)","harengus (11)","morhua (13)","platessa (5)",
"vulgaris (7)"),las=1, cex.axis=1.1)


#GRAPH 6
boxplot(DBM~Age*Species, main="Difference in Log of Biomass", col=(c("blue3", "lightblue")))
#summary(lmer(DBM~Age+(1+Age|Species)))
summary.aov(lm(DBM~Age*Species))
das<-(lm(DBM~Age*Species-Age-Species-1))
coef(das)
ci_das<-confint(das,level=.95)
lci_das<-ci_das[,1]
uci_das<-ci_das[,2]
mn_das<- coef(summary.lm(das))[,1]
elci_das<-(exp(lci_das)-1)*100
euci_das<-(exp(uci_das)-1)*100
emn_das<-(exp(mn_das)-1)*100
# And now I do the plot..
par(mar=c(5,5,5,13),cex=1)
# So first I pull out of the mn_BM_I_peak / mn_das object the first 4 taxonomic groups
#(you'll have to customize yourself) using the c(1:4) from the DI category
#old
plot(emn_das[c(1,3,5,7,9)],1:5, ylim=c(1,length(emn_das)+1),
     bty="U",pch=19,xlim=c(min(elci_das),max(euci_das)),ylab="",yaxt="n",xlab="", cex=1)
#plot(mn_BM_I_peak[c(1:4)],1:4,ylim=c(1,length(mn_BM_I_peak)+2),
#    bty="U",pch=19,xlim=c(min(LCI_BMI_peak),max(UCI_BMI_peak)),ylab="",yaxt="n",xlab="") 
# Then I get the ones that are NDD / age y
#young
points(emn_das[c(2,4,6,8,10)],6:10,pch=19, cex=1)
#ci
segments(elci_das[c(1,3,5,7,9)],1:5,euci_das[c(1,3,5,7,9)],1:5)
segments(elci_das[c(2,4,6,8,10)],6:10,euci_das[c(2,4,6,8,10)],6:10)
# Throw a line in to divide them
abline(h=c(5.5),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
axis(2,at=c(3,8,11),labels=c("Old","Young",""),las=1,tck=0, cex.axis=1.2)
mtext(side=1,"Loss of Biomass in %",line=3, cex=1.2)
axis(4,at=c(1:10),labels=c("aeglefinus","harengus","morhua","platessa","vulgaris",
               "aeglefinus","harengus","morhua","platessa","vulgaris"),las=1, cex.axis=1)

as.list(emn_das)

##################################################
###                                             ###
###           CSEE RAW BM (not logged) BRAW     ###
###                     all species             ###
###                                             ###
###################################################


rm(list=ls(all=T))
library(matlab)
library(lme4)
library(arm)
setwd("C:/Users/Emilie/Documents/CSEE")
ASD_MGMT<-read.csv("ASD_MGMT.csv", header=T)
names(ASD_MGMT)
ASD_MGMT$Stock.ID
levels(ASD_MGMT$Management)
ASD_MGMT$Species
geo.mean <- function(x,n) prod(x)^(1/n) 

unique.stocks <- as.character(unique(ASD_MGMT$Stock.ID))
num.stocks <- length(unique.stocks)
str(ASD_MGMT)



CSEE.BM.RAW <- zeros(num.stocks,8)

for(i in 1:num.stocks)  
{
  if(i != 69 && i != 72 && i != 125 && i != 19 && i !=18 && i != 9 && i != 29 && i !=41)
  {
    name <- unique.stocks[i]
    Order<-as.character(ASD_MGMT$Order[ASD_MGMT$Stock.ID == name ][1])
    MGMT<-as.character(ASD_MGMT$Management[ASD_MGMT$Stock.ID == name ] [1])
    Species<-as.character(ASD_MGMT$Species[ASD_MGMT$Stock.ID == name] [1])
    st <- ASD_MGMT[ASD_MGMT$Stock.ID == name,grep("BM",names(ASD_MGMT))] 
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
      BM_geo_LQS <- geo.mean(BM_LQ_start_tot,years) #geomean of youngest over first 5 years
      BM_geo_UQS <- geo.mean(BM_UQ_start_tot,years) #geomean of oldest over first 5 years
      BM_geo_LQN <- geo.mean(BM_LQ_now_tot,years) #geomean of youngest over last 5 years (now)
      BM_geo_UQN <- geo.mean(BM_UQ_now_tot,years) #geomean of oldest over last 5 years (now)
      BM_dec_LQ <- (BM_geo_LQN - BM_geo_LQS) 
      BM_dec_UQ <- (BM_geo_UQN - BM_geo_UQS) 
      percent_difference_old<-((BM_geo_UQN-BM_geo_UQS)/BM_geo_UQS)*100
      percent_difference_young<-((BM_geo_LQN-BM_geo_LQS)/BM_geo_LQS)*100 #geometric mean
      CSEE.BM.RAW[i,]<-c(name, Order, Species, MGMT, BM_dec_UQ, percent_difference_old, 
                  BM_dec_LQ, percent_difference_young)
      
    }
  }
}


CSEE.BM.RAW[,]
bad.stocks <- c(which(CSEE.BM.RAW[,6] == 0), attr(na.omit(CSEE.BM.RAW[,5]),"na.action"))
CSEE.BM.RAW <- CSEE.BM.RAW[-bad.stocks,]
CSEE.BM.RAW[,]
# save(GEO,file="GEO.rdata")
CSEE.BM.RAW <- rbind(CSEE.BM.RAW[,1:6],CSEE.BM.RAW[ ,c(1,2,3,4,7,8)]) #divide young and old
CSEE.BM.RAW <- as.data.frame(CSEE.BM.RAW,stringsAsFactors=F)
CSEE.BM.RAW$Age <- c(rep("Old",length(CSEE.BM.RAW[,1])/2),rep("Young",length(CSEE.BM.RAW[,1])/2))
colnames(CSEE.BM.RAW) <- c("name", "Order", "Species", "MGMT", "DBM", "PDBM", "Age")
str(CSEE.BM.RAW)
CSEE.BM.RAW$DBM <- as.numeric(CSEE.BM.RAW$DBM)
CSEE.BM.RAW$PDBM <- as.numeric(CSEE.BM.RAW$PDBM)
qqnorm(CSEE.BM.RAW$DBM)
qqline(CSEE.BM.RAW$DBM)

name<-CSEE.BM.RAW$name
Order<-CSEE.BM.RAW$Order
Species<-CSEE.BM.RAW$Species
MGMT<-CSEE.BM.RAW$MGMT
DBM<-CSEE.BM.RAW$DBM
PDBM<-CSEE.BM.RAW$PDBM
Age<-CSEE.BM.RAW$Age

oldflat<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Pleuronectiformes" & CSEE.BM.RAW$Age == "Old",]
sum(oldflat$DBM)
summary(oldflat$DBM)
youngflat<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Pleuronectiformes" & CSEE.BM.RAW$Age == "Young",]
sum(youngflat$DBM)
oldgad<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Gadiformes" & CSEE.BM.RAW$Age == "Old",]
sum(oldgad$DBM)
younggad<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Gadiformes" & CSEE.BM.RAW$Age == "Young",]
sum(younggad$DBM)
oldherr<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Clupeiformes" & CSEE.BM.RAW$Age == "Old",]
sum(oldherr$DBM)
youngherr<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Clupeiformes" & CSEE.BM.RAW$Age == "Young",]
sum(youngherr$DBM)
oldperc<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Perciformes" & CSEE.BM.RAW$Age == "Old",]
sum(oldperc$DBM)
youngperc<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Perciformes" & CSEE.BM.RAW$Age == "Young",]
sum(youngperc$DBM)
mean(youngperc$DBM)
summary(as.factor(CSEE.BM.RAW$Order))
summary(youngperc)

MEAN.BM.RAW <- zeros(8,4)
colnames(MEAN.BM.RAW) <- c("Order","MEAN", "MEDIAN", "Age")

summary(MEAN.BM.RAW)        

MEAN.BM.RAW[1,1]<-"Pleuronectiformes"
MEAN.BM.RAW[1,2]<-as.numeric(mean(oldflat$DBM))
MEAN.BM.RAW[1,3]<-as.numeric(median(oldflat$DBM))
MEAN.BM.RAW[1,4]<-as.character("Old")

MEAN.BM.RAW[2,1]<-"Pleuronectiformes"
MEAN.BM.RAW[2,2]<-as.numeric(mean(youngflat$DBM))
MEAN.BM.RAW[2,3]<-as.numeric(median(youngflat$DBM))
MEAN.BM.RAW[2,4]<-as.character("Young")

MEAN.BM.RAW[3,1]<-"Gadiformes"
MEAN.BM.RAW[3,2]<-as.numeric(mean(oldgad$DBM, na.rm=T))
MEAN.BM.RAW[3,3]<-as.numeric(median(oldgad$DBM, na.rm=T))
MEAN.BM.RAW[3,4]<-as.character("Old")

MEAN.BM.RAW[4,1]<-"Gadiformes"
MEAN.BM.RAW[4,2]<-as.numeric(mean(younggad$DBM, na.rm=T))
MEAN.BM.RAW[4,3]<-as.numeric(median(younggad$DBM, na.rm=T))
MEAN.BM.RAW[4,4]<-as.character("Young")

MEAN.BM.RAW[5,1]<-"Clupeiformes"
MEAN.BM.RAW[5,2]<-as.numeric(mean(oldherr$DBM))
MEAN.BM.RAW[5,3]<-as.numeric(median(oldherr$DBM))
MEAN.BM.RAW[5,4]<-as.character("Old")

MEAN.BM.RAW[6,1]<-"Clupeiformes"
MEAN.BM.RAW[6,2]<-as.numeric(mean(youngherr$DBM))
MEAN.BM.RAW[6,3]<-as.numeric(median(youngherr$DBM))
MEAN.BM.RAW[6,4]<-as.character("Young")

MEAN.BM.RAW[7,1]<-"Perciformes"
MEAN.BM.RAW[7,2]<-as.numeric(mean(oldperc$DBM))
MEAN.BM.RAW[7,3]<-as.numeric(median(oldperc$DBM))
MEAN.BM.RAW[7,4]<-as.character("Old")

MEAN.BM.RAW[8,1]<-"Perciformes"
MEAN.BM.RAW[8,2]<-as.numeric(mean(youngperc$DBM))
MEAN.BM.RAW[8,3]<-as.numeric(median(youngperc$DBM))
MEAN.BM.RAW[8,4]<-as.character("Young")

#write.csv(MEAN.BM.RAW, file="MEAN_MEDIAN_BM_RAW.csv")
MEAN.BM.RAW<-as.matrix(MEAN.BM.RAW)


oldflat<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Pleuronectiformes" & CSEE.BM.RAW$Age == "Old",]
summary(oldflat$PDBM)
median(oldflat$PDBM)
youngflat<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Pleuronectiformes" & CSEE.BM.RAW$Age == "Young",]
mean(youngflat$PDBM)
median(youngflat$PDBM)
oldgad<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Gadiformes" & CSEE.BM.RAW$Age == "Old",]
mean(oldgad$PDBM, na.rm=T)
median(oldgad$PDBM, na.rm=T)
younggad<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Gadiformes" & CSEE.BM.RAW$Age == "Young",]
mean(younggad$PDBM, na.rm=T)
median(younggad$PDBM, na.rm=T)
oldherr<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Clupeiformes" & CSEE.BM.RAW$Age == "Old",]
mean(oldherr$PDBM,na.rm=T)
median(oldherr$PDBM,na.rm=T)
youngherr<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Clupeiformes" & CSEE.BM.RAW$Age == "Young",]
mean(youngherr$PDBM)
median(youngherr$PDBM,na.rm=T)
oldperc<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Perciformes" & CSEE.BM.RAW$Age == "Old",]
mean(oldperc$PDBM)
median(oldperc$PDBM)
youngperc<-CSEE.BM.RAW[CSEE.BM.RAW$Order == "Perciformes" & CSEE.BM.RAW$Age == "Young",]
mean(youngperc$PDBM)
median(youngperc$PDBM)

#REMOVE SCORPAENIFORMES
CSEE.BM.RAW$Order
rm <- c(which(CSEE.BM.RAW$Order == "Scorpaeniformes")) #eerst samenvoegen
CSEE.BM.RAW<- CSEE.BM.RAW[-rm,] #dan eruit halen
CSEE.BM.RAW$Order

CSEE.BM.RAW$MGMT<- as.factor(CSEE.BM.RAW$MGMT)
CSEE.BM.RAW$Age <- as.factor(CSEE.BM.RAW$Age)
CSEE.BM.RAW$Order <- as.factor(CSEE.BM.RAW$Order)


# remove NAFO, ICCAT #
rm <- c(which(CSEE.BM.RAW$MGMT == "NAFO"))
CSEE.BM.RAW<- CSEE.BM.RAW[-rm,]
CSEE.BM.RAW$MGMT
summary(CSEE.BM.RAW)
write.csv(CSEE.BM.RAW, file="CSEE_BM_RAW.csv")


#load new raw bm data
CSEE.BM.RAW<-read.csv("CSEE_BM_RAW.csv", header=T)
name<-CSEE.BM.RAW$name
Order<-CSEE.BM.RAW$Order
Species<-CSEE.BM.RAW$Species
MGMT<-CSEE.BM.RAW$MGMT
DBM<-CSEE.BM.RAW$DBM
PDBM<-CSEE.BM.RAW$PDBM
Age<-CSEE.BM.RAW$Age


summary(lm(DBM~Order))
summary.aov(lm(DBM~Order))
boxplot(DBM~Order)

summary(lm(DBM~MGMT))
summary.aov(lm(DBM~MGMT))
boxplot(DBM~MGMT)
summary(MGMT)

###################################################
###                                             ###
###         CSEE RAW BM Limited Species         ###
###                                             ###
###                                             ###
###################################################

cod<-subset(CSEE.BM.RAW, CSEE.BM.RAW$Species == "morhua")
had<-subset(CSEE.BM.RAW, CSEE.BM.RAW$Species=="aeglefinus")
her<-subset(CSEE.BM.RAW, CSEE.BM.RAW$Species=="harengus")
sole<-subset(CSEE.BM.RAW, CSEE.BM.RAW$Species=="vulgaris")
plai<-subset(CSEE.BM.RAW, CSEE.BM.RAW$Species=="platessa")

cod<-as.data.frame(cod)
had<-as.data.frame(had)
her<-as.data.frame(her)
sole<-as.data.frame(sole)
plai<-as.data.frame(plai)

SPRAW <- rbind(cod,had,her,sole,plai)
write.csv(SPRAW, file="SPPBMRAW.csv")
summary(SPRAW)

#gives 26 morhua, 24 aegle, 22 hareng, 14 vulg, 10 plat /2
# load new dataframe containing only those species

rm(list=ls(all=T))
SPLIM<-read.csv("SPPBMRAW.csv", header=T)
summary(SPLIM)
DBM<-SPLIM$DBM
PDNUM<-SPLIM$PDBM
Order<-SPLIM$Order
name<-SPLIM$name
Species<-SPLIM$Species
MGMT<-SPLIM$MGMT
Age<-SPLIM$Age

summary(lm(DBM~Species))
summary.aov(lm(DBM~Species))
#Graph 7
boxplot(DBM~Species, main="Difference in Biomass for Selected Species", 
        col=(c("lightpink", "salmon", "lightpink","palevioletred4", "palevioletred4")))

summary(lm(DBM~Age*Species))
summary(lmer(DBM~Age+(1+Age|Species)))
summary.aov(lm(DBM~Age*Species))

#### backtransform logs in r
log(10000)
exp(9.21034)

log(93282.83)
exp(11.44339)

###################################################
###                                             ###
###               NUMBER LOGGED                 ###
###                  Numberrs                   ###
###                                             ###
###################################################

rm(list=ls(all=T))
library(matlab)
library(lme4)
library(arm)

ASD_MGMT<-read.csv("ASD_MGMT.csv", header=T)
names(ASD_MGMT)
ASD_MGMT$Stock.ID
levels(ASD_MGMT$Management)
ASD_MGMT$Species
geo.mean <- function(x,n) prod(x)^(1/n) 

unique.stocks <- as.character(unique(ASD_MGMT$Stock.ID))
num.stocks <- length(unique.stocks)
names(ASD_MGMT)

CSEE <- zeros(num.stocks,8)

for(i in 1:num.stocks)  
{
  if(i != 69 && i != 72 && i != 125 && i != 19 && i !=18 && i != 9 && i != 29 && i !=41)
  {
    name <- unique.stocks[i]
    Order<-as.character(ASD_MGMT$Order[ASD_MGMT$Stock.ID == name ][1])
    MGMT<-as.character(ASD_MGMT$Management[ASD_MGMT$Stock.ID == name ] [1])
    Species<-as.character(ASD_MGMT$Species[ASD_MGMT$Stock.ID == name] [1])
    st <- ASD_MGMT[ASD_MGMT$Stock.ID == name,grep("Num",names(ASD_MGMT))] 
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
      NUM_LQ_start <- st[1:years,LQ+1]
      NUM_UQ_start <-  st[1:years,UQ+1]
      NUM_LQ_now <- st[(ts_len-years+1):ts_len,LQ+1]
      NUM_UQ_now <- st[(ts_len-years+1):ts_len,UQ+1]
      NUM_LQ_start_tot <- rowSums(NUM_LQ_start)
      NUM_UQ_start_tot <- rowSums(NUM_UQ_start)
      NUM_LQ_now_tot <- rowSums(NUM_LQ_now)
      NUM_UQ_now_tot <- rowSums(NUM_UQ_now)
      NUM_geo_LQS <- log(geo.mean(NUM_LQ_start_tot,years)) #geomean of youngest over first 5 years
      NUM_geo_UQS <- log(geo.mean(NUM_UQ_start_tot,years)) #geomean of oldest over first 5 years
      NUM_geo_LQN <- log(geo.mean(NUM_LQ_now_tot,years)) #geomean of youngest over last 5 years (now)
      NUM_geo_UQN <- log(geo.mean(NUM_UQ_now_tot,years)) #geomean of oldest over last 5 years (now)
      NUM_dec_LQ <- (NUM_geo_LQN - NUM_geo_LQS) 
      NUM_dec_UQ <- (NUM_geo_UQN - NUM_geo_UQS) 
      percent_difference_old<-((NUM_geo_UQN-NUM_geo_UQS)/NUM_geo_UQS)*100
      percent_difference_young<-((NUM_geo_LQN-NUM_geo_LQS)/NUM_geo_LQS)*100 #geometric mean
      CSEE[i,]<-c(name, Order, Species, MGMT, NUM_dec_UQ, percent_difference_old,NUM_dec_LQ, percent_difference_young)
      }
  }
}
i
CSEE[,]
bad.stocks <- c(which(CSEE[,1] == 0), attr(na.omit(CSEE[,5]),"na.action"))
CSEE <- CSEE[-bad.stocks,]

# save(GEO,file="GEO.rdata")
CSEE.ALL <- rbind(CSEE[,1:6],CSEE[,c(1,2,3,4,7,8)]) #divide young and old
CSEE.ALL <- as.data.frame(CSEE.ALL,stringsAsFactors=F)
CSEE.ALL$Age <- c(rep("Old",length(CSEE.ALL[,1])/2),rep("Young",length(CSEE.ALL[,1])/2))
colnames(CSEE.ALL) <- c("name", "Order", "Species", "MGMT", "DNUM", "PDNUM", "Age")
str(CSEE.ALL)
CSEE.ALL$DNUM <- as.numeric(CSEE.ALL$DNUM)
CSEE.ALL$PDNUM <- as.numeric(CSEE.ALL$PDNUM)
qqnorm(CSEE.ALL$DNUM)
qqline(CSEE.ALL$DNUM)

oldflat<-CSEE.ALL[CSEE.ALL$Order == "Pleuronectiformes" & CSEE.ALL$Age == "Old",]
z<-sum(-oldflat$DNUM)
exp(z)

mean(oldflat$DNUM)
median(oldflat$DNUM)
youngflat<-CSEE.ALL[CSEE.ALL$Order == "Pleuronectiformes" & CSEE.ALL$Age == "Young",]
mean(youngflat$DNUM)
median(youngflat$DNUM)
oldgad<-CSEE.ALL[CSEE.ALL$Order == "Gadiformes" & CSEE.ALL$Age == "Old",]
mean(oldgad$DNUM)
median(oldgad$DNUM)
younggad<-CSEE.ALL[CSEE.ALL$Order == "Gadiformes" & CSEE.ALL$Age == "Young",]
mean(younggad$DNUM, na.rm=T)
median(younggad$DNUM, na.rm=T)
oldherr<-CSEE.ALL[CSEE.ALL$Order == "Clupeiformes" & CSEE.ALL$Age == "Old",]
mean(oldherr$DNUM)
median(oldherr$DNUM)
youngherr<-CSEE.ALL[CSEE.ALL$Order == "Clupeiformes" & CSEE.ALL$Age == "Young",]
mean(youngherr$DNUM)
median(youngherr$DNUM)
oldperc<-CSEE.ALL[CSEE.ALL$Order == "Perciformes" & CSEE.ALL$Age == "Old",]
mean(oldperc$DNUM)
median(oldperc$DNUM)
youngperc<-CSEE.ALL[CSEE.ALL$Order == "Perciformes" & CSEE.ALL$Age == "Young",]
mean(youngperc$DNUM)
median(youngperc$DNUM)



oldflat<-CSEE.ALL[CSEE.ALL$Order == "Pleuronectiformes" & CSEE.ALL$Age == "Old",]
mean(oldflat$PDNUM)
median(oldflat$PDNUM)
youngflat<-CSEE.ALL[CSEE.ALL$Order == "Pleuronectiformes" & CSEE.ALL$Age == "Young",]
mean(youngflat$PDNUM)
median(youngflat$PDNUM)
oldgad<-CSEE.ALL[CSEE.ALL$Order == "Gadiformes" & CSEE.ALL$Age == "Old",]
mean(oldgad$PDNUM)
median(oldgad$PDNUM)
younggad<-CSEE.ALL[CSEE.ALL$Order == "Gadiformes" & CSEE.ALL$Age == "Young",]
mean(younggad$PDNUM, na.rm=T)
median(younggad$PDNUM, na.rm=T)
oldherr<-CSEE.ALL[CSEE.ALL$Order == "Clupeiformes" & CSEE.ALL$Age == "Old",]
mean(oldherr$PDNUM)
median(oldherr$PDNUM)
youngherr<-CSEE.ALL[CSEE.ALL$Order == "Clupeiformes" & CSEE.ALL$Age == "Young",]
mean(youngherr$PDNUM)
median(youngherr$PDNUM)
oldperc<-CSEE.ALL[CSEE.ALL$Order == "Perciformes" & CSEE.ALL$Age == "Old",]
mean(oldperc$PDNUM)
median(oldperc$PDNUM)
youngperc<-CSEE.ALL[CSEE.ALL$Order == "Perciformes" & CSEE.ALL$Age == "Young",]
mean(youngperc$PDNUM)
median(youngperc$PDNUM)


# so i will remove the scorpaeniformes, bc it may change the model
#REMOVE SCORPAENIFORMES
CSEE.ALL$Order
rm <- c(which(CSEE.ALL$Order == "Scorpaeniformes")) #eerst samenvoegen
CSEE.ALL<- CSEE.ALL[-rm,] #dan eruit halen
CSEE.ALL$Order
CSEE.ALL$MGMT<- as.factor(CSEE.ALL$MGMT)
CSEE.ALL$Age <- as.factor(CSEE.ALL$Age)
CSEE.ALL$Order <- as.factor(CSEE.ALL$Order)
summary(CSEE.ALL)
CSEE.ALL[,1]
# remove NAFO, ICCAT #
rm <- c(which(CSEE.ALL$MGMT == "NAFO"))
CSEE.ALL<- CSEE.ALL[-rm,]
CSEE.ALL$MGMT
summary(CSEE.ALL)

rm <- c(which(CSEE.ALL$MGMT == "ICCAT"))
CSEE.ALL<- CSEE.ALL[-rm,]
CSEE.ALL$MGMT
summary(CSEE.ALL)

#to make the stocks that I use for number equal to the stocks I use for BM
#i need to take out 7 stocks:
# "ALPLAICBSAIf", "ALPLAICBSAIm", "AMPL4T", "PCODBSAI", "PCODGA", 
#"POLL4VWX5Zc","WPOLLGA"

rmALF<-c(which(CSEE.ALL$name == "ALPLAICBSAIf"))
CSEE.ALL<- CSEE.ALL[-rmALF,]
rmALM<-c(which(CSEE.ALL$name =="ALPLAICBSAIm"))
CSEE.ALL<- CSEE.ALL[-rmALM,]
rmAM<-c(which(CSEE.ALL$name =="AMPL4T"))
CSEE.ALL<- CSEE.ALL[-rmAM,]
rmPCI<-c(which(CSEE.ALL$name =="PCODBSAI"))
CSEE.ALL<- CSEE.ALL[-rmPCI,]
rmPCA<-c(which(CSEE.ALL$name =="PCODGA"))
CSEE.ALL<- CSEE.ALL[-rmPCA,]
rmPO<-c(which(CSEE.ALL$name =="POLL4VWX5Zc"))
CSEE.ALL<- CSEE.ALL[-rmPO,]
rmWP<-c(which(CSEE.ALL$name =="WPOLLGA"))
CSEE.ALL<- CSEE.ALL[-rmWP,]
summary(CSEE.ALL)
CSEE.ALL[,1]
summary(CSEE.ALL)
write.csv(CSEE.ALL, file="CSEE_ALL_NUM.csv")

# models #

rm(list=ls(all=T))
CSEE_ALL_NUM<-read.csv("CSEE_ALL_NUM.csv", header=T)
# not necessary, no zeroes bad<- c(which(CSEE_ALL_NUM[,7] == 0), attr(na.omit(CSEE_ALL_NUM[,7]),"na.action"))
#CSEE_ALL_NUM <- CSEE_ALL_NUM[-bad,]
#CSEE_ALL_NUM

name<-CSEE_ALL_NUM$name
Order<-CSEE_ALL_NUM$Order
Species<-CSEE_ALL_NUM$Species
MGMT<-CSEE_ALL_NUM$MGMT
DNUM<-CSEE_ALL_NUM$DNUM
PDNUM<-CSEE_ALL_NUM$PDNUM
Age<-CSEE_ALL_NUM$Age


# lets peak at the data.."name", "Order", "Species", "MGMT", "DNUM", "PDNUM", "Age")
#GRAPH N1
lm.mod.1 <- lm(DNUM~Age,CSEE_ALL_NUM)
summary(lm.mod.1)
anova(lm.mod.1)
summary.aov(lm.mod.1)
boxplot(DNUM~Age,main="Difference in Log of Number * 1000", col=(c("blue3", "lightblue")))

#FANCY
dna<-(lm(DNUM~Age-1))
coef(dna)
ci_dna<-confint(dna,level=.95)
lci_dna<-ci_dna[,1]
uci_dna<-ci_dna[,2]
mn_dna<- coef(summary.lm(dna))[,1]
elci_dna<-(exp(lci_dna)-1)*100
euci_dna<-(exp(uci_dna)-1)*100
emn_dna<-(exp(mn_dna)-1)*100
par(mar=c(5,5,5,13),cex=1)
#old
plot(emn_dna[c(1)],1.5, ylim=c(1,length(emn_dna)+1),
     bty="U",pch=19,xlim=c(min(elci_dna),max(euci_dna)),ylab="",yaxt="n",xlab="", cex=1)
#young
points(emn_dna[c(2)],2.5,pch=19, font=2)
# Now I add the CI's
segments(elci_dna[c(1)],1.5,euci_dna[c(1)],1.5)
segments(elci_dna[c(2)],2.5,euci_dna[c(2)],2.5)
# Throw a line in to divide them
abline(h=c(2),lty=2,lwd=0.5,col="blue")
# Add some text to the margins on the next 3 lines...
axis(2,at=c(1.5,2.5),labels=c("Old","Young"),las=1,tck=0, cex.axis=1.2)
mtext(side=1,"Loss of Numbers in %",line=3, cex=1.2)

#Graph N2
lm.mod.2 <-lm(DNUM~MGMT, CSEE_ALL_NUM)
summary(lm.mod.2)
summary.aov(lm.mod.2)
boxplot(DNUM~MGMT,main="Difference in Log of Number * 1000", col=(c("coral","seagreen", "aliceblue","cadetblue")))
#plot(lm.mod.2)

#graph 2 fancy
summary(MGMT)
dm<-(lm(DNUM~MGMT-1))
coef(dm)
ci_dm<-confint(dm,level=.95)

lci_dm<-ci_dm[,1]
uci_dm<-ci_dm[,2]
mn_dm<- coef(summary.lm(dm))[,1]

elci_dm<-(exp(lci_dm)-1)*100
euci_dm<-(exp(uci_dm)-1)*100
emn_dm<-(exp(mn_dm)-1)*100

par(mar=c(5,5,5,13),cex=1)
#mgmt
plot(emn_dm[c(1:4)],1.5:4.5, ylim=c(1,length(emn_dm)+1),
     bty="U",pch=19,xlim=c(min(elci_dm),max(euci_dm)),ylab="",yaxt="n",xlab="", cex=1)
#points(mn_BM_I_peak[c(5:8)],6:9,pch=19)
# Now I add the CI's
segments(elci_dm[c(1:4)],1.5:4.5,euci_dm[c(1:4)],1.5:4.5)
# Throw a line in to divide them
abline(h=c(2,3,4),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
mtext(side=1,"Loss of Numbers in %",line=3, cex=1.2)
axis(4,at=c(1.5,2.5,3.5,4.5),labels=c("DFO (8)","DFO+NOAA (6)","ICES (104)","NOAA (64)"),las=1, cex.axis=1.2)


#Graph N3
lm.mod.3 <-lm(DNUM~Order, CSEE_ALL_NUM)
summary(lm.mod.3)
summary.aov(lm.mod.3)
boxplot(DNUM~Order, main="Difference in Log of Number * 1000", col=(c("salmon", "lightpink", "hotpink", "palevioletred4")))

dno<-(lm(DNUM~Order-1))
coef(dno)
ci_dno<-confint(dno,level=.95)

lci_dno<-ci_dno[,1]
uci_dno<-ci_dno[,2]
mn_dno<- coef(summary.lm(dno))[,1]

elci_dno<-(exp(lci_dno)-1)*100
euci_dno<-(exp(uci_dno)-1)*100
emn_dno<-(exp(mn_dno)-1)*100

par(mar=c(5,5,5,13),cex=1)
#order
plot(emn_dno[c(1:4)],1.5:4.5, ylim=c(1,length(emn_dno)+1),
     bty="U",pch=19,xlim=c(min(elci_dno),max(euci_dno)),ylab="",yaxt="n",xlab="", cex=1)
#points(mn_BM_I_peak[c(5:8)],6:9,pch=19)
# Now I add the CI's
segments(elci_dno[c(1:4)],1.5:4.5,euci_dno[c(1:4)],1.5:4.5)
# Throw a line in to divide them
abline(h=c(2,3,4),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
mtext(side=1,"Loss of Numbers in %",line=3, cex=1.2)
axis(4,at=c(1.5,2.5,3.5,4.5),labels=c("Clupeiformes","Gadiformes","Perciformes","Pleuronectiformes"),las=1, cex.axis=1.2)
as.list(emn_dno)


#GRAPH N4
lm.mod.4 <-lm(DNUM~Age*Order)
summary(lm.mod.4)
summary.aov(lm.mod.4) #or use anova(lm.mod.4) gives same result
boxplot(DNUM~Age*Order, main="Difference in Log Number * 1000",  col=(c("blue3", "lightblue")))

DNAO<-(lm(DNUM~Age*Order-Age-Order-1))
coef(DNAO)
ci_DNAO<-confint(DNAO,level=.95)

lci_DNAO<-ci_DNAO[,1]
uci_DNAO<-ci_DNAO[,2]
mn_DNAO<- coef(summary.lm(DNAO))[,1]

elci_DNAO<-(exp(lci_DNAO)-1)*100
euci_DNAO<-(exp(uci_DNAO)-1)*100
emn_DNAO<-(exp(mn_DNAO)-1)*100


# And now I do the plot..
par(mar=c(5,5,5,13),cex=1)
#old
OLD<-c(1,3,5,7)
plot(emn_DNAO[c(1,3,5,7)],1:4, ylim=c(1,length(emn_DNAO)+1),
     bty="U",pch=19,xlim=c(min(elci_DNAO),max(euci_DNAO)),ylab="",yaxt="n",xlab="", cex=1)
#young
YOUNG<-c(2,4,6,8)
points(emn_DNAO[c(2,4,6,8)],6:9,pch=19, cex=1)
# Now I add the CI's
segments(elci_DNAO[c(1,3,5,7)],1:4,euci_DNAO[c(1,3,5,7)],1:4)
segments(elci_DNAO[c(2,4,6,8)],6:9,euci_DNAO[c(2,4,6,8)],6:9)
# Throw a line in to divide them
abline(h=c(5),lty=2,lwd=0.5,col="blue")
# Add some text to the margins on the next 3 lines...
axis(2,at=c(2.5,7.5,11.5),labels=c("Old","Young",""),las=1,tck=0, cex.axis=1.2)
mtext(side=1,"Loss of Numbers in %",line=3, cex=1.2)
axis(4,at=c(1,2,3,4,6,7,8,9),labels=c("Clupeiformes","Gadiformes","Perciformes","Pleuronectiformes",
                                      "Clupeiformes","Gadiformes","Perciformes","Pleuronectiformes"),las=1, cex.axis=1.2)
as.list(emn_DNAO)


#lm.mod.a <-lm(DNUM~Age*MGMT,CSEE_ALL_NUM)
#summary(lm.mod.a)
#best.mod <- step(lm.mod.a)
#anova(lm.mod.a)
#summary.aov(lm.mod.a)
#summary(Species)
#boxplot(DNUM~Age*Species)

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

#mod.NUM.log <- lmer(DNUM ~ Age + (1+Age|Order),CSEE.ALL) 
#summary(mod.NUM.log)
#anova(mod.NUM.log)
#ranef(mod.NUM.log)
#fixef(mod.NUM.log)
#tot<-(ranef(mod.NUM.log))$Order + fixef(mod.NUM.log)
#tot
#qqnorm(CSEE.ALL$PDNUM)
#qqline(CSEE.ALL$PDNUM)

##################################################################
###                                                            ###
###             Log Number for selected species                ###
###                           SPNumberrlogg                    ###
###                                                            ###
##################################################################

#using CSEE_ALL_NUM.csv

cod<-subset(CSEE_ALL_NUM,CSEE_ALL_NUM$Species == "morhua")
had<-subset(CSEE_ALL_NUM,CSEE_ALL_NUM$Species=="aeglefinus")
her<-subset(CSEE_ALL_NUM,CSEE_ALL_NUM$Species=="harengus")
sole<-subset(CSEE_ALL_NUM,CSEE_ALL_NUM$Species=="vulgaris")
plai<-subset(CSEE_ALL_NUM,CSEE_ALL_NUM$Species=="platessa")

cod<-as.data.frame(cod)
had<-as.data.frame(had)
her<-as.data.frame(her)
sole<-as.data.frame(sole)
plai<-as.data.frame(plai)

SPPNUMLOG <- rbind(cod,had,her,sole,plai)
write.csv(SPPNUMLOG, file="SPPNUMLOG.csv")

#load new dataframe containing only those species

rm(list=ls(all=T))
SPLIM<-read.csv("SPPNUMLOG.csv", header=T)
summary(SPLIM)

SPLIMDNUM<-SPLIM$DNUM
SPLIMPDNUM<-SPLIM$PDNUM
SPLIMOrder<-SPLIM$Order
SPLIMname<-SPLIM$name
SPLIMSpecies<-SPLIM$Species
SPLIMMGMT<-SPLIM$MGMT
SPLIMAGE<-SPLIM$Age

#GRAPH N5
summary(lm(SPLIMDNUM~SPLIMSpecies))
summary.aov(lm(SPLIMDNUM~SPLIMSpecies))
boxplot(SPLIMDNUM~SPLIMSpecies, main="Difference in Log of Number * 1000 for Selected Species", 
        col=(c("lightpink", "salmon", "lightpink","palevioletred4", "palevioletred4")))

#GRAPH N5
DNS<-(lm(SPLIMDNUM~SPLIMSpecies-1))
coef(DNS)
ci_DNS<-confint(DNS,level=.95)

lci_DNS<-ci_DNS[,1]
uci_DNS<-ci_DNS[,2]
mn_DNS<- coef(summary.lm(DNS))[,1]

elci_DNS<-(exp(lci_DNS)-1)*100
euci_DNS<-(exp(uci_DNS)-1)*100
emn_DNS<-(exp(mn_DNS)-1)*100


par(mar=c(5,5,5,13),cex=1)
#old
plot(emn_DNS[c(1:5)],1.5:5.5, ylim=c(1,length(emn_DNS)+1),
     bty="U",pch=19,xlim=c(min(elci_DNS),max(euci_DNS)),ylab="",yaxt="n",xlab="", cex=1)
#points(mn_BM_I_peak[c(5:8)],6:9,pch=19)
# Now I add the CI's
segments(elci_DNS[c(1:5)],1.5:5.5,euci_DNS[c(1:5)],1.5:5.5)
# Throw a line in to divide them
abline(h=c(2,3,4,5),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
mtext(side=1,"Loss of Numbers in %",line=3, cex=1.2)
axis(4,at=c(1.5:5.5),labels=c("aeglefinus (12)","harengus (11)","morhua (13)","platessa (5)",
                              "vulgaris (7)"),las=1, cex.axis=1.1)

as.list(emn_DNS)


#GRAPH N6
summary(lm(SPLIMDNUM~SPLIMAGE*SPLIMSpecies))
summary.aov(lm(SPLIMDNUM~SPLIMAGE*SPLIMSpecies))
boxplot(SPLIMDNUM~SPLIMAGE*SPLIMSpecies, main="Difference in Log of Number * 1000 for Selected Species"
        , col=(c("blue3", "lightblue")))

das<-(lm(SPLIMDNUM~SPLIMAGE*SPLIMSpecies-SPLIMAGE-SPLIMSpecies-1))
coef(das)
ci_das<-confint(das,level=.95)
lci_das<-ci_das[,1]
uci_das<-ci_das[,2]
mn_das<- coef(summary.lm(das))[,1]
elci_das<-(exp(lci_das)-1)*100
euci_das<-(exp(uci_das)-1)*100
emn_das<-(exp(mn_das)-1)*100
# And now I do the plot..
par(mar=c(5,5,5,13),cex=1)
# So first I pull out of the mn_BM_I_peak / mn_das object the first 4 taxonomic groups
#(you'll have to customize yourself) using the c(1:4) from the DI category
#old
plot(emn_das[c(1,3,5,7,9)],1:5, ylim=c(1,length(emn_das)+1),
     bty="U",pch=19,xlim=c(min(elci_das),max(euci_das)),ylab="",yaxt="n",xlab="", cex=1)
#plot(mn_BM_I_peak[c(1:4)],1:4,ylim=c(1,length(mn_BM_I_peak)+2),
#    bty="U",pch=19,xlim=c(min(LCI_BMI_peak),max(UCI_BMI_peak)),ylab="",yaxt="n",xlab="") 
# Then I get the ones that are NDD / age y
#young
points(emn_das[c(2,4,6,8,10)],6:10,pch=19, cex=1)
#ci
segments(elci_das[c(1,3,5,7,9)],1:5,euci_das[c(1,3,5,7,9)],1:5)
segments(elci_das[c(2,4,6,8,10)],6:10,euci_das[c(2,4,6,8,10)],6:10)
# Throw a line in to divide them
abline(h=c(5.5),lty=2,lwd=0.5,col="blue")
#abline(v=0,lty=2,lwd=0.5,col="grey")
# Add some text to the margins on the next 3 lines...
axis(2,at=c(3,8,11),labels=c("Old","Young",""),las=1,tck=0, cex.axis=1.2)
mtext(side=1,"Loss of Numbers in %",line=3, cex=1.2)
axis(4,at=c(1:10),labels=c("aeglefinus","harengus","morhua","platessa","vulgaris",
                           "aeglefinus","harengus","morhua","platessa","vulgaris"),las=1, cex.axis=1)

as.list(emn_das)

boxplot(SPLIMDNUM~SPLIMAGE)
boxplot(SPLIMDNUM~SPLIMSpecies)
boxplot(SPLIMDNUM~SPLIMSpecies*SPLIMMGMT)
boxplot(SPLIMDNUM~SPLIMOrder*SPLIMAGE)
boxplot(SPLIMDNUM~SPLIMOrder)
boxplot(SPLIMDNUM~SPLIMOrder*SPLIMMGMT*SPLIMAGE)
boxplot(SPLIMDNUM~SPLIMSpecies*SPLIMMGMT)
boxplot(SPLIMDNUM~SPLIMAGE*SPLIMSpecies)
boxplot(SPLIMPDNUM~SPLIMSpecies)


# remove outlier
#rm(list=ls(all=T))
#SPLIM<-read.csv("SPPoutlierremoved.csv", header=T)
#summary(SPLIM)

##################################################################
###                                                            ###
###                         RAW  Number                        ###
###                         CSEENUMRRAW                        ###
###                                                            ###
##################################################################


rm(list=ls(all=T))
library(matlab)
library(lme4)
library(arm)
detach("package:nlme")

setwd("C:/Users/Emilie/Documents/CSEE")

ASD_MGMT<-read.csv("ASD_MGMT.csv", header=T)
names(ASD_MGMT)
ASD_MGMT$Stock.ID
levels(ASD_MGMT$Management)
ASD_MGMT$Species
geo.mean <- function(x,n) prod(x)^(1/n) 

unique.stocks <- as.character(unique(ASD_MGMT$Stock.ID))
num.stocks <- length(unique.stocks)
names(ASD_MGMT)

CSEENUMRAW <- zeros(num.stocks,8)

for(i in 1:num.stocks)  
{
  if(i != 69 && i != 72 && i != 125 && i != 19 && i !=18 && i != 9 && i != 29 && i !=41)
  {
    name <- unique.stocks[i]
    Order<-as.character(ASD_MGMT$Order[ASD_MGMT$Stock.ID == name ][1])
    MGMT<-as.character(ASD_MGMT$Management[ASD_MGMT$Stock.ID == name ] [1])
    Species<-as.character(ASD_MGMT$Species[ASD_MGMT$Stock.ID == name] [1])
    st <- ASD_MGMT[ASD_MGMT$Stock.ID == name,grep("Num",names(ASD_MGMT))] 
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
      NUM_LQ_start <- st[1:years,LQ+1]
      NUM_UQ_start <-  st[1:years,UQ+1]
      NUM_LQ_now <- st[(ts_len-years+1):ts_len,LQ+1]
      NUM_UQ_now <- st[(ts_len-years+1):ts_len,UQ+1]
      NUM_LQ_start_tot <- rowSums(NUM_LQ_start)
      NUM_UQ_start_tot <- rowSums(NUM_UQ_start)
      NUM_LQ_now_tot <- rowSums(NUM_LQ_now)
      NUM_UQ_now_tot <- rowSums(NUM_UQ_now)
      NUM_geo_LQS <- geo.mean(NUM_LQ_start_tot,years) #geomean of youngest over first 5 years
      NUM_geo_UQS <- geo.mean(NUM_UQ_start_tot,years) #geomean of oldest over first 5 years
      NUM_geo_LQN <- geo.mean(NUM_LQ_now_tot,years) #geomean of youngest over last 5 years (now)
      NUM_geo_UQN <- geo.mean(NUM_UQ_now_tot,years) #geomean of oldest over last 5 years (now)
      NUM_dec_LQ <- (NUM_geo_LQN - NUM_geo_LQS) 
      NUM_dec_UQ <- (NUM_geo_UQN - NUM_geo_UQS) 
      percent_difference_old<-((NUM_geo_UQN-NUM_geo_UQS)/NUM_geo_UQS)*100
      percent_difference_young<-((NUM_geo_LQN-NUM_geo_LQS)/NUM_geo_LQS)*100 #geometric mean
      CSEENUMRAW[i,]<-c(name, Order, Species, MGMT, NUM_dec_UQ, percent_difference_old,NUM_dec_LQ, percent_difference_young)
    }
  }
}

CSEENUMRAW[,]
summary(CSEENUMRAW)
bad.stocks <- c(which(CSEENUMRAW[,6] == 0), attr(na.omit(CSEENUMRAW[,5]),"na.action"))
CSEENUMRAW <- CSEENUMRAW[-bad.stocks,]
CSEENUMRAW <- rbind(CSEENUMRAW[,1:6],CSEENUMRAW[,c(1,2,3,4,7,8)]) #divide young and old
CSEENUMRAW <- as.data.frame(CSEENUMRAW,stringsAsFactors=F)
CSEENUMRAW$Age <- c(rep("Old",length(CSEENUMRAW[,1])/2),rep("Young",length(CSEENUMRAW[,1])/2))
colnames(CSEENUMRAW) <- c("name", "Order", "Species", "MGMT", "DNUM", "PDNUM", "Age")


rmPO<-c(which(CSEENUMRAW$name =="POLL4VWX5Zc"))
CSEENUMRAW<- CSEENUMRAW[-rmPO,]

summary(CSEENUMRAW)
CSEENUMRAW$DNUM <- as.numeric(CSEENUMRAW$DNUM)
CSEENUMRAW$PDNUM <- as.numeric(CSEENUMRAW$PDNUM)
qqnorm(CSEENUMRAW$DNUM)
qqline(CSEENUMRAW$DNUM)
summary(CSEENUMRAW)


oldflat<-CSEENUMRAW[CSEENUMRAW$Order == "Pleuronectiformes" & CSEENUMRAW$Age == "Old",]
summary(oldflat)
mean(oldflat$DNUM)
youngflat<-CSEENUMRAW[CSEENUMRAW$Order == "Pleuronectiformes" & CSEENUMRAW$Age == "Young",]
mean(youngflat$DNUM)
summary(youngflat)
oldgad<-CSEENUMRAW[CSEENUMRAW$Order == "Gadiformes" & CSEENUMRAW$Age == "Old",]
mean(oldgad$DNUM)
median(oldgad$DNUM)
summary(oldgad)

younggad<-CSEENUMRAW[CSEENUMRAW$Order == "Gadiformes" & CSEENUMRAW$Age == "Young",]
mean(younggad$DNUM, na.rm=T)
median(younggad$DNUM, na.rm=T)
summary(younggad)

oldherr<-CSEENUMRAW[CSEENUMRAW$Order == "Clupeiformes" & CSEENUMRAW$Age == "Old",]
mean(oldherr$DNUM)
median(oldherr$DNUM)
summary(oldherr)

youngherr<-CSEENUMRAW[CSEENUMRAW$Order == "Clupeiformes" & CSEENUMRAW$Age == "Young",]
mean(youngherr$DNUM)
median(youngherr$DNUM)
summary(youngherr)

oldperc<-CSEENUMRAW[CSEENUMRAW$Order == "Perciformes" & CSEENUMRAW$Age == "Old",]
mean(oldperc$DNUM)
median(oldperc$DNUM)
summary(oldperc)

youngperc<-CSEENUMRAW[CSEENUMRAW$Order == "Perciformes" & CSEENUMRAW$Age == "Young",]
mean(youngperc$DNUM)
median(youngperc$DNUM)
summary(youngperc)

#percentage difference numbers raw(nonlogged)
oldflat<-CSEENUMRAW[CSEENUMRAW$Order == "Pleuronectiformes" & CSEENUMRAW$Age == "Old",]
mean(oldflat$PDNUM)
median(oldflat$PDNUM)
youngflat<-CSEENUMRAW[CSEENUMRAW$Order == "Pleuronectiformes" & CSEENUMRAW$Age == "Young",]
mean(youngflat$PDNUM)
median(youngflat$PDNUM)
oldgad<-CSEENUMRAW[CSEENUMRAW$Order == "Gadiformes" & CSEENUMRAW$Age == "Old",]
mean(oldgad$PDNUM)
median(oldgad$PDNUM)
younggad<-CSEENUMRAW[CSEENUMRAW$Order == "Gadiformes" & CSEENUMRAW$Age == "Young",]
mean(younggad$PDNUM, na.rm=T)
median(younggad$PDNUM, na.rm=T)
oldherr<-CSEENUMRAW[CSEENUMRAW$Order == "Clupeiformes" & CSEENUMRAW$Age == "Old",]
mean(oldherr$PDNUM)
median(oldherr$PDNUM)
youngherr<-CSEENUMRAW[CSEENUMRAW$Order == "Clupeiformes" & CSEENUMRAW$Age == "Young",]
mean(youngherr$PDNUM)
median(youngherr$PDNUM)
oldperc<-CSEENUMRAW[CSEENUMRAW$Order == "Perciformes" & CSEENUMRAW$Age == "Old",]
mean(oldperc$PDNUM)
median(oldperc$PDNUM)
youngperc<-CSEENUMRAW[CSEENUMRAW$Order == "Perciformes" & CSEENUMRAW$Age == "Young",]
mean(youngperc$PDNUM)
median(youngperc$PDNUM)

summary(CSEENUMRAW)

#remove scorp, iccat etc then
#REMOVE SCORPAENIFORMES
summary(CSEENUMRAW)
rm <- c(which(CSEENUMRAW$Order == "Scorpaeniformes")) #eerst samenvoegen
CSEENUMRAW<- CSEENUMRAW[-rm,] #dan eruit halen
CSEENUMRAW$Order

CSEENUMRAW$MGMT<- as.factor(CSEENUMRAW$MGMT)
CSEENUMRAW$Age <- as.factor(CSEENUMRAW$Age)
CSEENUMRAW$Order <- as.factor(CSEENUMRAW$Order)


# remove NAFO, ICCAT #
rm <- c(which(CSEENUMRAW$MGMT == "NAFO"))
CSEENUMRAW<- CSEENUMRAW[-rm,]
CSEENUMRAW$MGMT
summary(CSEENUMRAW)


#to make the stocks that I use for number equal to the stocks I use for BM
#i need to take out 7 stocks:
# "ALPLAICBSAIf", "ALPLAICBSAIm", "AMPL4T", "PCODBSAI", "PCODGA", 
#"POLL4VWX5Zc","WPOLLGA"

rmALF<-c(which(CSEENUMRAW$name == "ALPLAICBSAIf"))
CSEENUMRAW<- CSEENUMRAW[-rmALF,]
rmALM<-c(which(CSEENUMRAW$name =="ALPLAICBSAIm"))
CSEENUMRAW<- CSEENUMRAW[-rmALM,]
rmAM<-c(which(CSEENUMRAW$name =="AMPL4T"))
CSEENUMRAW<- CSEENUMRAW[-rmAM,]
rmPCI<-c(which(CSEENUMRAW$name =="PCODBSAI"))
CSEENUMRAW<- CSEENUMRAW[-rmPCI,]
rmPCA<-c(which(CSEENUMRAW$name =="PCODGA"))
CSEENUMRAW<- CSEENUMRAW[-rmPCA,]
#rmPO<-c(which(CSEENUMRAW$name =="POLL4VWX5Zc"))
#CSEENUMRAW<- CSEENUMRAW[-rmPO,]
rmWP<-c(which(CSEENUMRAW$name =="WPOLLGA"))
CSEENUMRAW<- CSEENUMRAW[-rmWP,]
summary(CSEENUMRAW)
CSEENUMRAW[,1]
summary(CSEENUMRAW)

write.csv(CSEENUMRAW, file="CSEENUMRAW.csv")


###################################################
###                                             ###
###         CSEE RAW NUM Limited Species        ###
###                       NUMRAWSPLIM           ###
###                                             ###
###################################################

#reload the limited species data
rm(list=ls(all=T))
CSEENUMRAW<-read.csv("CSEENUMRAW.csv", header=T)
summary(CSEENUMRAW)

name<-CSEENUMRAW$name
Order<-CSEENUMRAW$Order
Species<-CSEENUMRAW$Species
MGMT<-CSEENUMRAW$MGMT
DNUM<-CSEENUMRAW$DNUM
PDNUM<-CSEENUMRAW$PNUM
Age<-CSEENUMRAW$Age


cod<-subset(CSEENUMRAW, CSEENUMRAW$Species == "morhua")
had<-subset(CSEENUMRAW, CSEENUMRAW$Species=="aeglefinus")
her<-subset(CSEENUMRAW, CSEENUMRAW$Species=="harengus")
sole<-subset(CSEENUMRAW, CSEENUMRAW$Species=="vulgaris")
plai<-subset(CSEENUMRAW, CSEENUMRAW$Species=="platessa")

cod<-as.data.frame(cod)
had<-as.data.frame(had)
her<-as.data.frame(her)
sole<-as.data.frame(sole)
plai<-as.data.frame(plai)

SPNRAW <- rbind(cod,had,her,sole,plai)
write.csv(SPNRAW, file="SPPNUMRAW.csv")
summary(SPNRAW)

#gives 26 morhua, 24 aegle, 22 hareng, 14 vulg, 10 plat /2 (want young and old)
#load new dataframe containing only those species

rm(list=ls(all=T))
SPNLIM<-read.csv("SPPNUMRAW.csv", header=T)
summary(SPNLIM)
DNUM<-SPNLIM$DNUM
PDNUM<-SPNLIM$PDNUM
Order<-SPNLIM$Order
name<-SPNLIM$name
Species<-SPNLIM$Species
MGMT<-SPNLIM$MGMT
Age<-SPNLIM$Age

summary(lm(DNUM~Species))
summary.aov(lm(DNUM~Species))
#Graph N7
boxplot(DNUM~Species, main="Difference in Number for Selected Species", 
        col=(c("lightpink", "salmon", "lightpink","palevioletred4", "palevioletred4")))

summary(lm(DBM~Age*Species))
summary(lmer(DBM~Age+(1+Age|Species)))
summary.aov(lm(DBM~Age*Species))
