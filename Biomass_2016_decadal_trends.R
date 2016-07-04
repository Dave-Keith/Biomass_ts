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
      #db[[name]]$total_ratio   <- db[[name]]$total/max(db[[name]]$total,na.rm=T)
      #db[[name]]$total_stan   <- scale(log(db[[name]]$total))
      #db[[name]]$total_offset   <- max(db[[name]]$total,na.rm=T)
      
      db[[name]]$old   <- rowSums(subset(db[[name]],select = c(which(names(db[[name]]) %in% var.names.old))),na.rm=T)
      # Replace 0's with half minimum value in time series just so log works 
      if(any(db[[name]]$old == 0))db[[name]]$old[db[[name]]$old == 0] <- 0.5*min(db[[name]]$old[db[[name]]$old > 0])
      #db[[name]]$old_ratio   <- db[[name]]$old/max(db[[name]]$old,na.rm=T)
      #db[[name]]$old_stan   <- scale(log(db[[name]]$old))
      #db[[name]]$old_offset   <- max(db[[name]]$old,na.rm=T)
      
      db[[name]]$young <- rowSums(subset(db[[name]],select = c(which(names(db[[name]]) %in% var.names.young))),na.rm=T)
      # Replace 0's with half minimum value in time series just so log works 
      if(any(db[[name]]$young == 0))db[[name]]$young[db[[name]]$young == 0] <- 0.5*min(db[[name]]$young[db[[name]]$young > 0])
      #db[[name]]$young_ratio   <- db[[name]]$young/max(db[[name]]$young,na.rm=T)
      #db[[name]]$young_offset   <- max(db[[name]]$young,na.rm=T)
      #db[[name]]$young_stan   <- scale(log(db[[name]]$young))
      # Now we can remove the individual age data as that's not especially necessary to keep
      db[[name]] <- subset(db[[name]],select = -c(which(names(db[[name]]) %in% var.names)))
    } # end if(length(cols) > 0 && length(rows) > 0)
  } # end if(is.element(name,c("ALPLAICBSAIm","ARFLOUNDBSAIm","GHALBSAIm","NRSOLEEBSAIm","YSOLEBSAIm"))==F)
} # end for(i in 1:num.stocks)  
# Make sure the columns all the same length
range(lapply(db,ncol))

db <- do.call("rbind",db)

#############  End Section 1 Data Processing#############  End Section 1 Data Processing#############  End Section 1 Data Processing

##############  Section 2, the decadal models....
# So the goal here is to get the change in biomass decade by decade for each stock and look to see what's been happening
# The results can be used to get the trends over the decades by stock and overall...
names(db)
stocks <- unique(db$Stock.ID)
res.lst <- NULL
decades <- seq(1940,2010,by=10)

#tapply(db$Year,db$Stock.ID,min)
# Run this for each populations
pdf(paste("D:/Github/Current_papers/Biomass_ts/Figures/decadal_fits.pdf",sep=""),onefile=T)
par(mfrow=c(1,3),mar=c(2,2,4,1))
for(i in 1:length(stocks))
{
  # Set up the results list...
  res.lst[[i]] <- data.frame(years = decades[1:7], Stock.ID = dat$Stock.ID[1:7],Order = dat$Order[1:7],
                             Genus = dat$Genus[1:7], Species = dat$Species[1:7],Management = dat$Management[1:7],
                             r.tot = rep(NA,7),r.young = rep(NA,7),r.old = rep(NA,7),
                             r.se.tot = rep(NA,7),r.se.young = rep(NA,7),r.se.old=rep(NA,7),num.years = rep(NA,7),
                             pred.total.first=NA,pred.total.last=NA,pred.old.first=NA,pred.old.last=NA,
                             pred.young.first=NA,pred.young.last=NA,
                             y.total.first=NA,y.total.last=NA,y.old.first=NA,y.old.last=NA,
                             y.young.first=NA,y.young.last=NA)
  # Get the data.
  dat <- db[db$Stock.ID == stocks[i],]
  # Now loop through each decade
  for(j in 1:(length(decades)-1))
  {
    # Get the current decade
    cur.dec <- seq(decades[j],(decades[j+1]-1),1)
    # If there are any data for the j'th decade do the calculations
    if(any(dat$Year %in% cur.dec)==T)
    {
      # Get the number of years of data for the decade
      res.lst[[i]]$num.years[j] <- length(which(dat$Year %in% cur.dec))
      # This gets the data for the end of the previous decade when there is just one data point
      # we use this as the starting point of the model like we have a break-point regression
      # This way we are now using all our data for the stocks that start counts in the year before the start of the decade.
      if(res.lst[[i]]$num.years[j] == 1)
      {
        ender.total <- log(dat$total[which(dat$Year %in% cur.dec)])
        ender.old   <- log(dat$young[which(dat$Year %in% cur.dec)])  
        ender.young <- log(dat$old[which(dat$Year %in% cur.dec)])
      }# end if(res.lst[[i]]$num.years[j] == 1)
      # Run if we have more than 1 year of data.
      if(res.lst[[i]]$num.years[j] > 1)
      {
        # Get the x's and y's
        x <- 1:res.lst[[i]]$num.years[j]
        y.total <- log(dat$total[which(dat$Year %in% cur.dec)])
        y.young <- log(dat$young[which(dat$Year %in% cur.dec)])
        y.old <- log(dat$old[which(dat$Year %in% cur.dec)])
        # If we have data for j = 1 (i.e. the 1940's) then we need to run this...
        if(j == 1)
        {
          mod.total <-lm(y.total~x)
          mod.old <-lm(y.old~x)
          mod.young <-lm(y.young~x)
          # Get the predictions
          pred.total <- predict(mod.total)
          pred.old <- predict(mod.old)
          pred.young <- predict(mod.young)
          
          ender.total <- y.total[length(y.total)]
          ender.old <-   y.old[length(y.old)]
          ender.young <- y.young[length(y.young)]
          
          res.lst[[i]]$r.tot[j] <- coef(mod.total)[2]
          res.lst[[i]]$r.old[j] <- coef(mod.old)[2]
          res.lst[[i]]$r.young[j] <- coef(mod.young)[2]
          # For our NCEAS work we used 2*se, here I'm just using the actual SE...
          res.lst[[i]]$r.se.tot[j]  <- summary(mod.total)$coefficients[2,2]
          res.lst[[i]]$r.se.old[j]  <- summary(mod.old)$coefficients[2,2]
          res.lst[[i]]$r.se.young[j]  <- summary(mod.young)$coefficients[2,2]
          res.lst[[i]]$pred.total.last[j]   <- pred.total[length(pred.total)]
          res.lst[[i]]$pred.total.first[j]  <- pred.total[1]
          res.lst[[i]]$pred.old.last[j]   <- pred.old[length(pred.total)]
          res.lst[[i]]$pred.old.first[j]  <- pred.old[1]
          res.lst[[i]]$pred.young.last[j]   <- pred.young[length(pred.total)]
          res.lst[[i]]$pred.young.first[j]  <- pred.young[1]
          # Get the first and last datapoints
          res.lst[[i]]$y.total.last[j]   <- y.total[length(y.total)]
          res.lst[[i]]$y.total.first[j]  <- y.total[1]
          res.lst[[i]]$y.old.last[j]   <- y.old[length(y.old)]
          res.lst[[i]]$y.old.first[j]  <- y.old[1]
          res.lst[[i]]$y.young.last[j]   <- y.young[length(y.young)]
          res.lst[[i]]$y.young.first[j]  <- y.young[1]
          
          # Make the plots
          plot(exp(y.total)~x,log="y",
               main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
               cex.main=0.7,bty="L")
          lines(exp(pred.total)~x)
          plot(exp(y.old)~x,log="y",
               main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
               cex.main=0.7,bty="L")
          lines(exp(pred.old)~x)
          plot(exp(y.young)~x,log="y",
               main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
               cex.main=0.7,bty="L")
          lines(exp(pred.young)~x)
        }# end if(j == 1)
        
        #If we didn't have data to get a population growth rate for the previous decade we run this bit of code
        if(j > 1)
        { 
        if(is.na(res.lst[[i]]$r.tot[j-1]) == T) 
        {
          # This will only run if we had 1 year of data in the previous decade, we start the regression using that value
          if(is.na(res.lst[[i]]$num.years[j-1])==F)
          {
            mod.total <- lm(I(y.total-ender.total) ~  x +0)
            mod.old <- lm(I(y.old-ender.old) ~  x +0)
            mod.young <- lm(I(y.young-ender.young) ~  x +0)
            
            pred.total <- predict(mod.total)+ender.total
            pred.old <- predict(mod.old)+ender.old
            pred.young <- predict(mod.young)+ender.young
            
            res.lst[[i]]$r.tot[j] <- coef(mod.total)[1]
            res.lst[[i]]$r.old[j] <- coef(mod.old)[1]
            res.lst[[i]]$r.young[j] <- coef(mod.young)[1]
            # For our NCEAS work we used 2*se, here I'm just using the actual SE...
            res.lst[[i]]$r.se.tot[j]  <- summary(mod.total)$coefficients[1,2]
            res.lst[[i]]$r.se.old[j]  <- summary(mod.old)$coefficients[1,2]
            res.lst[[i]]$r.se.young[j]  <- summary(mod.young)$coefficients[1,2]
            res.lst[[i]]$pred.total.last[j]   <- pred.total[length(pred.total)]
            res.lst[[i]]$pred.total.first[j]  <- pred.total[1]
            res.lst[[i]]$pred.old.last[j]   <- pred.old[length(pred.old)]
            res.lst[[i]]$pred.old.first[j]  <- pred.old[1]
            res.lst[[i]]$pred.young.last[j]   <- pred.young[length(pred.young)]
            res.lst[[i]]$pred.young.first[j]  <- pred.young[1]
            # Get the first and last datapoints
            res.lst[[i]]$y.total.last[j]   <- y.total[length(y.total)]
            res.lst[[i]]$y.total.first[j]  <- y.total[1]
            res.lst[[i]]$y.old.last[j]   <- y.old[length(y.old)]
            res.lst[[i]]$y.old.first[j]  <- y.old[1]
            res.lst[[i]]$y.young.last[j]   <- y.young[length(y.young)]
            res.lst[[i]]$y.young.first[j]  <- y.young[1]
            
            # Make the plots
            plot(exp(y.total)~x,log="y",
                 main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
                 cex.main=0.7,bty="L")
            lines(exp(pred.total)~x)
            plot(exp(y.old)~x,log="y",
                 main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
                 cex.main=0.7,bty="L")
            lines(exp(pred.old)~x)
            plot(exp(y.young)~x,log="y",
                 main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
                 cex.main=0.7,bty="L")
            lines(exp(pred.young)~x)
          } # end if(res.lst[[i]]$num.years[j-1] == 1)
          
          # This runs if we didn't have any data in the previous decade
          if(is.na(res.lst[[i]]$num.years[j-1])==T)
          {
            mod.total <-lm(y.total~x)
            mod.old <-lm(y.old~x)
            mod.young <-lm(y.young~x)
            # Get the predictions
            pred.total <- predict(mod.total)
            pred.old <- predict(mod.old)
            pred.young <- predict(mod.young)
            
            ender.total <- pred.total[length(pred.total)]
            ender.old <- pred.old[length(pred.old)]
            ender.young <- pred.young[length(pred.young)]
            
            res.lst[[i]]$r.tot[j] <- coef(mod.total)[2]
            res.lst[[i]]$r.old[j] <- coef(mod.old)[2]
            res.lst[[i]]$r.young[j] <- coef(mod.young)[2]
            # For our NCEAS work we used 2*se, here I'm just using the actual SE...
            res.lst[[i]]$r.se.tot[j]  <- summary(mod.total)$coefficients[2,2]
            res.lst[[i]]$r.se.old[j]  <- summary(mod.old)$coefficients[2,2]
            res.lst[[i]]$r.se.young[j]  <- summary(mod.young)$coefficients[2,2]
            res.lst[[i]]$pred.total.last[j]   <- pred.total[length(pred.total)]
            res.lst[[i]]$pred.total.first[j]  <- pred.total[1]
            res.lst[[i]]$pred.old.last[j]   <- pred.old[length(pred.old)]
            res.lst[[i]]$pred.old.first[j]  <- pred.old[1]
            res.lst[[i]]$pred.young.last[j]   <- pred.young[length(pred.young)]
            res.lst[[i]]$pred.young.first[j]  <- pred.young[1]
            # Get the first and last datapoints
            res.lst[[i]]$y.total.last[j]   <- y.total[length(y.total)]
            res.lst[[i]]$y.total.first[j]  <- y.total[1]
            res.lst[[i]]$y.old.last[j]   <- y.old[length(y.old)]
            res.lst[[i]]$y.old.first[j]  <- y.old[1]
            res.lst[[i]]$y.young.last[j]   <- y.young[length(y.young)]
            res.lst[[i]]$y.young.first[j]  <- y.young[1]
            
            # Make the plots
            plot(exp(y.total)~x,log="y",
                 main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
                 cex.main=0.7,bty="L")
            lines(exp(pred.total)~x)
            plot(exp(y.old)~x,log="y",
                 main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
                 cex.main=0.7,bty="L")
            lines(exp(pred.old)~x)
            plot(exp(y.young)~x,log="y",
                 main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
                 cex.main=0.7,bty="L")
            lines(exp(pred.young)~x)
          } # end if(res.lst[[i]]$num.years[j-1] != 1)
          
        } # end if(is.na(res.lst[[i]]$r.tot[j-1]) ==T) 
        
        # Now if we had data from before we do the breakpoint regression...
        if(is.na(res.lst[[i]]$r.tot[j-1]) ==F) 
        {
          # Get the start point for the regression.
          ender.total <- as.numeric(res.lst[[i]]$pred.total.last[j-1])
          ender.old   <- as.numeric(res.lst[[i]]$pred.old.last[j-1])
          ender.young <- as.numeric(res.lst[[i]]$pred.young.last[j-1])
        #} # end } # end if(is.na(res.lst[[i]]$r.tot[j-1]) ==F) 

        # Now this essentially uses the last predicted point from previous regression as the intercept for
        # the subsequent regression, kinda cool!  What I'm actually doing is a regression with a 0 intercept
        # and subtracting off the "intercept" using the last point from the previous regression.
        # Get the x's and y's
        mod.total <- lm(I(y.total-ender.total) ~  x +0)
        mod.old <- lm(I(y.old-ender.old) ~  x +0)
        mod.young <- lm(I(y.young-ender.young) ~  x +0)
        
        ### ENDED LAST NIGHT SOMETHING IS WRONG HERE!!!!
        pred.total <- predict(mod.total)+ender.total
        pred.old <- predict(mod.old)+ender.old
        pred.young <- predict(mod.young)+ender.young
        
        res.lst[[i]]$r.tot[j] <- coef(mod.total)[1]
        res.lst[[i]]$r.old[j] <- coef(mod.old)[1]
        res.lst[[i]]$r.young[j] <- coef(mod.young)[1]
        # For our NCEAS work we used 2*se, here I'm just using the actual SE...
        res.lst[[i]]$r.se.tot[j]  <- summary(mod.total)$coefficients[1,2]
        res.lst[[i]]$r.se.old[j]  <- summary(mod.old)$coefficients[1,2]
        res.lst[[i]]$r.se.young[j]  <- summary(mod.young)$coefficients[1,2]
        res.lst[[i]]$pred.total.last[j]   <- pred.total[length(pred.total)]
        res.lst[[i]]$pred.total.first[j]  <- pred.total[1]
        res.lst[[i]]$pred.old.last[j]   <- pred.old[length(pred.old)]
        res.lst[[i]]$pred.old.first[j]  <- pred.old[1]
        res.lst[[i]]$pred.young.last[j]   <- pred.young[length(pred.young)]
        res.lst[[i]]$pred.young.first[j]  <- pred.young[1]
        # Get the first and last datapoints
        res.lst[[i]]$y.total.last[j]   <- y.total[length(y.total)]
        res.lst[[i]]$y.total.first[j]  <- y.total[1]
        res.lst[[i]]$y.old.last[j]   <- y.old[length(y.old)]
        res.lst[[i]]$y.old.first[j]  <- y.old[1]
        res.lst[[i]]$y.young.last[j]   <- y.young[length(y.young)]
        res.lst[[i]]$y.young.first[j]  <- y.young[1]
        
        # Make the plots
        plot(exp(y.total)~x,log="y",
             main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
             cex.main=0.7,bty="L")
        lines(exp(pred.total)~x)
        plot(exp(y.old)~x,log="y",
             main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
             cex.main=0.7,bty="L")
        lines(exp(pred.old)~x)
        plot(exp(y.young)~x,log="y",
             main=paste(dat$Stock.ID[1],"-",dat$Species[1],"-",min(cur.dec),"-",max(cur.dec),sep=""),
             cex.main=0.7,bty="L")
        lines(exp(pred.young)~x)
        
        
      } # end if(is.na(res.lst[[i]]$r.tot[j-1]) ==F) 
      } # end if(res.lst[[i]]$num.years[j] > 1)
    } # end if(any(dat$Year %in% cur.dec)==T)
    } # end if(j > 1)  

    
  } # end for(j in 1:(length(decades)-1))
} # end for(i in 1:length(stocks))
dev.off()
res.db <- do.call("rbind",res.lst)

# Re-arrange the data so we have the data set up to use for a model to compare different r values.
names(res.db)
tst1 <- res.db[,c(-grep("young",names(res.db)),-grep("old",names(res.db)))]
tst1 <- cbind(tst1,rep("total",nrow(tst1)))

tst2 <- res.db[,c(-grep("young",names(res.db)),-grep("total",names(res.db)),-grep("tot",names(res.db)))]
tst2 <- cbind(tst2,rep("old",nrow(tst1)))

tst3 <- res.db[,c(-grep("old",names(res.db)),-grep("total",names(res.db)),-grep("tot",names(res.db)))]
tst3 <- cbind(tst3,rep("young",nrow(tst1)))

colnames(tst1) <- c("years","Stock.ID","Order","Genus","species","Management","r","se","num.years","pred.first",
                    "pred.last","y.first","y.last","lambda","age") 
colnames(tst2) <- colnames(tst1)
colnames(tst3) <- colnames(tst1)

mod.db <- rbind(tst1,tst2,tst3)

# So now that I have the data let's start looking at it...
par(mfrow=c(3,1))
plot(res.db$r.tot)
plot(res.db$r.young)
plot(res.db$r.old)

par(mfrow=c(3,1))
plot(res.db$r.tot~res.db$years)
plot(res.db$r.old~res.db$years)
plot(res.db$r.young~res.db$years)

par(mfrow=c(3,1))
plot(res.db$r.tot~res.db$Management)
plot(res.db$r.old~res.db$Management)
plot(res.db$r.young~res.db$Management)

par(mfrow=c(3,1))
plot(res.db$r.tot~res.db$Order)
plot(res.db$r.old~res.db$Order)
plot(res.db$r.young~res.db$Order)


median(res.db$r.tot,na.rm=T)
median(res.db$r.young,na.rm=T)
median(res.db$r.old,na.rm=T)

aggregate(r.tot~years+Order,res.db,FUN=mean)
aggregate(r.young~years+Order,res.db,FUN=mean)
aggregate(r.old~years+Order,res.db,FUN=mean)

aggregate(r.tot~years,res.db,FUN=median)
aggregate(r.young~years,res.db,FUN=median)
aggregate(r.old~years,res.db,FUN=median)

aggregate(r.tot~years+Management,res.db,FUN=mean)
aggregate(r.young~years+Management,res.db,FUN=mean)
aggregate(r.old~years+Management,res.db,FUN=mean)

apply(res.db[res.db$Management == "DFO",c(25:26)],2,function(x) prod(x,na.rm=T))

res.db$lambda.tot <-  exp(res.db$r.tot)
res.db$lambda.young <-  exp(res.db$r.young)
res.db$lambda.old <-  exp(res.db$r.old)
median(aggregate(lambda.tot~Stock.ID,res.db,FUN=prod)$lambda.tot,na.rm=T)
median(aggregate(lambda.young~Stock.ID,res.db,FUN=prod)$lambda.young,na.rm=T)
median(aggregate(lambda.old~Stock.ID,res.db,FUN=prod)$lambda.old,na.rm=T)


mod.1 <- lm(r~age,mod.db)
summary(mod.1)
AIC(mod.1)

mod.2 <- lm(r~age*years,mod.db)
summary(mod.2)
AIC(mod.1,mod.2)

mod.3 <- lm(r~age*years*Order,mod.db)
summary(mod.3)
AIC(mod.1,mod.2,mod.3)

mod.1a <- lm(r.tot~as.factor(years)*Order,res.db)
summary(mod.1a)
