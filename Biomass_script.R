#########################################################################################################
###
###             Here I'm processing the data to get it how I want it for the analyses and figures...


rm(list=ls(all=T))
#ASD <- read.table("/home/keithdm/Documents/PhD/R_programs/Age_structured_database/ASD_Final.csv",
#                  header=T,sep =",")#,na.strings="-9999")
load("C:/Users/Emilie/Dropbox/ASD/Analyses/Processed_DB.RData")
#load("D:/Dropbox/ASD/Analyses/Processed_DB.RData")
# This ran the below loop to tidy up the database into something usable!
setwd("D:/Dropbox/ASD/Analyses")
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
    # Now pick the Biomass for this stock
    st <- BM.calc[stocks == name,] 
#BM.calc is the Num.0:Num.31 for all the stocks, so including those for which there is no Num available
#BM.calc[stocks == name,] is te biomass in numbers for stock i=5
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
    BM_geo_LQS <- geo.mean(BM_LQ_start_tot,years)
    BM_geo_UQS <- geo.mean(BM_UQ_start_tot,years)
    BM_geo_LQN <- geo.mean(BM_LQ_now_tot,years)
    BM_geo_UQN <- geo.mean(BM_UQ_now_tot,years)
    # Now quickly calculate the total decline, negative if decline, positive is growth...
    BM_dec_LQ <- (BM_geo_LQN - BM_geo_LQS) 
    BM_dec_UQ <- (BM_geo_UQN - BM_geo_UQS) 
    # Can see both decline, about 7400 tonnes less in the Youngest age classes, and 8900 fewer tonnes in the oldest now
 
    # Now quickly calculate the percent decline, negative if decline, positive is growth...
    BM_per_dec_LQ <- (BM_geo_LQN - BM_geo_LQS) / BM_geo_LQS
    BM_per_dec_UQ <- (BM_geo_UQN - BM_geo_UQS) / BM_geo_UQS
    # Can see both decline, but BM of UQ declines by 17% more than LQ, cool beans!
    
#    } # End the if loop
    
 #      } # End the for loop
#save(ASD_FR,file="F_ssb.RData")

