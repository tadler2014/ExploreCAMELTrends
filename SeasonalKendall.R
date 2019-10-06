### trendAnalyis_VANR_v2.R 

# Orginally Created by Kristen Underwood
# 12 February 2017

# Reconfigured by Thomas Adler
# November 25th 2018

# Background Resources:
# FROM "DETECTING TRENDS..." ch 16 stored in c:R ACRWC\Trend Analysis
# Helsel & Hirsch, 1992
# EPA Trend Mon Tech Note

# Part 0 - Housekeeping (install pkgs, read in data, define model parameters)
# Part 1 - Pre-Filtering for Gaps and Lack of Data (>5 years and <1/3rd gaps)
# Part ? - Examine Q record for stationarity (MannKendall, acf, slope)
# Part 2 - Simple Trend Test: Not Adjusted for X
#      2a - Parametric: Regress C on Time
#      2b - Nonparametric: Mann Kendall C on Time
# Part 3 - Trend Adjusted for X
#      3a - Parametric:  C Regressed on X (Flow) and Time
#      3b - Mixed: Regress Conc on Discharge
#             First: Parametric: regress C (LogConc) on X (Flow) to get residuals (R)
#             Then: Nonparametric: Mann Kendall trend test of R on Time
#             Alternatively: t-Test on Before and After gap in record of flow-adjusted C (residuals)
#      3c - Nonparametric: Flow Adjusted Conc (R) on Time 

#------------------- Housekeeping ----------------------------------------------
# Set working directory
# setwd("C:/R ACRWC/Trend Analysis")

# Install Packages
#install.packages("trend")  #for running Mann-Kendall (requires input vector as time series, and I think of regularly spaced intervals)
install.packages("plyr") #for running counts exceeding thresholds by groupings
install.packages("Kendall") #for running Mann-Kendall
install.packages("zyp") #for estimating slope of MK trend
install.packages('xts')
install.packages("ggplot2")

library(ggplot2)
library(zoo)
#library(trend)
library(plyr)
library(Kendall)
library(zyp)
library(xts)

#------------------- Data Read-In and Define Inputs ----------------------------------------------
# Read in Data
myData <- read.csv("Input Data/CAMELSchem.csv")

# Choose variable(s) of interest for trend analysis
v1 <- as.numeric(as.character(myData$doc))
vname = "DOC [mg/L]" # Define variable name

# Choose weighting/adjusting variable(s) of interest for trend analysis (i.e. flow)
wv1 = as.numeric(as.character(myData$q_inst))
wvname = "Instanteous Flow [Q]" # Define weighting variable name

# Define Date Column
Date <- as.Date(myData$sample_start_dt)

#------------------- Create Data Frame ----------------------------------------------
# Create dataframe for selected variable
myData1 <- data.frame("gauge_id" = myData$gauge_id,"date" = Date, "v1" = v1,"wv1" = wv1)
# Remove missing data
myData1 <- subset(myData1, !is.na(myData1$v1))
# Assign Year and Month to each entry
myData1$yr <- strftime(myData1$date, "%Y")  #append column containing year
myData1$mo <- strftime(myData1$date, "%m")  #append column containing month
# Identify sites by Gauge ID and total number of sites by Gauge ID
Gauges=as.numeric(unique(myData1$gauge_id))
numGauges=length(Gauges)
# Order data frame by date 
myData1 <- myData1[order(myData1$date),] 
# Remove temp files
rm(Date, v1, wv1)
#------------------- Filtering for Gaps and Lack of Data ----------------------------------------------
# Filter 1: At least 5 years of monthly data
preFilter1<- data.frame("Gauge_ID"=as.numeric(),"Years_of_data"=as.numeric())
for (i in 1:numGauges){
  #Select Gauge/Site to analyze
  g <- Gauges[i] 
  #Create dataset for selected gauge
  C <- subset(myData1, myData1$gauge_id==g, select=c(yr, mo, v1))# sample date
  #Determine amount of years in data set
  nCyrs=length((unique(C$yr)))
  #Store results of loop in prefilter data table
  preFilter1[i,] <- c(g, nCyrs)
}  
Filter1 <- subset(preFilter1, preFilter1$Years_of_data >= 5)
rm(C, nCyrs)


# Filter 2: No more than 1/3rd gap in data
preFilter2 <- data.frame("Gauge_ID"=as.numeric(),"BiggestGap"=as.numeric(),"OneThirdTotalDays"=as.numeric())
for (i in 1:nrow(Filter1)){
  #Select Gauge/Site to analyze
  g <- Filter1$Gauge_ID[i] 
  #Create dataset for selected gauge
  D <- subset(myData1, myData1$gauge_id==g, select=c(date,v1,yr,mo))
  D_diff <- as.numeric(abs(difftime(tail(D$date, -1),head(D$date, -1),units = "days"))) #Gap between each date
  D_diff_tot <- as.numeric(abs(difftime(D$date[1],D$date[nrow(D)],units = "days")))
  #Determine the largest gap
  maxDdiff=max(D_diff) 
  #Determine 1/3 of total data set
  oneTtot= D_diff_tot / 3 
  #If the gap is larger than 1/3rd the data
  while (maxDdiff > oneTtot){ 
    #Find the location of the gap
    l <- which.max(D_diff) 
    #Find the length of the data set below gap
    llow=length((unique(D$yr[1:l]))) 
    #Find the length of the data set above gap
    lhigh=length((unique(D$yr[(l+1):length(D$yr)]))) 
    #If both sides of the gap are less than 5 years, move on
    if (llow<5 & lhigh < 5){ 
    break;
    }
    #If the amount of data above the gap is greater than the amount below
    if (lhigh >= llow){ 
      #Only hold on to data above gap
      myData1 <- myData1[!(myData1$gauge_id == g & myData1$date < D$date[(l+1)]),]
      #Take subset of the data above gap
      D <- D[(l+1):length(D$date),] 
    }
    #If the amount of data below the gap is greater than the amount above
    if (llow > lhigh){ 
      #Only hold on to data below gap
      myData1 <- myData1[!(myData1$gauge_id == g & myData1$date > D$date[(l)]),]
      #Take subset of the data above gap and rerun gap assessment
      D <- D[1:l,] 
    }
    #Rerun gap assessment 
    D_diff <- as.numeric(abs(difftime(tail(D$date, -1),head(D$date, -1),units = "days"))) #Gap between each date
    D_diff_tot <- as.numeric(abs(difftime(D$date[1],D$date[nrow(D)],units = "days")))
    maxDdiff=max(D_diff)
    oneTtot = D_diff_tot / 3
      
  }
  #Store results of loop in prefilter data table
  preFilter2[i,] <- c(g, maxDdiff, oneTtot)
}
Filter2 <- subset(preFilter2, preFilter2$OneThirdTotalDays >= preFilter2$BiggestGap)
rm(i, D, D_diff, D_diff_tot, g, l, lhigh, llow, maxDdiff, oneTtot)


#-----------------LOWESS application to flow & DOC ---------------------------------
# Define smoothing span (Suggested value = 2/3)
# Smoothing span value is any where from 0 (Jimi Hendrix) to 1 (Carlos Santana) 
mySpan = 2/3 
# Potential to come up with code that determines the correct smoothing span?

for (i in 1:nrow(Filter2)){
  # Identify gauge_ID
  g <- Filter2$Gauge_ID[i] 
  # Create dataset for selected gauge
  temp <- subset(myData1, myData1$gauge_id==g & myData1$wv1>0, select=c(gauge_id,date,v1,wv1,yr,mo))
  # Input dataset variables into flow-adjuster function
  myResiduals <- flow.adjuster(temp$v1,temp$wv1,mySpan)
  # Reorder main data frame by flow (aka the independent value of lowess regression)
  temp <- temp[order(log(temp$wv1)),]
  # Add residual to data frame (residual should already be ordered by independent value)
  temp$v1Residuals <- myResiduals
  # Reorder temporary data by date
  temp <- temp[order(temp$date),]
  # Bind data frames with each loop
  if (i==1){
    myData2 <- temp
  } else {
    myData2 <- rbind(myData2,temp)
  }
  rm(temp, myResiduals,g,i)
}
  
#-----------------Seasonal Tests (NEW) -----------------------------------

Summary_data <- data.frame(Gauge=numeric(),tau=numeric(),sl=numeric(),S=numeric(),D=numeric(),varS=numeric(),firstDate=factor(),lastDate=factor(),stringsAsFactors = FALSE)
for (i in 1:nrow(Filter2)){
  g <- Filter2$Gauge_ID[i] #Select Gauge/Site to analyze
  sk <- subset(myData2, myData2$gauge_id==g)# sample date
  #Aggregate variable on months and year and get median
  x=aggregate(v1Residuals ~ mo + yr , sk , median )
  x$Date <- as.yearmon(paste(x$yr, x$mo), "%Y %m")
  #Convert dataframe to time series
  z <- read.zoo(x[c("Date", "v1Residuals")], FUN = as.yearmon, format = "%b %Y")
  tt <- as.ts(z)
  #Plot the data
  autoplot(z,main=g) + scale_x_yearmon()
  #Calculate the Kendall stats
  res<-SeasonalMannKendall(tt) #Seasonal Kendall
  #Determine dates of first and last values
  NonNAindex <- which(!is.na(sk$v1Residuals))
  firstNonNA <- min(NonNAindex)
  lastNonNA <-max(NonNAindex)
  firstDate <- sk$date[firstNonNA]
  lastDate <- sk$date[lastNonNA]
  temp <- data.frame(Gauge=g,tau = res$tau,sl=res$sl,S=res$S,D=res$D,varS=res$varS,firstDate=firstDate,lastDate=lastDate)
  Summary_data <- rbind(temp, Summary_data)
}
  
  








#-----------------Seasonal Tests (NEW.2) -----------------------------------
library(EnvStats)
kendallSeasonalTrendTest(v1 ~ mo + yr,
                         data = sk)




#-----------------Seasonal Tests-----------------------------------
Summary_data <- data.frame(Site=character(),Gauge=numeric(),tau=numeric(),sl=numeric(),S=numeric(),D=numeric(),varS=numeric(),firstDate=factor(),lastDate=factor(),stringsAsFactors = FALSE)
for (i in 1:nrow(Filter2)){
#-----Define Model Parameters ----
g <- Filter2$Gauge_ID[i] #Select Gauge/Site to analyze
myFileNameRoot <- subset(A, A$gauge_id==g, select=gauge_name)
myFileNameRoot <- myFileNameRoot[1,]
myFileNameRoot <- as.character(myFileNameRoot)
myWshd <- as.character(g)

#-----------Create Watershed (Gauge) Dataset-------------
# Choose data
B_date <- subset(A, A$gauge_id==g, select=sample_start_dt)# sample date
B_date <- as.Date(B_date[,1])
B_flow <- subset(A, A$gauge_id==g, select=q_1) #sample flow
B_DOC <- subset(A, A$gauge_id==g, select=doc) #sample DOC
B_gauge <- subset(A, A$gauge_id==g, select=gauge_id) #sample gauge
B_site <- subset(A, A$gauge_id==g, select=gauge_name)#sample site
# Compile data
B <- data.frame(B_site, B_gauge, B_date, B_flow, B_DOC)
rm(B_DOC,B_flow,B_gauge,B_site)
#View data
#str(B)
#head(B)
#tail(B)
#class(B)

# Assign Year and Month to each entry
B$yr <- strftime(B$B_date, "%Y")  #append column containing year
B$mo <- strftime(B$B_date, "%m")  #append column containing month

# first make sure WQ values are numeric, not factor
B$doc <- as.numeric(as.character(B$doc))
B$q_1 <- as.numeric(as.character(B$q_1))

#Log transform the WQ data
B$logQ_cfs <- log10(B$q_1)
B$log_doc <- log10(B$doc)

#Determine dates of first and last values
NonNAindex <- which(!is.na(B$doc))
firstNonNA <- min(NonNAindex)
lastNonNA <-max(NonNAindex)
firstDate <- B$B_date[firstNonNA]
lastDate <- B$B_date[lastNonNA]
#--------------Pre-Filtering of the data-----------------------------
#Count number of days by year that flows are above percentile
#Thresh<- quantile(B$q_1, probs = c(0.9), na.rm =TRUE) #90th percentile

#--------------Pre-processing--------------------------------------
#Aggregate DOC on months and year and get median
x=aggregate( doc ~ mo + yr , B , median )
x$Date <- as.yearmon(paste(x$yr, x$mo), "%Y %m")
#Aggregate flow on months and year and get median
y=aggregate( q_1 ~ mo + yr , B , median )
y$Date <- as.yearmon(paste(y$yr, y$mo), "%Y %m")
#Convert dataframe to time series
z <- read.zoo(x[c("Date", "doc")], FUN = as.yearmon, format = "%b %Y")
tt <- as.ts(z)

#Plot the data
autoplot(z,main=myFileNameRoot) + scale_x_yearmon()


#----------------Seasonal Kendall and Thiel-Sen--------------------------------
res<-SeasonalMannKendall(tt) #Seasonal Kendall


model.T = mblm(v1 ~ ndate, 
               data=x)

df <- data.frame(Site=myFileNameRoot,Gauge=g,tau = res$tau,sl=res$sl,S=res$S,D=res$D,varS=res$varS,firstDate=firstDate,lastDate=lastDate)

Summary_data <- rbind(df, Summary_data)
rm(df)
}



#-----------------Export-----------------------------------------
# Write CSV in R
write.csv(Summary_data, file = "Kendalldata.csv")











#################################################################
#Everything past this point is archive code from Kristen
#################################################################

#---------3c: Nonparametric: Flow Adjusted Conc (R) on Time 
# set LOWESS smoothing span
mySpan <- 0.67

# decide on one smoothing span
plot(B$logQ_cfs, B$log_doc, pch=16, col='gray50', xlab = "Log [Daily Mean Discharge (cfs)]", 
     ylab = paste0("Log [ ", Constituent, " ",DataType," ", Units," ]"))
lines(lowess(B$log_doc ~ B$logQ_cfs, f=mySpan), lwd=2, col=4)
myLowessFit <- lowess(B$log_doc ~ B$logQ_cfs, f=mySpan)
legend("topright", col = 4, lty=1, lwd = 2, cex = 1.2,
       legend = paste0("f = ", mySpan))
str(myLowessFit)

# extract the LLOWESS residuals
myResiduals_Lowess<-lmDataDF$y[order(lmDataDF$x)]-myLowessFit$y
# plot the residuals on LogConc- LogQ plot
par(mfrow=c(1,1), mar=c(4, 5, 1, 1), mgp = c(2, 0.5, 0), tck = -0.02,
    omi=c(0.5, 0.5, 0.75, 0.5))
plot(lmDataDF$x, myResiduals_Lowess, col='red', xlab = "Log [Daily Mean Discharge (cfs)]", 
     ylab = paste0("Residuals of LOWESS from \nLog [ ", Constituent, " ",DataType," ", Units," ] - Log Q Plot"))
abline(h=0,col='gray50', lty = 2)


# Or Use Loess function 
#loess.model <- loess(lmDataDF$y ~ lmDataDF$x, span=0.9, degree=1)
#loess.model
loess.model2 <- loess(B$log_doc ~ B$logQ_cfs, span=mySpan, degree=1)
loess.model2

plot(B$logQ_cfs, B$log_doc, pch=16, col='gray50', xlab = "Log [Daily Mean Discharge (cfs)]", 
     ylab = paste0("Log [ "DOC" ",DataType," ", mg/L," ]") )
#hat1 <- predict(loess.model)
#lines(lmDataDF$x[order(lmDataDF$x)], hat1[order(lmDataDF$x)], col="red", lwd=2)
hat2 <- predict(loess.model2)
lines(lmDataDF$x[order(lmDataDF$x)], hat2[order(lmDataDF$x)], col="blue", lwd=2)
legend("topright", col = 4, lty=1, lwd = 2, cex = 1.2,
       legend = paste0("f = ", mySpan))

# extract residuals from Loess model
myResiduals_Loess<-residuals(loess.model2)
# plot the residuals on LogConc- LogQ plot
plot(lmDataDF$x, myResiduals_Loess, col='red', xlab = "Log [Daily Mean Discharge (cfs)]", 
     ylab = paste0("Residuals of LOESS from \nLog [ ", Constituent, " ",DataType," ", Units," ] - Log Q Plot"))
abline(h=0,col='gray50', lty=2)

# Now reorder the residuals in Time to examine for trends
# Take original data set ordered in Time, and reorder by increasing Flow
# Append Residuals to reordered data set.
# Reorder appended dataset to be once again ordered in time
tempReorderDF <- lmDataDF[with(lmDataDF, order(x)), ]
tempReorderDF$r <- myResiduals_Loess
ReorderDF <- tempReorderDF[with(tempReorderDF, order(x2)), ]

# Run nonparametric trend test - MannKendall
# of residuals on time
# two-sided homogeinity test
# H0: S = 0 (no trend)
# HA: S != 0 (monotonic trend)
mkResults_Loess_residuals <- MannKendall(ReorderDF$r)

# Now to estimate the slope using the Sen estimator in the zyp package 
sen.slope_Loess_resid <- zyp.sen(r ~ x2, data=ReorderDF) 

#extract parameters for plotting 
MKPValue <- paste("2-sided pValue = ", round( mkResults_Loess_residuals[[2]][1],digits = 4) )
tauValue <- paste("tau Value = ", round( mkResults_Loess_residuals$tau,digits = 4) )
MKScore <- paste("Score = ", round(mkResults_Loess_residuals$S,digits = 4))
MKslope <- paste("MK Slope = ", round(sen.slope_Loess_resid$coef[2],digits = 4) )

# plot Loess residuals vs time
par(mfrow=c(1,1), mar=c(4, 5, 2, 1), mgp = c(2, 0.5, 0), tck = -0.02,
    omi=c(0.5, 0.5, 0.75, 0.5))
plot(ReorderDF$x2, ReorderDF$r, typ='p', col = 'red', main = "", xlab = "Time (# days since 1/1/1970)",
     ylab = paste0("Loess Residuals of ", Constituent, " LogC-LogQ Plot"))
abline(h=0,lty = 2, lwd=1, col='gray50')
lines(lowess(ReorderDF$x2, ReorderDF$r), lwd=3, col=2)

mtext(paste("n =", N), side = 3, line = 0, adj = c(0,0), col = "blue")
mtext(MKslope, side = 3, line = 1, adj = c(1,1), col = "blue")
mtext(MKScore, side = 3, line = 2, adj = c(1,1), col = "blue")
mtext(tauValue, side = 3, line = 3, adj = c(1,1), col = "blue")
mtext(MKPValue, side = 3, line = 4, adj = c(1,1), col = "blue")
legend("topright", col = c(2), lty=c(1), lwd = c(2), cex = 1.1,
       legend = c("f=0.8"))

#examine residuals for autocorrelation 
#myacf <- acf(ReorderDF$r)  #calls plot.acf by default which includes default 95% conf int
myacf <- acf(ReorderDF$r, lag.max = 180)  

# Print to console for cross chk  
summary(mkResults_Loess_residuals)
mkResults_Loess_residuals
sen.slope_Loess_resid$coef 
confint.zyp(sen.slope_Loess_resid)


########
# Alternative: examine diff in means of residuals Before and After
# with t-Tests
Resid_Bef <- ReorderDF[ which(ReorderDF$Treat == "Before"), ]
n_x3 <- length(Resid_Bef$x)
x3 <- Resid_Bef$r

Resid_Aft <- ReorderDF[ which(ReorderDF$Treat == "After"), ]
n_x4 <- length(Resid_Aft$x)
x4 <- Resid_Aft$r



#-----Alt 3b-1: Run two-sample Frequentist t-test ---------------------

# determine if mean of Gp 1 (Before) is different from mean of Gp 2 (After)

#set Confidence level and calculate alpha level
confLevel <- 0.95
alpha <- 1 - confLevel

# -----------two-sided t-test
# H0: true mean of Gp1 is = true mean of Gp2
# H1: true mean of Gp1 is != true mean of Gp2

my2G_Ftt_2sided <- t.test(x3, x4, conf.level = confLevel, alternative = "two.sided")
my2G_Ftt_2sided 

# If pVal < alpha
# reject H0 (that mean Gp1 is = mean of Gp2)
# so true mean of Gp1 may be less than or greater/equal to mean of Gp2

# -----------one-sided t-test
# H0: true mean of Gp1 is <= true mean of Gp2
# H1: true mean of Gp1 is > true mean of Gp2

my2G_Ftt_1sided <- t.test(x3, x4, conf.level = confLevel, alternative = "greater")
my2G_Ftt_1sided 

# If pVal < alpha
# reject H0 (that mean Gp1 is <= mean of Gp2 )
# so true mean of Gp1 may be greater than mean of GP2


# ----------- Box Plot of Residuals----
plotHeader <- c("Lewis Creek, LCR14")
#plotHeader <- c("Little Otter, MDC1.2")
ylabel <- c("Residuals of \nLogC-LogQ Regression")
mainTitle <- Constituent

b <- boxplot(x3, x4, plot=0) # Do the boxplot but do not show it
# Now b$n holds the counts for each factor, we're going to write them in names

par(mar=c(7,8,4.1,6))  #set print margins of the plot
# bottom margin = 6.4 to leave room for x label with vertical hatch mark labels
# rest of margins are default margins
boxplot(x3, x4,  xlab = "", ylab = ylabel, 
        names=c("Before", "After"), range = 0, las=2, cex=2, pch="*", col=(c("lightgrey","lightgrey")), 
        main = mainTitle, ylim = c(-2, 2), boxwex = 0.4)
mtext(paste(b$n), side=3, line=0.1, at=1:length(b$n), cex = 0.8)
mtext(paste("n ="), side=3, line=0.1, at=0.5)
mtext(plotHeader, side=3, line=3, cex=1.0)

#------Alt 3b-2: Run Nonparametric Wilcoxon Rank sum test---------------------

# determine if mean of Gp 1 (Downstream) is different from mean of Gp 2 (Upstream)

#set Confidence level and calculate alpha level 
confLevel <- 0.95
alpha <- 1 - confLevel

# -----------two-sided t-test
# H0: true mean of Gp1 is = true mean of Gp2
# H1: true mean of Gp1 is != true mean of Gp2 

my2G_Wilcox_2sided <- wilcox.test(x3, x4, conf.level=confLevel, alternative = "two.sided")
my2G_Wilcox_2sided

# -----------one-sided t-test
# H0: true mean of Gp1 is <= true mean of Gp2
# H1: true mean of Gp1 is > true mean of Gp2

my2G_Wilcox_1sided <- wilcox.test(x3, x4, conf.level=confLevel, alternative = "greater")
my2G_Wilcox_1sided



