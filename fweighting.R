# Flow Weighting Function

# Created by: Thomas Adler
# Last Updated: September 14th, 2019

# Input:
# v1 = variable of interest
# wv1 = weighting variable
# mySpan = smoothing span value

# Output:
# Residuals of lowess regression

# IMPORTANT:
# The output residuals are ordered by the flow. If adding them to a data frame make
# sure the data frame is also ordered by flow first.

# Example:
# Create dataset for selected gauge
#   d <- subset(df, df$gauge_id==g & df$wv1>0, select=c(gauge_id,date,v1,wv1,yr,mo))
# Input dataset variables into function
#   myResiduals <- flow.adjuster(d$v1,d$wv1,mySpan)
# Reorder main data frame by flow (aka the independent value of lowess regression)
#   D <- D[order(D$wv1.log),]
# Add residual to data frame (residual should already be ordered by independent value)
#   D$v1Residuals <- myResiduals
# Reorder main data by date
#   D <- D[order(D$date),]

flow.adjuster <- function(v1, wv1, mySpan){
  # Step 1: Transform Data
  v1.log <- log10(v1)
  wv1.log <- log10(wv1)
  
  # Step 2: Determine lowess fit (variable ~ weighting variable)
  myLowessFit <- lowess(v1.log~wv1.log, f=mySpan)

  # Plot Regression
  #plot(D$wv1.log,D$v1.log, main = as.character(g),xlab="Log(Daily Flow (Q))",ylab = "Log(DOC [mg/L])")
  #lines(myLowessFit$x,myLowessFit$y, lwd=2, col=9)
  
  # Step 3:Extract the LLOWESS residuals
  myResiduals<-v1.log[order(wv1.log)]-myLowessFit$y
  return(myResiduals)
  # Plot residuals to make sure there is no trend left
  #plot(wv1.log[order(wv1.log)], myResiduals, col='red', xlab = "Log [Daily Mean Discharge (cfs)]", 
  #ylab = paste0("Residuals of LOWESS from Log DOC-Q Plot"),main = g)
}



