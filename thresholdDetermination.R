# quantileDetermination.R
# K. Beck
# Jan. 15, 2013

# Load data with exprCountsLoad.R
options(scipen= 10)

# Objective: Use quantiles to determine the empirical threshold for calculating the dilution adjustment 
# which is fed to dilutioneffect.pl

# # # # # # 
# 
# FUNCTIONS
#
# # # # # # 

determineThreshold = 
    function (data, column = 2, q = 0.9995) {
        # Store the column of interest in a temporary variable
        tmp = as.numeric(data[,column])
        
        # Determine the value along the cumulative distribution function at the quantile specified
        qq = quantile(tmp, probs = q, names = FALSE)
        totalExpr = sum(tmp) # get total expression
        threshold = qq/totalExpr # calculate the threshold relating to the quantile
        
        cat("There are", length(tmp[tmp > qq]), "high abundance genes.\n")
        highGenes = tmp[tmp > qq]
        cat("The smallest of which is", min(highGenes), "\n")
        
        # Estimate the dilution adjustment
        aboveQuantile = sum(tmp[tmp > qq]) 
        adjFactor = 1/(1-(aboveQuantile/totalExpr))
        cat("The estimated diution adjustment is", adjFactor, "\n")
        
        cat("The quantile determined threshold is", threshold, "\n")
        
        return(threshold)
    }


# # # # # # # #

## Determine thresholds for both cow and human data at each lactation stage
# Cow
threshHarhayPrepub = sapply(2:7, function(x) determineThreshold(prepuberty, column = x, q= 0.9995))
threshHarhay = sapply(2:7, function(x) determineThreshold(lactation, column = x, q= 0.9995))

# Human
threshNRcolostrum = sapply(2:3, function(x) determineThreshold(NRcolostrum, column = x, q= 0.9995)) # dilution adjustment ~1.21
threshNRtransitional = sapply(2:5, function(x) determineThreshold(NRtransitional, column = x, q= 0.9995)) # dilution adjustment ~1.21
threshNRmature = sapply(2:7, function(x) determineThreshold(NRmature, column = x, q= 0.9995))

# These thresholds can be used as input for dilution_effect.pl for each replicate



# # # # # # # #

# TODO maybe report the high abundance genes for comparison between human lactation
intersect(highGenesHs$GeneID, highGenesHsC$GeneID) # 11 shared high abundance genes between colostrum and mature


########## Quantile and other statistics research ##########

# Explore quantile values
# returns the x value at this portion of the cumultive distribution curve
# quantile(lactation$mean, probs = c(0.05, 0.95, 0.999, .999908, 1), names=TRUE) # cow
# quantile(NRmature$mean, probs = c(0.05, 0.95, 0.999, .999908, 1), names=TRUE)  # human



## How similar are the distributions of non-lactating and lactating?
## Mann-Whitney-Wilcoxon Test
# Null hypothesis: non-lactating and lactating are identical populations
# if p < 0.05 we can reject the null hypothesis
# wilcox.test(x=NRcolostrum$mean,    y=NRmature$mean) # returns p < 2.2e-16
# wilcox.test(x=NRtransitional$mean, y=NRmature$mean) # returns p < 2.2e-16
# wilcox.test(x=NRcolostrum$mean,    y=NRtransitional$mean) # p-value = 0.01624....ask Danielle about this case???
# therefore non-lac and lac are nonidentical populations

# wilcox.test(x=prepuberty$mean,     y=lactation$mean) # returns p < 2.2e-16
# therefore colostrum and mature are nonidentical populations
# ref: http://www.r-tutor.com/elementary-statistics/non-parametric-methods/mann-whitney-wilcoxon-test

# This needs m and n values, but how do you determine those?
# plot(dwilcox(lactation$mean))

## Kolmogorov-Smirnov Tests
# ks.test(x=prepuberty$mean, y=lactation$mean) # p < 2.2e-16
# ks.test(x=NRcolostrum$mean, y=NRmature$mean) # p < 2.2e-16


# Description of statistical modeling
# d = probability density values (pdf = probability distribution function)
# returns height of the distribution curve at specified (input) x value(s)
# in some functions you can specify the mean and sd to tailor it to your curve
# otherwise it treats x as the z score
# p = cumulative probabilities (cdf = cumulative distribution function)
# returns the area below the curve up until the value specified
# i.e. pnorm(1) returns the cumulative probabilities up until 1 in the curve specified
# defaults to mean = 0 and sd = 1
# to find the area above that point pnorm(1, lower.tail=FALSE)
# q = quantile values
# gives you quantiles or critical values
# qnorm(0.95) # p = 0.05, one-tailed (upper) will return 1.644854
# qnorm(0.95, mean = 5, sd =1) # can shift the curve to fit your data like so
# r = random numbers from the distribution
# ref: http://ww2.coastal.edu/kingw/statistics/R-tutorials/prob.html
