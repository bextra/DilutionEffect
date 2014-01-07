# distributionModeling.R
# K. Beck
# January 23, 2013

## Explorative plots
# Absolute change non-lactating to lactating
plot(density(cca$Change),
     main = "Distribution of bovine absolute change")
plot(density(hha$Change),
     main = "Distribution of human absolute change")
# Log 2 fold change non-lactating to lactating
plot(density(cc$Change))
plot(density(hh$Change))
hist(cc$Change,
     breaks= 1000,
     main = "Histogram of log2 fold change in cow samples")
hist(hh$Change,
     breaks= 2000,
     main = "Histogram of log2 fold change in human samples")

# # # # # # # # # # # # # # # # # # # # # #
#
#      Notes on Statistical Modeling
# 
# # # # # # # # # # # # # # # # # # # # # #
## d = probability density values (pdf = probability distribution function)
# returns height of the distribution curve at specified (input) x value(s)
# in some functions you can specify the mean and sd to tailor it to your curve
# otherwise it treats x as the z score
## p = cumulative probabilities (cdf = cumulative distribution function)
# returns the area below the curve up until the value specified
# i.e. pnorm(1) returns the cumulative probabilities up until 1 in the curve specified
# defaults to mean = 0 and sd = 1
# to find the area above that point pnorm(1, lower.tail=FALSE)
## q = quantile values
# gives you quantiles or critical values
# qnorm(0.95) # p = 0.05, one-tailed (upper) will return 1.644854
# qnorm(0.95, mean = 5, sd =1) # can shift the curve to fit your data like so
## r = random numbers from the distribution
# ref: http://ww2.coastal.edu/kingw/statistics/R-tutorials/prob.html


# # # # # # # # # # # # # # #
#
# VGAM - Pareto Distribution
# 
# # # # # # # # # # # # # # #
## Explore Pareto functions
library("VGAM")
# location and k values are set to the vignette example
# no reason or explanation was provided in help materials
plot(lactation$mean, dpareto(lactation$mean, location = 3, shape = exp(1)),
     xlab = "values of x",
     ylab = "prob. of x",
     main = "Pareto Probability Density Split into 10 Equal Areas")
lines(lactation$mean, dpareto(lactation$mean, location = 1, shape=.5), col = "red")
lines(lactation$mean, dpareto(lactation$mean, location = 1, shape= .2), col = "blue")
# Plot "works", but doesn't provide any helpful information.


# Probability of each x value along a GPD (Generalized Pareto Distribution)
pp = ppareto(lactation$mean, location = 3, shape = exp(1)) # How do you determine location (alpha) and shape (k) values??
range(pp) # 0 and 1
qq = qpareto(pp, location = alpha, shape = k) # Errors out bec. p must be between 0 & 1...conflicts with range above

pareto1() # Pops out pareto formula
# Notes: After reading through help pages, vignettes, and the manual, the documentation
# for this package is unclear and missing crucial definitions of how to use the
# functions.  In addition, when seeking quantiles with qpareto with an input of 
# probabilites from ppareto() the qpareto function returns an error stating that values
# must be between 0 and 1. However range(ppareto()) returns 0 -1.
# For an alternative resource, explore POT package or search other GPD (Generalized Pareto
# Distribution) functions in R.
# plot(lactation$mean, dnorm(lactation$mean))


# # # # # # # # # # # # # # #
#
# POT - Pareto Distribution
# 
# # # # # # # # # # # # # # #
## Explore Peaks Over a Threshold
library("POT")

# Example from news review of package
# ref: http://pot.r-forge.r-project.org/docs/Rnews_2007-1.pdf
# page 34
data("ardieres") # Given test data
plot(density(ardieres))
# Extract independent events from the time series and select
# a suitable threshold for a good asymptotic approximation
tmp <- clust(ardieres, u = 0.85, tim.cond = 7/365,
             clust.max = TRUE) # u = threshold
# cluster stage should balance between bias and variance
# The threshold selection stage is a compromise between
# bias and variance. On one hand, if a too
# high threshold is selected, the bias decreases as the
# asymptotic approximation in equation (1) is good
# enough while the variance increases as there is not
# enough data above this threshold. On the other
# hand, by taking a lower threshold, the variance decreases
# as the number of observations is larger and
# the bias increases as the asymptotic approximation
# becomes poorer.
par(mfrow=c(2,2))
mrlplot(tmp[,"obs"], xlim = c(0.85, 17))
diplot(tmp, u.range = c(0.85, 17))
tcplot(tmp[,"obs"], u.range = c(0.85, 17))
par(mfrow=c(1,1))

events <- clust(ardieres, u = 5,
                tim.cond = 7/365, clust.max = TRUE)
obs <- events[,"obs"]
pwmu <- fitgpd(obs, thresh = 5, "pwmu")
mle <- fitgpd(obs, thresh = 5, "mle")

par(mfrow=c(2,2))
plot(pwmu, npy=2.63)

gpd.pfrl(mle, 0.995, range = c(24.54545, 101.31313),
         conf = 0.95)


## My TURN ##


# Manual suggests determining the threshold using
# tcplot
# mrlplot
# lmomplot
# exiplot
# diplot
# Goal: select enough events to reduce variance, but not too much as that you select
# events coming from a central part of the distribution and induce bias

# These plots *should* help you determine a threshold
par(mfrow=c(1,2))
tcplot(lactation$mean, u.range=c(0,180000))
range(lactation$mean)

# Another method to help determine threshold
# Examples show a clearly decreasing line with the threshold being at the x max
# But ours goes down and then up....so where is the actual threshold
# Also, how do you know which probability you should be using?
# In the example, it seems directed at a point but it is not explained as to why.
mrlplot(lactation$mean, u.range= c(0, quantile(lactation$mean, probs = 0.999908)),
        col = c("green", "black", "green"),
        nt = 200)
mrlplot(NRmature$mean, u.range= c(0, quantile(lactation$mean, probs = 0.999908)),
        col = c("green", "black", "green"),
        nt = 200)

mrlplot(lactation$mean, u.range= c(0, quantile(lactation$mean, probs = 0.995)),
        col = c("green", "black", "green"),
        nt = 200)

mrlplot(NRmature$mean, u.range= c(0, quantile(lactation$mean, probs = 0.995)),
        col = c("green", "black", "green"),
        nt = 200)

## L-moments plot ## Errors out, manual notes that this has poor performance on real data
lmomplot(lactation$mean, u.range= c(1, quantile(lactation$mean, probs = 0.999908)),
         identity=FALSE)
lmomplot(NRmature$mean, u.range= c(1, quantile(lactation$mean, probs = 0.999908)),
         identity=FALSE)


## These should help determine what the proper scale and shape are
# Need a threshold from above though
mom = fitgpd(lactation$mean, threshold= 1, "moments")
par(mfrow=c(2,2))
plot(mom)

par(mfrow=c(2,2))
mle = fitgpd(lactation$mean, threshold=1, est = "mle") # Maximum likelihood estimation
plot(mle)
# How do you interpret these plots?

# # # # # # # # # # # # # # #
#
# Exponential Distribution
# 
# # # # # # # # # # # # # # #

# Plot the distribution of the data
par(mfrow=c(1,1))
# Cow
plot(lactation$mean, dexp(NRmature$mean, rate = 1/mean(lactation$mean)), 
     main = "Exponential",
     ylab = "Density",
     xlab = "TPM"
)

# Get the probabilities for each point along the distribution and test quantiles
pp = pexp(lactation$mean, rate = 1/mean(lactation$mean)) # cow
qexp(0.999908, rate = 1/mean(lactation$mean))

pp = pexp(NRmature$mean, rate = 1/mean(NRmature$mean)) # human
qexp(0.999908, rate = 1/mean(NRmature$mean))

## Test quantiles for sample cut points and calculate possible dilution adjustment factors
## COW ##
qq = qexp(0.999999, rate = 1/mean(lactation$mean))
totalExpr = sum(lactation$mean)
aboveQQ = sum(subset(lactation$mean, lactation$mean > qq)) # 66 high abundance genes
adjFactor = 1/(1-(aboveQQ/totalExpr))

## HUMAN ##
qq = qexp(0.999999, rate = 1/mean(NRmature$mean))
totalExpr = sum(NRmature$mean)
aboveQQ = sum(subset(NRmature$mean, NRmature$mean > qq)) # 308 high abundance genes
adjFactor = 1/(1-(aboveQQ/totalExpr))

# Note:  for pexp and qexp
# lower.tail = TRUE # find prob of less than that number
# lower.tail = FALSE # find prob of greater than that number