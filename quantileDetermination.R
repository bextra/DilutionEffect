# quartileDetermination.R
# K. Beck
# Jan. 15, 2012

# Load data with exprCountsLoad.R
options(scipen= 10)


## Mann-Whitney-Wilcoxon Test
# Null hypothesis: non-lactating and lactating are identical populations
# if p < 0.05 we can reject the null hypothesis
wilcox.test(x=NRcolostrum$mean, y=NRmature$mean) # returns p < 2.2e-16
wilcox.test(x=NRtransitional$mean, y=NRmature$mean) # returns p < 2.2e-16
wilcox.test(x=NRcolostrum$mean, y=NRtransitional$mean) # p-value = 0.01624....ask Danielle about this case???
# therefore non-lac and lac are nonidentical populations

wilcox.test(x=prepuberty$mean, y=lactation$mean)
# returns p < 2.2e-16
# therefore colostrum and mature are nonidentical populations
# ref: http://www.r-tutor.com/elementary-statistics/non-parametric-methods/mann-whitney-wilcoxon-test

# This needs m and n values, but how do you determine those?
plot(dwilcox(lactation$mean))

## Kolmogorov-Smirnov Tests
ks.test(x=prepuberty$mean, y=lactation$mean)
ks.test(x=NRcolostrum$mean, y=NRmature$mean)




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

# returns the x value at this portion of the curve
quantile(rivers, probs=c(0.25, .5, .75, 1))
plot(density(rivers))
quantile(lactation$mean, probs = c(0.05, 0.95, .999908, 1), names=TRUE)
quantile(lactation$mean, probs = c(0.05, 0.9908, .9968, 1), names=TRUE)
quantile(lactation$mean, probs = c(0.05, 0.95), names=TRUE)
# ref: http://ww2.coastal.edu/kingw/statistics/R-tutorials/prob.html


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


### Explore Peaks Over a Threshold
library("POT")
demo(POT)