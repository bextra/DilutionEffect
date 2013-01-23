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


ks.test(x=prepuberty$mean, y=lactation$mean)
ks.test(x=NRcolostrum$mean, y=NRmature$mean)

plot(dwilcox())
