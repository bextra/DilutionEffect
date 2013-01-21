
plot(density(cca$Change))
plot(density(hha$Change))


plot(density(cc$Change))
plot(density(hh$Change))

hist(hh$Change,
     breaks= 2000)

hist(cc$Change,
     breaks= 1000)

# Pareto Distribution
library("VGAM")
alpha = 3
k = exp(1)
plot(rpareto(1000, location=alpha, shape=k))
plot(lactation$mean)

library("POT")
?tcplot
par(mfrow=c(1,2))
tcplot(lactation$mean)
tcplot(lactation$mean, u.range=c(0, 1000))

mrlplot(lactation$mean)

tcplot(lactation$mean)

lmomplot(lactation$mean, identify=FALSE)
x = runif(10000)
tcplot(x, u.range=c(0.9, 0.995))
tcplot(x)