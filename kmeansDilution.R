library("fpc")

tt = read.table("~/cowchange.txt", header=TRUE, sep="\t")

range(tt$Change)

# pam will be slow on anything greater than 200 pts
# use clara instead

# silhouette rules
# if s(i) = 1, cluster appropriately
# if s(i) < 0, missclassified
# if s(i) = 0, can't tell if classified correctly

pamx = pam(tt$Change, 3)

summary(pamx)

plot(pamx)

pamk(tt$Change, krange = 2:10)

summ
fit = kmeans(tt$Change, centers = 3)

aggregate(tt$Change, by=list(fit$cluster), FUN = mean)


wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")



