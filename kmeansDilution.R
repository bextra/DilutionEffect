# kmeansDilution.R
# K. Beck
# Jan. 8, 2012

# Objective: Determine if k-means clustering across the log2 fold change between
# non-lactating and lactating states can better describe the dilution effect

# Secondary Objective: Determine best clustering algorithm: PAM or CLARA or mixture

# Load packages
library("fpc") # uses pam and clara
library("cluster")

# Load data
tt = read.table("~/cowchange.txt", header=TRUE, sep="\t")

# Explore clustering algorithms
# pam will be slow on anything greater than 200 pts
# use clara instead

# silhouette rules
# if s(i) = 1, cluster appropriately
# if s(i) < 0, missclassified
# if s(i) = 0, can't tell if classified correctly


# # # # # # # # # # # # # # # # # # # # #  #
#
#                FUNCTIONS
#
# # # # # # # # # # # # # # # # # # # # #  #
printRelevant =
    function (pamobject) {
        cat("Mediods\n", pamobject$medoids, 
            "\nOverall average silhouette\n", pamobject$objective, 
            "\nAver. silhouette by cluster\n", pamobject$silinfo$clus.avg.widths)
    }


pamPlain = 
    function(data, k = 3) {
        pamx = pam(data, k)
        printRelevant(pamx)
        return(pamx)
    }


kRangeClara =
    function(data, range, silrank = "multiasw") {
        pamkx = pamk(data, krange = range, criterion=silrank, usepam=FALSE)
        cat("Optimal cluster size =", pamkx$nc, "\n")
        printRelevant(pamkx$pamobject)
        cat("\nCluster Info\n")
        pamkx$pamobject$clusinfo
        return(pamkx)
    }


claraPlain =
    function(data, k = 3) {
        clarax = clara(data, k = k, samples=50, pamLike=FALSE)        
        # set samples to be drawn from data set, 50 is recommended        
        printRelevant(clarax)
        # plot(clarax)
        return(clarax)
    }

#### ---- Test Partitioning Around Medoids (PAM) ---- ####
pamPlain(tt$Change)

# Notes: Really slow, silhouette for third cluster is low, but within range of 
# sample silhouette values in Reynolds et al (2006) Clustering Rules.

#### ---- Test Partitioning Around Medoids w Estim. n of Clusters ---- ####
set.seed(3)
kRangeClara(tt$Change, range= 2:10)

# Notes: Separates out the 16 most enriched genes and provides the highest silhouette

#### ---- Test Clustering Large Applications (CLARA) ---- ####
set.seed(3)
claraPlain(tt$Change)

# Notes: With k = 3, this should perform the same as pamk, but they aren't and I don't know why.


#### ---- Test k-means ---- ####
fit = kmeans(tt$Change, centers = 3)

aggregate(tt$Change, by=list(fit$cluster), FUN = mean)

# Notes: gives cluster at expected centroids but does not report a goodness
# metric like PAM... too basic for large data?
# Cluster means: -2.76, -0.80, 2.23


#### ---- Play with plots ---- ####
range(tt$Change)
options(scipen = 10)

plot(density(rnorm(n = 100, mean=pamkx$pamobject[[1]])))
plot(pamkx$pamobject$medoids)