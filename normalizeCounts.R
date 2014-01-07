# normalizeCounts.R

# Objective: normalize all human replicates to have the same total gene expression
    # Ideally, this will allow samples to be more closely compared and identify an accurate and robust adjustment threshold for dilution effect


# Load counts with exprCountsLoad.R


# find total numbers of genes with zero expression for each sample
# these should be similar in order to normalize data
sapply(2:8, function(x) length(which(NRmature[,x] == 0)))

# get the sum of each sample
sampleSums = sapply(2:8, function(x) sum(NRmature[,x]))

# determine what the minimum is
repMin = min(sampleSums)

# divide each other sample by the maximum and record that number
factorAllReps = sampleSums/repMin

options(scipen = 10)
# multiple the expression of each gene by that number
normalizedDF = data.frame(sapply(2:8, function(x) as.integer(NRmature[,x]/factorAllReps[x-1])))
normalizedDF = cbind(NRmature$GeneID, normalizedDF)
names(normalizedDF) = names(NRmature)

normalizedSums = sapply(2:8, function(x) sum(normalizedDF[,x]))

write.table(normalizedDF[,c("GeneID", "X179.C_m2_2")], file="NormMatureRep1.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(normalizedDF[,c("GeneID", "X174.X_m2_0")], file="NormMatureRep2.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(normalizedDF[,c("GeneID", "X359_m0_0")],   file="NormMatureRep3.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(normalizedDF[,c("GeneID", "X153.C_m0_0")], file="NormMatureRep4.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(normalizedDF[,c("GeneID", "X168.R_m2_2")], file="NormMatureRep5.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(normalizedDF[,c("GeneID", "X183.G_m2_2")], file="NormMatureRep6.txt", quote=FALSE, sep = "\t", row.names=FALSE)


# # # # extra code # # # #

# pull the row where this condition is true for the second column
NRmature[which(NRmature[,2] == 0), 2]

# test
testDF = data.frame(ids = c("A", "B", "C"),
                    rep1 = c(1, 5, 2),
                    rep2 = c(1, 1, 2),
                    rep3 = c(2, 1, 2),
                    rep4 = c(2, 2, 2),
                    rep5 = c(2, 1, 1))

sampleSums = sapply(2:6, function(x) sum(testDF[,x]))