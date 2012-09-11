# R script to generate zero counts for missing IDs

p <- read.table(file="uniq_prepubCount.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

l <- read.table(file="uniq_lactCount.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

comb <- merge(p, l, by="GeneID", all=TRUE, suffixes = c("P","L"))

pcomb <- as.data.frame(cbind(comb$GeneID, comb$CountP))
names(pcomb) <- c("GeneID", "Count")
write.table(pcomb, file="pcounts.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

lcomb <- as.data.frame(cbind(comb$GeneID, comb$CountL))
names(lcomb) <- c("GeneID", "Count")
write.table(lcomb, file="lcounts.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)


