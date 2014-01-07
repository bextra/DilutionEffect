# R script to select appropriate samples

# prepubertal mammary sample, Grab GeneID and Count
p <- read.table(file="/Users/kristenspencer/Work/1_Milk/Analysis-DGL/prepubResults.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

GeneID <- p[p$Library=="BGA146",]$TxID
Count <-  p[p$Library=="BGA146",]$Count

pcount <- as.data.frame(cbind(GeneID, Count))
write.table(pcount, file="prepubCount.txt", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

# lactating mammary sample, Grab GeneID and Count
l <- read.table(file="lactResults.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

GeneID <- l$TxID
Count <-  l$Count

lcount <- as.data.frame(cbind(GeneID, Count))
write.table(lcount, file="lactCount.txt", quote=FALSE, sep="\t", col.names=FALSE
, row.names=FALSE)



