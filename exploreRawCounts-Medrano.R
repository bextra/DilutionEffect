# Explore Medrano lab raw data
# Dec. 5, 2012

# Are genes in the Medrano set integers? i.e. tag counts
tt = read.csv("~/Work/1_Milk/RawData/HolsteinMedranoLab/H2106_D90RNA-Seq.csv")
is.integer(tt$Total.gene.reads)

# Load all replicates
flist = list.files("~/Work/1_Milk/RawData/HolsteinMedranoLab/") 
setwd("~/Work/1_Milk/RawData/HolsteinMedranoLab/")
listMedrano = sapply(flist, read.csv)

row.names(listMedrano)
colnames(listMedrano)

# Make a vector of the cleaned up file names to use as df headers
headers = gsub("(RNA-Seq.csv)", "", colnames(listMedrano), perl = TRUE)

# Get the gene IDs
FeatureIDs = matrix(unlist(listMedrano["Feature.ID", ]), nrow = 24580, ncol = 9)
FeatureIDs = FeatureIDs[,1]

# Get the expression counts for each file
ExpressionCounts = data.frame(matrix(unlist(listMedrano["Total.gene.reads", ]), nrow = 24580, ncol = 9))
names(ExpressionCounts) = headers

# Ultimate data frame
ExpressionCounts = (cbind(FeatureIDs, ExpressionCounts))

par(mfrow=c(3,3), mex = 0.2, oma = c(.3, .3, 1, .3))
mtext("Medrano Lab Raw Expression", side = 3, outer=TRUE)
pie(ExpressionCounts$H2137_D15, labels = "", main = "H2137_Day15", lty = 0)
pie(ExpressionCounts$H2141_D15, labels = "", main = "H2141_Day15", lty = 0)
pie(ExpressionCounts$H2144_D15, labels = "", main = "H2144_Day15", lty = 0)
pie(ExpressionCounts$H2106_D90, labels = "", main = "H2106_Day90", lty = 0)
pie(ExpressionCounts$H2180_D90, labels = "", main = "H2180_Day90", lty = 0)
pie(ExpressionCounts$H2338_D90, labels = "", main = "H2338_Day90", lty = 0)
pie(ExpressionCounts$H2089_D250, labels = "", main = "H2089_Day250", lty = 0)
pie(ExpressionCounts$H2137_D250, labels = "", main = "H2137_Day250", lty = 0)
pie(ExpressionCounts$H2144_D250, labels = "", main = "H2144_Day250", lty = 0)


 
 