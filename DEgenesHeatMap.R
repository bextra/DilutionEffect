source("http://bioconductor.org/biocLite.R")
library("BiocInstaller")
biocLite("AnnotationDbi")
biocLite("DESeq")
library("ggplot2")
library("RColorBrewer")
library("vcd")



DE_adj <- read.delim("~/Work/1_Milk/Diff_Expr_Genes/DE_genes_pseudoct_adj.txt")
View(DE_adj)
head(DE_adj)

d <- ggplot(DE_adj, aes(baseMeanA, adjLog2FoldChange)) + opts(legend.position = "none")
d + stat_binhex(bins = 10) # tiled w hexagons and color for frequency 
d + stat_bin2d(bins = 10)  # ditto but with square tiles

## Not complaining example w DEgenes
ggplot(data= DE_adj, aes(x=baseMeanA, y=baseMeanB)) + scale_x_log10('baseMeanA') + scale_y_log10('baseMeanB') + geom_tile()

ggplot(data= DE_adj, aes(x=baseMeanA, y=baseMeanB)) + scale_x_log10('baseMeanA') + scale_y_log10('baseMeanB') + geom_tile(aes(fill=log2FoldChange))
ggplot(data= DE_adj, aes(x=baseMeanA, y=log2FoldChange)) + scale_x_log10('baseMeanA') +  geom_tile(aes(fill=baseMeanB))
# kind of cool with geom_jitter

## Working Test Example
my.data <- data.frame(test2)
ggplot(my.data, aes(x=x, y=y, z=z, fill=z)) + scale_y_log10('y') + geom_tile()

## Working Test Example
ggplot(my.data, aes(x, y)) +
  geom_tile(aes(fill=z), colour="white") + 
  scale_fill_gradient(low="white", high="black",breaks = seq(min(my.data$z),max(my.data$z),length.out = 5)) + scale_x_log10('x') +    
  scale_y_log10('y') + opts(axis.ticks = theme_blank()) 

# citation("DESeq") #cite a package example
