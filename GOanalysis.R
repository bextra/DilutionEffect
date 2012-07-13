# Gene Ontology Analysis for Differentially Expressed Genes

#Download if needed and load relevant packages
source("http://bioconductor.org/biocLite.R")
library("BiocInstaller")
library("AnnotationDbi")
library("hgu95av2.db") # Affymetrix Human genome U95 Set annotation data for this chip
biocLite("org.Bt.eg.db") #  installs genome wide annotation for bovine
require("org.Bt.eg.db")  #  loads above
library("GO.db")


## Import data structure
DE_adj <- read.delim("~/Work/1_Milk/Diff_Expr_Genes/DE_genes_pseudoct_adj.txt")
View(DE_adj)
head(DE_adj)
DErefseqIDs <- DE_Adj$id #store the RefSeq IDs from DE output

## must convert RefSeq IDs to Entrez Gene IDs
mappedkeys(org.Bt.egREFSEQ2EG) #view all keys in this database
eghash <- as.list(org.Bt.egREFSEQ2EG[DErefseqIDs])
  #pulls the gene ids as well as mrna ids and stores them all in hash
  #
numREFSEQ2EG <- length(eghash) #check number of keys in hash

for (i )


#### HELP TOOLS ####
vignette("AnnotationDbi") # for help
# apropos("function you kind of know the name of") #search function

# AnnotationDbi Help Specific
available.dbschemas() # list potential databases by organisms
available.db0pkgs2()

#quality control step to see the contents of database
qcdata <- capture.output(org.Bt.eg())

