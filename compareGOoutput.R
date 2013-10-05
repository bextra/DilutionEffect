# compareGOoutput.R

setwd("Work/1_Milk/DilutionEffect/DAVID-GO-Output/Human-DefaultGOforswitchers/")

switchers = read.delim("upAdj_unique-GOdefaults.txt", header=TRUE, sep = "\t", na.string="NA")
unadjUp   = read.delim("upNorm-GOdefaults.txt", header=TRUE, sep = "\t", na.string="NA")

length(intersect(switchers$Term, unadjUp$Term))