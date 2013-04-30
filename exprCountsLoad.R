# ExprCountsLoad.R
# Kristen Beck
# Dec. 6, 2012

# # # # # # # # # # # # # # # # # 
#
# Functions to Load Data 
#
# # # # # # # # # # # # # # # # #
# Needs all numeric vectors in the data frame i.e. no labels
# But not necessarily a class of numeric
quickRowMean = 
    function(df) {
        allRowMeans = sapply(1:nrow(df), function (row) mean(as.integer(df[row,])))
    }


# Function takes a list of files that you want to make a data frame with
# Returns the data frame with the expression averaged across the rows
loadHarhayCounts = 
    function (flist, n) {
        # Usage: loadHarhayCounts(<filelist>, <number of replicates>
        
        listHarhay = sapply(flist, read.table)
        
        # Check point to make sure data loaded correctly
        row.names(listHarhay) = c("GeneID", "RawCounts")
        colnames(listHarhay)
        
        # Make a vector of the cleaned up file names to use as df headers
        headers = gsub("(.txt)", "", colnames(listHarhay), perl = TRUE)
        
        # Get the gene IDs
        GeneIDs = matrix(unlist(listHarhay["GeneID", ]), nrow = 15768, ncol = n) # get all of SORTED gene ids
        GeneIDs = GeneIDs[-1,1] # fix the header and just take one of the gene lists
        
        # Get the expression counts for each file
        ExpressionCounts = data.frame(matrix(unlist(listHarhay["RawCounts", ]), nrow = 15768, ncol = n), row.names = NULL, stringsAsFactors=FALSE)
        ExpressionCounts = droplevels.data.frame(ExpressionCounts[-1,], row.names=NULL, stringsAsFactors=FALSE)
        names(ExpressionCounts) = headers
        
        # Calculate the mean expression for each gene
        ExpressionCounts$mean = quickRowMean(ExpressionCounts)
        
        # Add in the gene IDs to make the ultimate data frame
        ExpressionCounts = (cbind(GeneIDs, ExpressionCounts, row.names = NULL))
        
        return(ExpressionCounts)
    }

# # # # # # # # # # # # # # # # # 
#
# Load Data Sets 
#
# # # # # # # # # # # # # # # # #

### 1. Harhay Data Set ###
# Use just the six replicates that we have already included in the analysis
pflist = list.files("~/Work/1_Milk/DilutionEffect/RNASeqReps-Harhay/", pattern = "^p.+[1-6].txt")
lflist = list.files("~/Work/1_Milk/DilutionEffect/RNASeqReps-Harhay/", pattern = "^l.+[1-6].txt")
setwd("~/Work/1_Milk/DilutionEffect/RNASeqReps-Harhay/")

# Load the data and get the mean expression value across the gene replicates
prepuberty = loadHarhayCounts(pflist, n = 6) # n = number of replicates
lactation = loadHarhayCounts(lflist, n = 6)

### 2. Nomsen-Rivers Data Set ###
CinciData = read.table("~/Work/1_Milk/DilutionEffect/RawData/NomsenRivers_gene_expression/master_all_counts.txt",
                header = TRUE)

# Group 1 samples (true colostrum) 
NRcolostrum = CinciData[,c("X183.G_c2_0", "X187.K_c2_0")]
NRcolostrum$mean = quickRowMean(NRcolostrum)
NRcolostrum = cbind(GeneID = CinciData[,"GeneID"], NRcolostrum)

# Group 2 samples (transitional)
NRtransitional = CinciData[,c("X160.J_c0_2", "X188.L_c2_0", "X179.C_c2_0", "X154.D_c0_0")]
NRtransitional$mean = quickRowMean(NRtransitional)
NRtransitional = cbind(GeneID = CinciData[,"GeneID"], NRtransitional)

# Group 3 samples (mature)
NRmature = CinciData[,c("X179.C_m2_2", "X174.X_m2_0", "X359_m0_0", "X153.C_m0_0", "X168.R_m2_2", "X183.G_m2_2")]
NRmature$mean = quickRowMean(NRmature)
NRmature = cbind(GeneID = CinciData[,"GeneID"], NRmature)

setwd("~/Work/1_Milk/DilutionEffect/RNASeqReps-NR/")

# Write group to a file without the column containing means
write.table(CinciData[,c("GeneID", "X183.G_c2_0")], file="ColostrumRep1.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(CinciData[,c("GeneID", "X187.K_c2_0")], file="ColostrumRep2.txt", quote=FALSE, sep = "\t", row.names=FALSE)

write.table(CinciData[,c("GeneID", "X160.J_c0_2")], file="TransitionalRep1.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(CinciData[,c("GeneID", "X188.L_c2_0")], file="TransitionalRep2.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(CinciData[,c("GeneID", "X179.C_c2_0")], file="TransitionalRep3.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(CinciData[,c("GeneID", "X154.D_c0_0")], file="TransitionalRep4.txt", quote=FALSE, sep = "\t", row.names=FALSE)

write.table(CinciData[,c("GeneID", "X179.C_m2_2")], file="MatureRep1.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(CinciData[,c("GeneID", "X174.X_m2_0")], file="MatureRep2.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(CinciData[,c("GeneID", "X359_m0_0")],   file="MatureRep3.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(CinciData[,c("GeneID", "X153.C_m0_0")], file="MatureRep4.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(CinciData[,c("GeneID", "X168.R_m2_2")], file="MatureRep5.txt", quote=FALSE, sep = "\t", row.names=FALSE)
write.table(CinciData[,c("GeneID", "X183.G_m2_2")], file="MatureRep6.txt", quote=FALSE, sep = "\t", row.names=FALSE)

write.table(NRcolostrum[,c("GeneID", "mean")], file="AveragedColostrum.txt", quote=FALSE, sep= "\t", row.names= FALSE)
write.table(NRtransitional[,c("GeneID", "mean")], file="AveragedTransitional.txt", quote=FALSE, sep= "\t", row.names= FALSE)
write.table(NRmature[,c("GeneID", "mean")], file="AveragedMature.txt", quote=FALSE, sep= "\t", row.names= FALSE)

write.table(lactation[,c("GeneIDs", "mean")], file="AveragedLactation.txt", quote=FALSE, sep= "\t", row.names= FALSE)
write.table(prepuberty[,c("GeneIDs", "mean")], file="AveragedPrepuberty.txt", quote=FALSE, sep= "\t", row.names= FALSE)