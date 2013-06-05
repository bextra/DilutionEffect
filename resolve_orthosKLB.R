###############################################################################
# Adapted from D. Lemay by K. Beck
# 5/24/2013
# resolve multiple human orthologs by taking one with highest sequence identity
###############################################################################

# Biomart convert to human homologs
# Note: Do not deviate from this file format. In Biomart, make sure that columns are in this order!
# Note: Do not deviate from this file format. In Biomart, make sure that columns are in this order!
# Ensembl Gene ID   Ensembl Transcript ID   Description   % Identity with respect to query gene   Homology Type   Human Ensembl Gene ID
# ENSBTAG00000004043    ENSBTAT00000005286	Fc gamma 2 receptor  [Source:RefSeq peptide;Acc:NP_001001138]			
# ENSBTAG00000044129	ENSBTAT00000061027	Alpha-N-acetylgalactosaminide alpha-2,6-sialyltransferase 6  [Source:UniProtKB/Swiss-Prot;Acc:Q08E15]	92	ortholog_one2one	ENSG00000160408

setwd("~/Work/1_Milk/DilutionEffect/DAVID-GO-Input/1_Pre-processing/BtHomologConversion/")
#########################################################################
## Function:
## For ortholog types that are ortholog_one2one, the look up table is has simply the one Human ENSG ID for each ENSBTAT id
## For ortholog types that are ortholog_one2many or many2many, there are multiple mappings of a ENSBTAT ID to Human ENSG IDs 
## and we need to keep just the Human ENSG with the highest sequence identity
########################################################################
selectBestOrtholog =
    function (annot) {
        names(annot) <- c("EnsID", "EnsTransID", "X.Identity", "Homology.Type", "Human.EnsID") # rename columns for easier handling
        
        annot.one2one = data.frame(annot[which((annot$Homology.Type=="ortholog_one2one")|(annot$Homology.Type == "apparent_ortholog_one2one")), ], row.names=NULL)
        annot.many = data.frame(annot[which((annot$Homology.Type=="ortholog_many2many") | (annot$Homology.Type == "ortholog_one2many")), ], row.names=NULL)
        # Note: there are some genes with blank Homology.Type, but those genes also do not have human orthologs and therefore
        # do not need to be processed
        
        num.many = length(unique(annot.many$EnsTransID)) # collapse by EnsTransID because there are more unique terms compared to the EnsID (gene)
        
        # Extract multiple mapping orthologs and highest X identity for each ENSBTAT
        transcript.xidentity = as.matrix(sapply(unique(annot.many$EnsTransID), function(x) {
            # get the names and row number for transcripts that map to multiple orthologs in the main data set
            dupLocations = which(annot$EnsTransID == x) 
            
            # go over each list entry and get the max value
            ptLocation = annot$X.Identity[dupLocations] # get percentage X identity for duplicated ids
            max.val = max(ptLocation) # make a vector of all of the maximum, if there is a tie duplicates remain
        }))
        
        # Convert and cleanup list matrix into a data frame
        transcript.xidentity = data.frame(cbind(row.names(transcript.xidentity), transcript.xidentity[,1]), row.names = NULL, stringsAsFactors=FALSE)
        droplevels.data.frame(transcript.xidentity)
        names(transcript.xidentity) = c("EnsTransID", "X.Identity")
        
        
        # Retrieve the row number for ortholog term w maximum X.Identity
        maxLoc = as.matrix(sapply(1:num.many, function(x)
            which((annot.many$EnsTransID == transcript.xidentity$EnsTransID[x]) & (annot.many$X.Identity == transcript.xidentity$X.Identity[x])) # duplicates are in same cell
        ))
        
        # Process maximum X.Identity ties by taking the first ortholog term
        maxLoc = sapply(1:num.many, function(x) {
            if (length(maxLoc[[x]]) > 1)
                test = maxLoc[[x,1]][1]
            else 
                test = maxLoc[[x]]
        })
        
        # Pull out the multiple mapping orthologs with the highest identity
        merged.many = annot.many[maxLoc,]
        row.names(merged.many) = NULL # reset row names
        
        # Merge the single mapping orthologs with the higest identity of the multiple mapping orthologs
        compactOrthologs = rbind(annot.one2one, merged.many)
        return(compactOrthologs)
    }

runRefPred = function(ref.file, pred.file, export.file) {
    
    annotPRED = read.delim(file=pred.file, header=TRUE, sep="\t", row.names = NULL, stringsAsFactors = FALSE)
    annotPRED = selectBestOrtholog(annot=annotPRED)
    
    annotREF = read.delim(file=ref.file, header=TRUE, sep="\t", row.names = NULL, stringsAsFactors = FALSE)
    annotREF = selectBestOrtholog(annot=annotREF)
    
    allHomologs = rbind(annotREF, annotPRED) # merge predicted and referenced IDs
    
    # only export unique human Ensembl gene IDs
    write.table(unique(allHomologs$Human.EnsID), file = export.file, sep="\t", row.names=FALSE, quote=FALSE)
}


#########################################################################
##
## Run on differentially expressed gene sets
## 
########################################################################
# Background
runRefPred(ref.file="OLBackgroundREF.txt", pred.file="OLBackgroundPRED.txt", export.file="../../BovineQuantThresh-FINAL/HumanHomologs/All16KHsHomologsEmsembl.txt")

# Adjusted Up and Down
runRefPred(ref.file="OLadj_downallREF.txt", pred.file="OLadj_downallPRED.txt", export.file="../../BovineQuantThresh-FINAL/HumanHomologs/adjdownallHsHomologs.txt")
runRefPred(ref.file="OLadj_upallREF.txt", pred.file="OLadj_upallPRED.txt", export.file="../../BovineQuantThresh-FINAL/HumanHomologs/adjupallHsHomologs.txt")

# Normal Up and Down
runRefPred(ref.file="OLnorm_downallREF.txt", pred.file="OLnorm_downallPRED.txt", export.file="../../BovineQuantThresh-FINAL/HumanHomologs/normdownallHsHomologs.txt")
runRefPred(ref.file="OLnorm_upallREF.txt", pred.file="OLnorm_upallPRED.txt", export.file="../../BovineQuantThresh-FINAL/HumanHomologs/normupallHsHomologs.txt")



#### Top 500 Expressed Genes ####
adj500HomologsREF = selectBestOrtholog(annot= read.delim(file="OL-Top500AdjLacREF.txt", header=TRUE, sep="\t", row.names = NULL, stringsAsFactors = FALSE))
write.table(unique(adj500HomologsREF$Human.EnsID), file = "adj500homologs.txt", sep="\t", row.names=FALSE, quote=FALSE)

lac500HomologsREF = selectBestOrtholog(annot= read.delim(file="OL-Top500LacREF.txt", header=TRUE, sep="\t", row.names = NULL, stringsAsFactors = FALSE))
write.table(unique(lac500HomologsREF$Human.EnsID), file = "lac500homologs.txt", sep="\t", row.names=FALSE, quote=FALSE)

prepub500HomologsREF = selectBestOrtholog(annot= read.delim(file="OL-Top500PrepubREF.txt", header=TRUE, sep="\t", row.names = NULL, stringsAsFactors = FALSE))
write.table(unique(prepub500HomologsREF$Human.EnsID), file = "prepub500homologs.txt", sep="\t", row.names=FALSE, quote=FALSE)

# For predicted Ensembl gene IDs since there are so few add manually to gene lists before import to GO

#### The problem may be that a very small portion of the most abundant bovine genes are mapping to human homologs and 
#### therefore the GO data is missing the most important upgregulated and highly expressed transcripts



### Switchers - newly upregulated after adjustment ###
runRefPred(ref.file="OL-upAdj_uniqueREF.txt", pred.file="OL-upAdj_uniquePRED.txt", export.file="Homologs-upAdj_unique.txt")
