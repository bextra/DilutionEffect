# GOanalysis.R
# K. Beck
# Jan. 2, 2012

library("BiocInstaller")
biocLite("biomaRt")
library("biomaRt")

# # # # # # # # # # # # # # # # # 
#
# Get Human Orthologs
#
# # # # # # # # # # # # # # # # #
testIDs = c("NM_181029", "NM_174294", "NM_173929", "NM_174378", "NM_174828",
"NM_174508")

# set up which mart you are using first then list the attributes

listAttributes()
listFilters()
goids = getBM(attributes=c(""))

# probably want goTerm instead...
# where can you enter which type...i.e. bpfat or mffat
# how do you filter by benjamini
# example from biomaRt manual
goids = getBM(attributes=c('entrezgene','go_id'), 
              filters='entrezgene', 
              values=entrez, 
              mart=ensembl)

# # # # # # # # # # # # # # # # # 
#
# Get GO Annotation
#
# # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # 
#
# Convert Gene IDs
#
# # # # # # # # # # # # # # # # #
## RefSeq mRNA to Ensembl Gene IDs


# show all available bioMart databases
listMarts()

# show data sets with an ensembl database
listDatasets(ensembl)

# select the data base we want to work with
ensemblBT = useMart("ensembl", dataset="btaurus_gene_ensembl")
ensemblHS = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# alternatively
# ensembl = useDataset("ensembl", dataset="btaurus_gene_ensembl")


listDatasets(ensembl)
listFilters(ensembl)
# ensembl_gene_id
# ensembl_transcript_id
# refseq_mrna

getBM()
# attributes = vector one wants to retrieve i.e. the output
# filters = vector of filters i.e. the type of input
# values = vector of values for the filters 
# i.e. IDs you want converted
# mart = object with class of Mart, created by useMart()

# Convert Harhay's RefSeq to match Nomsen-River's Ensembl
convertedRStoEns = getBM(attributes= c("refseq_mrna", "ensembl_gene_id"), 
                         filters = "refseq_mrna",
                         values = lactation$GeneIDs,
                         mart = ensemblBT)
# only partial could be converted- 7261 genes
# ENSBTAG prefix is for Ensembl gene IDs for Bos taurus

# Convert Nomsen-Rivers' Ensembl to Match Harhay's RefSeq IDs
convertedEnstoRS = getBM(attributes= c("refseq_mrna"),
                         filters = "ensembl_gene_id",
                         values = CinciData$GeneID,
                         mart = ensemblHS)
# only partial could be converted
## Note: ENSG prefix is for Ensembl gene IDs that are human

?subset
options(scipen = 10)
head(sort(NRmature$mean, decreasing=TRUE), 100)

barplot(NRmature$mean)