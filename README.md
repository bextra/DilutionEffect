## Dilution Effect Analysis Pipeline ver 1.0

## A. Author
**Kristen Beck**  
PhD Candidate, [Korf Lab](www.korflab.ucdavis.edu)  
UC Davis Genome Center  
kristenbeck527 [at] gmail [dot] com

## B. Objective
Utilize this package to determine the quantitative threshold denoting high abundance transcripts in RNA-Seq data. Then using this data calculate and apply a dilution adjustment factor to gene expression values.

## C. Pipeline overview

1. Determine threshold for high abundance genes from gene expression data
2. Use threshold to calculate and apply dilution adjustment


## D. Important notes
1. Must be Unix based OS
2. Download RStudio [here](http://www.rstudio.com). It will be your friend when running R scripts.
3. Store expression files in the Data directory

## E. Quick start
TODO add notes here

## F. Install this project
Clone repository to directory of choice from Unix/Linux terminal  
``$ git clone https://github.com/kbeck527/DilutionEffect.git``  

or Download ZIP of repository (link at right)

## F. Input
Gene expression data in the following format:  
`Identifier<TAB>ExpressionValue<TAB>ExpressionValue`  

For example:  
```
GeneID	Rep1	Rep2	mean
GeneA	37		42		39.5
GeneB	529		544		536.5
GeneC	468		462		465
```
**Requirements:** 
- Store expression data files in ``Data/``  
- Final column with mean value across all replicates is required (must be named "mean" as well). The Mann-Whitney-Wilcoxon and Kolmogorov-Smirnov tests are completed using the mean expression value.
- Gene IDs must identical and in the same order for both expression data sets.
- Gene IDs may not be duplicated within a single expression data set.

## G. Detailed steps
###1. Determine threshold to separate high abundance genes  
***Option 1:***  Use RStudio  
i. Open `DilutionAdjModel.Rproj` in RStudio.  
This will load all necessary R scripts and data files into a interactive GUI with a console where you can run the scripts.  
	*Note:* `DilutionAdjModel.Rproj` is configured to run with sample data included in this repository. Update variables for each expression data set (file names as a character string). 
ii. Source script from console pane  
``> source(thresholdDetermination.R)``  
	  
***Option 2:***  Use terminal  
Run `thresholdDetermination.R` from the terminal with the following arguments:  
``$ Rscript thresholdDetermination.R <expressionData1.txt> <expressionData2.txt>``  

###2. Use quantile to get dilution adjustment factor
Using the main input with all replicates in one file, run the following:
``$ DilutionRerunner.pl <expression data> <threshold data> <n of reps>``  
This step will run ``dilution_effect.pl`` to adjust the expression data based on the threshold provided from ``thresholdDetermination.R``. It will also prepare files for subsequent differential expression analysis. Algorithmic details of dilution adjustment are provided in publication.

###3. Some stuff



## H. Advanced
### 1. Override default quantile value.  
Quantile is set to 0.9995 by default to isolate only a small subset of genes with the highest expression values. For other data sets, a different quantile may be more appropriate. This can be accomplished from RStudio by storing a different floating point number in the variable `q_defined` before sourcing `thresholdDetermination.R` i.e.  
``> q_defined = 0.9900 # or any desired number``  
Alternatively, you can set this parameter from the command line.  
``$ Rscript thresholdDetermination.R <file1.txt> <file2.txt> <quantile>``  

### 2. Other stuff
Note sample data provided are human data described in the Beck et al. pubmed id. can be retrieved from the Gene Expression Omnibus at link <> using the GEO Series accession number GSE45669.


## Old instructions (don't look at me Swan)  
1. Extract raw data from Danielle's files either as count (bovine) or FPKM (human)
2. Get thresholds
3. DilutionRerunner.pl
	Runs dilution_effect.pl and pseduocounter.pl - pseudocounter removed header
	Splits file into two files of adjusted and unadjusted
4. Manually move headers to top of file in human data (sorting put them at the bottom)
4. Move dilution outputs to folders separated by species
5. Pseudocount baseline/control sample (prepuberty and colostrum)
		pseudocounter.pl pcountsRep1.txt > psct_pcountsRep1.txt
		pseudocounter.pl ColostrumRep2.txt > psct_ColostrumRep2.txt
6. Move those to the corresponding species folders
7. Run diffexpr_Beck.R on previously indicated folders
8. Follow instructions in Preprocessing/README.txt
	Split DE output into up or down regulated genes
	Add header back to down regulated genes made by tail
	This will be GO input
9. Optional: Split files into unique up or down between adj and unadj
	Analyze gene list
	or put into GO analysis
		
