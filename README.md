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

## E. Install this project
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
  
This step runs ``dilution_effect.pl`` to adjust the expression data based on the threshold provided from ``thresholdDetermination.R`` on all replicates provided. Note: header is removed when pseudocounting to not interfere with sorting of gene names.

**Algorithm description** The Dilution Adjustment Model uses the threshold provided to determine an index separating "high" and "low" abundance genes. Genes that are above the index are considered high abundance. Their expression values are summed and used in the calculation of the dilution adjustment factor (formula provided in publication).  The adjustment factor is then multiplied to low abundance genes to create the adjusted data set. The adjusted data set is then used for subsequent analysis i.e. a pairwise comparison to another physiological state or a comparison to the raw unadjusted expression data.  Further details and an application of this model are provided here:  
Beck, K. et al. "Adjusting RNA-Seq data to account for the effect of highly abundant transcripts: a case study in milk production" BMC Bioinformatics. 2014 (submitted)


## H. Advanced
### 1. Override default quantile value.  
Quantile is set to 0.9995 by default to isolate only a small subset of genes with the highest expression values. For other data sets, a different quantile may be more appropriate. This can be accomplished from RStudio by storing a different floating point number in the variable `q_defined` before sourcing `thresholdDetermination.R` i.e.  
``> q_defined = 0.9900 # or any desired number``  
Alternatively, you can set this parameter from the command line.  
``$ Rscript thresholdDetermination.R <file1.txt> <file2.txt> <quantile>``  


		
