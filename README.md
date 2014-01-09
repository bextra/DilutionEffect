## Dilution Effect Analysis Pipeline ver 1.0

## A. Author
Kristen Beck (kristenbeck527[at]gmail[dot]com)

## B. Objective
Utilize this package to determine the quantitative threshold denoting high abundance transcripts in RNA-Seq data. Then using this data calculate and apply a dilution adjustment factor to gene expression values.

## C. Pipeline overview

1. Obtain expression data
2. Determine threshold for high abundance genes
3. Use threshold to calculate and apply dilution adjustment
4. Optional: Run dilution adjustment on multiple files

## D. Important notes
1. Must be Unix based OS

## E. Input
Gene expression data in the following format:
`Identifier<TAB>ExpressionValue`

```
GeneID	Expr
ENSG_00000003137	37
ENSG_00000002549	5329
ENSG_00000003056	468

```
Note: each replicate must be in it's own file
	



## Old instructions
1. Extract raw data from Danielle's files either as count (bovine) or FPKM (human)
2. Get thresholds
3. DilutionRerunner.pl
	Runs dilution_effect.pl and pseduocounter.pl
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
		
