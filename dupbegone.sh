#!/bin/bash
#can use this for Dilution-adjusted Gene Expression Project

#declare USAGE variable
USAGE= "This program will work with RNA-seq formated *.txt typically in two columns."
#print variable on a screen
echo $USAGE

# cuts off the first column, replaces spaces from 
# the end of gene name and the dup from the gene name, sorts, sends the duplicates
# to a new file.txt, the file now contains the names of all the duplicates
cut -f 1 sample4Kristen.txt | sed -e 's/ //g' -e 's/_dup[0-9]*//' | sort | uniq -d > duplicate_ids.txt 

#finds all the duplicate ids in the raw data file and send that list to another
# file which will list the duplicates with the frequency
grep -f duplicate_ids.txt sample4Kristen.txt > sample_duplicates.txt

# finds all the unique genes and their frequency and creates a new file
grep -v -f duplicate_ids.txt sample4Kristen.txt > sample_unique.txt 

#check your work, unique + duplicate = sample line numbers
wc -l *.txt