#!/bin/bash
# rerunner.sh
# This script will run a program repeatedly with multiple thresholds.


for filename in *.$1
# * is for anything and $1 is for anything from command line i.e. anyfile.txt

foreach thresh ("-a 0.05" "-a 0.10" "-a 0.15" "-a 0.20" "-a 0.25" "-a 0.30")
 	runprog < lcounts_Rep1.txt.in > lcounts_Rep1.txt.{$thresh}.out
done


do
 	open dilution_effect.pl
 	echo -a $thresh

done

# write as two nested for loops: one will loops over files and the other will 
# loop over the thresholds
