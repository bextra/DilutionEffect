#!/usr/bin/perl
# DilutionRerunner.pl
use strict; use warnings;

# This assumes that RefSeq IDs are not duplicated

# run the dilution adjusted script on lactation samples w empirical threshold
for (my $i = 1; $i <= 10; $i++) {
	`dilution_effect.pl -a 0.009237 FullSet/lcountsRep$i.txt > FullSet/Adj_lcountsRep$i.txt`;
	my $adj_file = "Adj_lcountsRep$i.txt";
	`pseudocounter.pl FullSet/$adj_file > FullSet/pseuct_$adj_file`;
}

# run the dilution adjusted script on pre-puberty samples w empirical threshold
for (my $i = 1; $i <= 10; $i++) {
	
	`dilution_effect.pl -a 0.009237 FullSet/pcountsRep$i.txt > FullSet/Adj_pcountsRep$i.txt`;
	# there should be no pcounts that can be adjusted, so use line below
	`pseudocounter.pl FullSet/pcountsRep$i.txt > FullSet/pseuct_pcountsRep$i.txt`;
	
	# if there are, double check data and use line below to pseudoct
	# my $adj_file = "Adj_pcountsRep$i.txt"; # these files should be empty
	# `pseudocounter.pl FullSet/$adj_file > FullSet/pseuct_$adj_file`;
}


# sort check



__DATA__
Empirical Threshhold	0.009237