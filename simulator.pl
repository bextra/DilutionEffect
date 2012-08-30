#!/usr/bin/perl
# simulator.ok
# K. Beck
use strict; use warnings;

die "process id is $$";

### Make the simulated rep 1-3 P and L ###
for (my $i = 1; $i < 3; $i++) {
	simulate_reps ($i); #TO DO: make it only replicate 3
	## TO DO: figure out a way to append the file name
	
}
# TO DO: transliterate addon to file name to process id

#Store the file name in a list maybe?

# # # # # # # # # # #
#
# SUBROUTINES
#
# # # # # # # # # # #

sub simulate_reps{
	my ($i) = @_;

	# get Raw Counts from appropriate samples
	# Two columns: GeneID \t Count
	# R --no-save --no-restore < select_samples.R
	
	# De-duplicate if multiple entries with same gene ID 
	#~/work/mammary_clusters_project/perl_scripts/uniquify_1stcolumn.pl < prepubCount.txt > uniq_prepubCount.txt
	#~/work/mammary_clusters_project/perl_scripts/uniquify_1stcolumn.pl < lactCount.txt > uniq_lactCount.txt
	
	# If identifier missing in one of the samples, it is count=0 (can use merge)
	#R --no-save --no-restore < add_zero_counts.R
	
	#sed 's/NA/0/' pcounts.txt > pcounts_Results.txt
	#sed 's/NA/0/' lcounts.txt > lcounts_Results.txt
	
	# Use R script to generate replicates: simulate_reps.R
	# For each gene, generate two additional replicates using a poisson distribution with lambda=Count
	
	print `sed 's/SAMPLE/pcounts/g' simulate_reps_GENERIC.R > simulate_reps_pcounts.R`;
	print `R --no-save --no-restore < simulate_reps_pcounts.R`;
	print `sed 's/SAMPLE/lcounts/g' simulate_reps_GENERIC.R > simulate_reps_lcounts.R`;
	print `R --no-save --no-restore < simulate_reps_lcounts.R`;
	
	print "This is dollar i $i";

}
__END__
append files to make notes that they are triplicates together



###Adjust those reps###

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

### Calculate the DEseq ###

### Repeat 1000 times ###




__END__
ask keith about the return at the end of a new line

make a shell script to run n of 3 simulations, adjust each and complete DEseq

there should be 1000 of these triplicates total


post process
must be able to deal with NA and NaN from DEseq errors

count the number of genes that have an adj p value less than 0.05