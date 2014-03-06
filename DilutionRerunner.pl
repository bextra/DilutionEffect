#!/usr/bin/perl
# DilutionRerunner.pl
use strict; use warnings;

# Objective: run the Dilution Adjusted Model on expression data with empirical 
# threshold generated from thresholdDetermination.R

## Note: This script assumes that RefSeq IDs are not duplicated ##

die "Usage: DilutionRerunner.pl <expression data> <threshold data> <n of reps> \n" unless @ARGV == 3;

my $exprFile = $ARGV[0];
my $reps = $ARGV[2];

# # # # # # # # # # # #
# 
# Parse expression data
#
# # # # # # # # # # # #
 
# Remove the three character file extension from the input
my $exprFileRoot = substr($exprFile, 0, -4); 

# Separate the replicates into their own files
print "Splitting $reps replicates into individual files.\n";
for (my $i = 2; $i <= $reps + 1; $i++) {
	my $rep_tmp = $i - 1;
	`cut -f 1,$i $exprFile > $exprFileRoot$rep_tmp.txt`
}


# # # # # # # # # # # #
# 
# Parse threshold data
#
# # # # # # # # # # # #

my @thresh;
my $count = 0;

open(my $threshFile, "<$ARGV[1]") or die "error opening $ARGV[1] for reading\n";

print "Retrieving threshold values.\n";
while(<$threshFile>){
	my $line = $_;
	next if ($line =~ m/Threshold/i); # skip the header

	my @tmpArray = split(/\s+/, $line); # get the threshold value
	push(@thresh, $tmpArray[1]);
	
	$count++; 	# count the thresholds seen
	if ($count == $reps) { last;} 
	# exits while loop before using the threshold for the mean expression value

}

close $threshFile;

# print array of thresholds if needed
#for (my $i = 0; $i < scalar(@thresh); $i++) {
#	print "$thresh[$i]\n";
#}


# # # # # # # # # # # #
# 
# Apply Dilution Adjustment Model
#
# # # # # # # # # # # #

# Run the dilution_effect.pl on all replicates
run_dilution($reps, $exprFileRoot, @thresh);


# # # # # # # # # # #
#
# SUBROUTINES
#
# # # # # # # # # # #

sub run_dilution {
	my ($reps, $exprFileRoot, @thresh) = @_;

	for (my $i = 1; $i <= $reps; $i++) {

		# run dilution_effect.pl with thresholds
		my $adj_file = join("_", "$exprFileRoot$i", "dilutionOUT.txt"); 
		`./dilution_effect.pl -a $thresh[$i-1] $exprFileRoot$i.txt > $adj_file`; 

		# pseudocount data in preparation of differential expression analysis
		my $pseuct_file = join("_", "$exprFileRoot$i", "dilutionOUT_pseudoct.txt");
		`./pseudocounter.pl $adj_file > $pseuct_file`; 
		
		# prepare file names for final output
		my $final_root = "$exprFileRoot$i";
		$final_root =~ s/Data/Results/;
		my $normEnding = "_norm.txt";
		my $adjEnding  = "_adj.txt";
				
		# separate dilution adjusted and unadjusted data for use with diffexpr_Beck.R 
		`sort $pseuct_file | cut -f 1,2 > $final_root$normEnding`;
		`sort $pseuct_file | cut -f 1,3 > $final_root$adjEnding`;
		
		# clean up unnecessary files
		`rm $adj_file`;
		`rm $pseuct_file`;
	}

}

