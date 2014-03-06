#!/usr/bin/perl
# DilutionRerunner.pl
use strict; use warnings;

# Objective: run the Dilution Adjusted Model expression data with empirical 
# threshold generated from thresholdDetermination.R

## Note: This script assumes that RefSeq IDs are not duplicated ##

die "Usage: DilutionRerunner.pl <expression data> <threshold data> <n of reps> \n" unless @ARGV == 3;

my $exprFile = $ARGV[0];
my $reps = $ARGV[2];

## Parse expression data so each replicate is in its own file

# Remove the three character file extension from the input
my $exprFileRoot = substr($exprFile, 0, -4); 

# Separate the replicates into their own files
print "Splitting $reps replicates into individual files.\n";
for (my $i = 2; $i <= $reps + 1; $i++) {
	my $rep_tmp = $i - 1;
	`cut -f 1,$i $exprFile > $exprFileRoot$rep_tmp.txt`
}



my @thresh;
## Process threshold file
my $count = 0;
open(my $threshFile, "<$ARGV[1]") or die "error opening $ARGV[1] for reading\n";
print "Retrieving threshold values.\n";
while(<$threshFile>){
	my $line = $_;
	next if ($line =~ m/Threshold/i); # skip the header

	my @tmpArray = split(/\s+/, $line); # get the threshold value
	push(@thresh, $tmpArray[1]);
	
	# count the thresholds seen and exit loop when all reps have been covered
	$count++; 
	if ($count == $reps) { last;} # exits while loop before using the threshold for the mean expression value
}
close $threshFile;

# print array of thresholds if needed
#for (my $i = 0; $i < scalar(@thresh); $i++) {
#	print "$thresh[$i]\n";
#}


# then everything else should work
## Harhay lactating cow RNA-seq count data ##
## my $reps = 6; # replicates

# update this
## my $path = "~/Work/1_Milk/DilutionEffect/RNASeqReps-Harhay"; # path to directory of data files


## my $fileroot = "lcountsRep"; # beginning of file without replicate number or extension
# thresholds returned from thresholdDetermination.R

# this will be from above

run_dilution($reps, $exprFileRoot, @thresh);




# notes need to fix file names in sub routine

# passing directory/file
# but adj file needs it to be just file

# may have to go back to the original way with path and file name specified

# # # # # # # # # # #
#
# SUBROUTINES
#
# # # # # # # # # # #

sub run_dilution {
	my ($reps, $exprFileRoot, @thresh) = @_;

	for (my $i = 1; $i <= $reps; $i++) {

		my $adj_file = join("_", "$exprFileRoot$i", "dilutionOUT.txt"); 
		`./dilution_effect.pl -a $thresh[$i-1] $exprFileRoot$i.txt > $adj_file`; # everything ok to here


		my $pseuct_file = join("_", "$exprFileRoot$i", "dilutionOUT_pseudoct.txt");
		`./pseudocounter.pl $adj_file > $pseuct_file`; 
		
		my $final_root = "$exprFileRoot$i";
		$final_root =~ s/Data/Results/;
		my $normEnding = "_norm.txt";
		my $adjEnding  = "_adj.txt";
				
		
		`sort $pseuct_file | cut -f 1,2 > $final_root$normEnding`;
		`sort $pseuct_file | cut -f 1,3 > $final_root$adjEnding`;
		
		`rm $adj_file`;
		`rm $pseuct_file`;
	}

}

