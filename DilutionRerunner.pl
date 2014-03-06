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


__END__

# then everything else should work

## Harhay lactating cow RNA-seq count data ##
my $reps = 6; # replicates

# update this
my $path = "~/Work/1_Milk/DilutionEffect/RNASeqReps-Harhay"; # path to directory of data files

# update this
my $fileroot = "lcountsRep"; # beginning of file without replicate number or extension
# thresholds returned from thresholdDetermination.R

# this will be from above
my @thresh = ('0.001734695', '0.001785946', '0.001719321', '0.001771564', '0.001733576', '0.001736869');
run_dilution($reps, $path, $fileroot, @thresh);







# # # # # # # # # # #
#
# SUBROUTINES
#
# # # # # # # # # # #

sub run_dilution {
	my ($reps, $path, $fileroot, @thresh) = @_;

	for (my $i = 1; $i <= $reps; $i++) {
		my $adj_file = "dilOUT_$fileroot$i.txt";
		
		`dilution_effect.pl -a $thresh[$i-1] $path/$fileroot$i.txt > $adj_file`;
		`pseudocounter.pl $adj_file > pseuct_$adj_file`;
		
		my $final_file = "$fileroot$i.txt";
		
		`sort pseuct_$adj_file | cut -f 1,2 > norm_$final_file`;
		`sort pseuct_$adj_file | cut -f 1,3 > adj_$final_file`;
		
		`rm $adj_file`;
		`rm pseuct_$adj_file`;
	}

}





__DATA__
Empirical Threshhold	0.009237

__END__
# #OLD CODE untouched
# 
# for (my $i = 1; $i <= 6; $i++) {
# 	`dilution_effect.pl -a 0.009237 FullSet/lcountsRep$i.txt > FullSet/Adj_lcountsRep$i.txt`;
# 	my $adj_file = "Adj_lcountsRep$i.txt";
# 	`pseudocounter.pl FullSet/$adj_file > FullSet/pseuct_$adj_file`;
# }
# 
# # run the dilution adjusted script on pre-puberty samples w empirical threshold
# for (my $i = 1; $i <= 10; $i++) {
# 	
# 	`dilution_effect.pl -a 0.009237 FullSet/pcountsRep$i.txt > FullSet/Adj_pcountsRep$i.txt`;
# 	# there should be no pcounts that can be adjusted, so use line below
# 	`pseudocounter.pl FullSet/pcountsRep$i.txt > FullSet/pseuct_pcountsRep$i.txt`;
# 	
# 	# if there are, double check data and use line below to pseudoct
# 	# my $adj_file = "Adj_pcountsRep$i.txt"; # these files should be empty
# 	# `pseudocounter.pl FullSet/$adj_file > FullSet/pseuct_$adj_file`;
# }
# 
# # sort check

#OLD CODE untouched




## Nommsen-Rivers mature lactating RNA-seq FPKM data ##
$reps = 6; # replicates
$path = "~/Work/1_Milk/DilutionEffect/RNASeqReps-NR";
$fileroot = "MatureRep"; # beginning of file without replicate number or extension
# thresholds returned from thresholdDetermination.R
@thresh = ('0.001588403', '0.001685752', '0.001579439', '0.001796879', '0.001846405', '0.001697971');
run_dilution($reps, $path, $fileroot, @thresh);


