#!/usr/bin/perl
# simulator.ok
# K. Beck
use strict; use warnings;

#print "process id is $$";


### Make the simulated rep 1-3 P and L ###
for (my $i = 1; $i <= 3; $i++) {
	simulate_reps ($i);
	# print `rename 's/PID/$i/g' PID*`;
	## TO DO: figure out a way to append the file name
	# maybe transliterate the file name
	# but how do you know what the name is being renamed
	
#	my $A = ($seq) =~ tr/A/A/;

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
	## if you don't have the results de-duplicated and in the correct form
	## then uncomment below, otherwise proceed
	
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

	# fill in for code directory where simulate scripts exist	
	my $code_dir = "/Users/kristenspencer/Work/Code/Milk/Simulate1000K_Code";

	print `sed 's/SAMPLE/pcounts/g' $code_dir/simulate_reps_GENERIC.R > $code_dir/simulate_reps_pcounts.R`;
	print `R --no-save --no-restore < $code_dir/simulate_reps_pcounts.R`;
	print `sed 's/SAMPLE/lcounts/g' $code_dir/simulate_reps_GENERIC.R > $code_dir/simulate_reps_lcounts.R`;
	print `R --no-save --no-restore < $code_dir/simulate_reps_lcounts.R`;

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



__END__
# Stella's rename script
#!/usr/bin/perl

use strict; use warnings;

my ($file, $del, $rename) = @ARGV;
die "usage: $0 <file> <rename what from name> <rename to>\n" unless @ARGV >= 2;

my ($keep1, $keep2) = $file =~ /^(.*)$del(.*)$/;
#print "$file\n$keep1\n$keep2\n";

die "$del does not exists in the file name\n" unless defined($keep1) or defined($keep2);
$rename = "" if not defined($rename);
print "new file name: $keep1$rename$keep2\n";
my $file2 = $keep1 . $rename . $keep2;

my $cmd = "mv $file $file2";
print "$cmd\n";

system($cmd) == 0 or die "renaming file failed: $!\n";