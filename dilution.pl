#!/usr/bin/perl
#
# dilution.pl
#
# A script to generate a dilution effect, to boost expression counts of genes whose expression
# is otherwise 'diluted' by a small number of very highly expressed genes.
#
# Author: Keith Bradnam, Genome Center, UC Davis…based on code originally written by Kristen Beck
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use List::Util qw(sum);
use warnings FATAL => 'all';
use Getopt::Long;

###################################################### 
#              Command-line options                  #
###################################################### 

my ($quantile, $help, $verbose, $normalize, $keep_means) = (0.9995, undef, undef, undef, undef);

my $usage = "
$0 <counts file for condition A in TSV format> <counts file for condition B in TSV format>

optional arguments:
        --quantile   : quantile to determine 'high abundance' genes <default = $quantile>, range 0.0001–0.9999
        --verbose    : turn on extra output (useful for debugging)
        --normalize  : for each replicate, normalize raw counts by total counts for that replicate <default = off>
        --keep_means : keep the intermediate files that contain mean expression across all replicates
        --help
";

GetOptions (
	"quantile=f" => \$quantile,
	"help"       => \$help,
	"verbose"    => \$verbose,
	"normalize"  => \$normalize,
	"keep_means" => \$keep_means
);


die $usage if ($help);
die $usage if ($quantile and ($quantile < 0.0001 or $quantile > 0.9999));
die $usage unless @ARGV == 2;

my ($counts_file_A, $counts_file_B) = @ARGV;





################################################################################
# Global variables
################################################################################

# Will store raw count for each gene in a hash
# hash will use primary key 'A' or 'B' for data from each condition followed by 
# hash keys of 'raw', 'adjusted' and 'pseucounted' for each type of data
my %count_data;

# to get means of normalized counts we will also need total counts for each replicate
# this will also be used later on, so can store in a global hash.
my %total_counts;

# also will want to have counts sorted in an array for each replicate, again tied
# to a primary key representing the condition (A or B)
my %sorted_counts;

# want to also store all replicate names from header row of counts files
my %replicate_names;




################################################################################
# Process counts file to determine mean expression values for each replicate
################################################################################

print "\n# Process counts file to determine mean expression values for each replicate\n";

my ($number_of_replicates_A, $counts_file_with_means_A) = process_counts_file($counts_file_A, 'A');
my ($number_of_replicates_B, $counts_file_with_means_B) = process_counts_file($counts_file_B, 'B');





#################################################################
# Test for difference in average expression and goodness of fit
#################################################################

print "\n# Testing for difference in mean expression levels and goodness of fit between both conditions\n";
print "# Using normalized count data\n" if ($normalize);

my $R_script = "test_difference_in_average_expression.R";
my $command = "Rscript $R_script $counts_file_with_means_A $counts_file_with_means_B";
system($command) && die "Can't run $command";

# can now remove the temp files that contained mean values
unlink($counts_file_with_means_A) unless ($keep_means);
unlink($counts_file_with_means_B) unless ($keep_means);



####################################
# Apply dilution factor
####################################

print "\n# Applying dilution factor to $counts_file_A\n";
apply_dilution_factor($counts_file_A, $number_of_replicates_A, 'A');

print "\n# Applying dilution factor to $counts_file_B\n";
apply_dilution_factor($counts_file_B, $number_of_replicates_B, 'B');

exit;




####################################
#
# S U B R O U T I N E S
#
####################################

sub process_counts_file{

	my ($input_file, $condition) = @_;

	# form output file name and trim any text/tsv file extension
	my $output_file = $input_file;
	$output_file =~ s/\.(txt|tsv|text)$//;
	$output_file .= "_mean_expression.tsv";

	# will determine the number of replicates from the input file
	my $number_of_replicates;

	open(my $in, "<", $input_file)  or die "Can't open $input_file\n";
	
	# read header line separately
	my $header = <$in>;
	chomp($header);
	@{$replicate_names{$condition}} = split(/\t/, $header);

	# now loop over file to process count data and place in hashes
	while(my $line = <$in>){
		chomp($line);

		# Store raw counts in a hash, tied to gene name. Will use this later on
		my ($gene, @f) = split(/\s+/, $line);
		push(@{$count_data{$condition}{'raw'}{$gene}}, @f);

		# push counts into %sorted_counts hash
		for (my $i = 0; $i < @f; $i++){
			push(@{$sorted_counts{$condition}{$i+1}}, $f[$i]);
		}

		# want the number of replicates for later on
		$number_of_replicates = @f;	

		# make running total of counts for each replicate	
		for (my $i = 1; $i <= @f; $i++){
			$total_counts{$condition}{$i} += $f[$i-1];
		}
	}
	close($in);

	# now need to sorted the sorted counts arrays
	foreach my $key (sort {$a <=> $b} keys %{$sorted_counts{$condition}}){		
		@{$sorted_counts{$condition}{$key}} = sort {$a <=> $b} @{$sorted_counts{$condition}{$key}};
	}
	

	# now we have total counts for each replicate we can
	# generate normalized counts from which we can calculate the mean

	open(my $out, ">", $output_file) or die "Can't open $output_file\n";
	print $out "GeneID\tmean\n";
	
	foreach my $gene (sort keys %{$count_data{$condition}{'raw'}}){
			
		# calculate means based on raw counts or normalized counts (if --normalize specified)
		my $mean;
		my @raw_counts = @{$count_data{$condition}{'raw'}{$gene}};

		if($normalize){
			my @normalized_counts;	
			for (my $i = 0; $i < @raw_counts; $i++){
				my $normalized_count = $raw_counts[$i] / $total_counts{$condition}{$i+1};
				push(@normalized_counts, $normalized_count);
			}
			$mean = sprintf ("%.11f", sum(@normalized_counts) / $number_of_replicates) 
		} else{			
			$mean = sprintf ("%.4f",  sum(@raw_counts)        / $number_of_replicates);
		}

		# write output to file
		print $out "$gene\t$mean\n";	
		
	}	
	close($out);
	
	return($number_of_replicates, $output_file);
}





sub apply_dilution_factor{
	my ($counts_file, $number_of_replicates, $condition) = @_;


	####################################
	# Determine and extract thresholds
	####################################
	
	print "\tCalculating thresholds using quantile = $quantile\n";
	


	print "\tDetermining high abundance genes\n";
	for (my $i = 1; $i <= $number_of_replicates; $i++){
		
		# take integer of quantile expression value and look up euqivalent expression at this position
		my $pos = int($quantile * @{$sorted_counts{$condition}{$i}});
		my $threshold_index = ${$sorted_counts{$condition}{$i}}[$pos];
		
		print "\n\t## $counts_file, replicate $i: ${$replicate_names{$condition}}[$i]\n" if ($verbose);
		print "\t## Using quantile $quantile, high abundance genes are those with expression >= $threshold_index\n" if ($verbose);
		
		# Of the total counts for any given replicate, we want to know how many of those
		# represent 'high abundance' or 'low abundance' counts. Will be able to infer the 
		# latter by subtracting the former from the total count
		my $high_abundance_total = 0;
		my $number_of_genes = @{$sorted_counts{$condition}{$i}};
		my $high_abundance_gene_count = $number_of_genes - $pos;
		

		## Calculate sum expression for high abundance genes
		for (my $j = $number_of_genes-1; $j > 0; $j--){
			my $raw_count = ${$sorted_counts{$condition}{$i}}[$j];
			
			if ($raw_count >= $threshold_index){
				print "\t\tGene: ", $j+1, ", expression: $raw_count\n" if ($verbose);
				$high_abundance_total += $raw_count;
			} else{
				# can stop after we have looked at $number_of_genes genes with highest values
				last;
			}
		}
		
		if ($verbose){
			my $low_abundance_total = $total_counts{$condition}{$i} - $high_abundance_total;
			print "\tTotal expression: $total_counts{$condition}{$i}\n";
			print "\tHigh abundance expression: $high_abundance_total\n";
			print "\tLow abundance expression: $low_abundance_total\n";
			print "\tHigh abundance genes: $high_abundance_gene_count\n";
		}
		

		## Calculate dilution proportion and adjustment factor
		my $dilution_proportion = sprintf("%.6f", $high_abundance_total / $total_counts{$condition}{$i});
		my $adjustment_factor = sprintf("%.6f", 1/(1-$dilution_proportion));
		print "\tDilution proportion: $dilution_proportion\n" if ($verbose);
		print "\tAdjustment factor: $adjustment_factor\n" if ($verbose); 

		
		## Adjust low abundance genes and print all data
		foreach my $gene (keys %{$count_data{$condition}{'raw'}}){
			my $raw_count = ${$count_data{$condition}{'raw'}{$gene}}[$i-1];

			my $new_count = $raw_count;
			if ($raw_count <= $threshold_index) {
				$new_count = sprintf("%.0f", $raw_count * $adjustment_factor);
			} 

			# add to adjusted and pseudocounted part of the hash
			$count_data{$condition}{'adjusted'}{$gene}[$i-1]      = $new_count + 1;
			$count_data{$condition}{'pseudocounted'}{$gene}[$i-1] = $raw_count + 1;
		}
	}

	print "\tWriting the adjusted and pseudocounted data to files\n";

	# form output file names
	my $threshold_file = $counts_file;
	my $adjusted_file = $counts_file;
	my $pseudocounted_file = $counts_file;
	$threshold_file =~ s/\.(txt|tsv|text)$/_thresholds\.txt/;
	$adjusted_file =~ s/\.(txt|tsv|text)$/_adjusted\.tsv/;
	$pseudocounted_file =~ s/\.(txt|tsv|text)$/_pseudocounted\.tsv/;


	# now need to write new output files
	open(my $out, ">", $adjusted_file) or die "Can't write to $adjusted_file\n";
	my $header = join("\t", @{$replicate_names{$condition}});
	print $out "$header\n";
	
	foreach my $gene (sort keys %{$count_data{$condition}{'adjusted'}}){
		my $values = join("\t", @{$count_data{$condition}{'adjusted'}{$gene}});
		print $out "$gene\t$values\n";
	}
	close($out);
	
	open($out, ">", $pseudocounted_file) or die "Can't write to $pseudocounted_file\n";
	print $out "$header\n";
	foreach my $gene (sort keys %{$count_data{$condition}{'pseudocounted'}}){
		my $values = join("\t", @{$count_data{$condition}{'pseudocounted'}{$gene}});
		print $out "$gene\t$values\n";
	}
	close($out);
	
	
}