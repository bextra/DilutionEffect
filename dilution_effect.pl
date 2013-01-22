#!/usr/bin/perl
# dilution_effect.pl
# K. Beck
use strict; use warnings;
use Getopt::Std;
use vars qw ($opt_h $opt_v $opt_a $opt_d);
getopts('hva:d');

# # # # # # # # # #
# 
# OPTIONS
#
# # # # # # # # # #
my $VERSION = "1.0";
if ($opt_h) {
print "Usage: dilution_effect.pl [options] <arguments...>
	Options:
	-h	help
	-v	version
	-a	dilution-adjusted gene expression data file by expression level with
		  specified threshold from command line i.e. -a 5 <file>
	*NOTE: Gene IDs must be unique. Rmeove duplicates with dupbegone.pl if needed\n";
	exit;
} elsif ($opt_v) {
	print "Version: ", $VERSION, "\n"; 
	exit;
} elsif ((not defined ($opt_a)) & (not defined ($opt_d))) {
	die "Please specify threshold method: adjust with specified threshold (-a) or determine threshold (-d)";
}

my ($filename) = @ARGV; # renames ARGV

# Define global variables
my $header;
my $total_freq;
my %genes_to_frequency; # gene ids and their frequencies
my %genes_to_abundance; # gene ids and their abundance as a proportion
						# of total frequency of all genes

## Get data from file (ref to hash) and cumulative frequency across all genes
load_data($filename);
printf STDERR "Cumulative total expression (all genes): %.4f\n", $total_freq;

# THRESHOLD QUALITY CONTROL CHECKS:
my $threshold;
if ($opt_a) {
	die "Usage: $0 -a <threshold 0.0 to 1.0> <file>\n" unless (($opt_a > 0) and ($opt_a <=1));
	# Need one floating point number between 0 and 1 as command line argument 
	$threshold = $opt_a; # Renames threshold for easy calling
	print STDERR "Specified threshold: $threshold\n";
}

# Determine threshold if option selected
if ($opt_d) {
	$threshold = determine_threshold();
	print STDERR" The determined threshold is $threshold\n";
}


exit; # move eventually

# Adjust frequency of low abundance transcripts
adjust_abundance();

warn "Other arguments were: @ARGV\n\n";


# # # # # # # # # # #
#
# SUBROUTINES
#
# # # # # # # # # # #
sub load_data {
	$header = <>;
	
	# Store frequencies and find frequency total
	while (my $line = <>) {
		next if ($line =~ m/^tracking_id|GeneID/i); # Skips first line
		
		# Store gene IDs and freq in a hash
		my ($gene, $raw_freq) = split(/\s+/, $line);

		$genes_to_frequency{$gene} = $raw_freq;
		
		#takes the total of all frequencies
		$total_freq += $raw_freq;
	}
}

sub adjust_abundance {
	## Use threshold to determine threshold index
	my $threshold_index = $total_freq * $threshold;
	printf STDERR "Threshold index: %.4f\n", $threshold_index;
	
	# Init variables for expression class totals
	my ($high_abundance_total, $low_abundance_total) = (0, 0);

	## Sum counts from low and high abundance genes separately	
	while (my ($gene, $raw_freq) = each %genes_to_frequency) {
		if ($raw_freq >= $threshold_index) {
			$high_abundance_total += $raw_freq;
		} else {
			$low_abundance_total += $raw_freq;
		}
	}
	printf STDERR "Low  abundance total expression: %.4f\n", $low_abundance_total;
	if ($high_abundance_total > 0) { 
		printf STDERR "High abundance total expression: %.4f\n", $high_abundance_total;
	}
	
	## Return data if no genes in high abundance class (typically non-lactating or control)
	if ($high_abundance_total == 0) {
		print "$header"; # Restore original header
		while (my ($gene, $raw_freq) = each %genes_to_frequency) {
			print "$gene\t$genes_to_frequency{$gene}\n"
		}
		die "No values are above set threshold: $threshold\n";
	} 
	
	## Calculate dilution proportion and adjustment factor
	my $dilution_proportion = $high_abundance_total / $total_freq;
	my $adj_factor = 1/(1-$dilution_proportion);
	printf STDERR "Dilution proportion:\t%.6f\n", $dilution_proportion;
	printf STDERR "Adjustment factor:\t%.6f\n", $adj_factor; 

	# Init adj. freq. and high genes counter
	my $dil_adjusted_freq;
	my $highgenes;
	
	# Format header for output file
	my($field1, $field2) = split(" ", $header); # Break up input header
	print STDOUT "$field1\t$field2\tAdj_$field2\n";
	
	## Adjust low abundance genes and print all data
	while (my ($gene, $raw_freq) = each %genes_to_frequency) {
		if ($raw_freq <= $threshold_index) {
			$dil_adjusted_freq = $raw_freq * $adj_factor; # Adjust low abundance
			printf STDOUT "$gene\t%d\t%d\n", $raw_freq, $dil_adjusted_freq;
			# Rounds to the nearest digit
		} else {
			# Print to file if redirected
			printf STDOUT "$gene\t%d\t%d\n", $raw_freq, $raw_freq;
			
			# Print to screen
			# TODO: to screen these will be duplicated...what's a workaround?
			printf STDERR "$gene\t%d\n", $raw_freq;
			$highgenes++;
		}
	}
	print STDERR "There were $highgenes high abundance genes\n";	
}

sub determine_threshold {
	my $tmp;
	
	# Sort gene IDs in order based on expression value- largest to smallest
	my @keys = sort { $genes_to_frequency{$b} <=> $genes_to_frequency{$a} } keys %genes_to_frequency;
	
	# Correct for data sets with odd number of genes
	my $true_length = int((scalar @keys) / 2) * 2;

	# Get the difference between the sorted expression values
	for (my $i = 0; $i < $true_length; $i++) {
		#print "$keys[$i]\t$genes_to_frequency{$keys[$i]}\n"; # checkpoint
 		
 		my $dif = $genes_to_frequency{$keys[$i]} - $genes_to_frequency{$keys[$i+1]};
 		print "$i\t$keys[$i]\t$genes_to_frequency{$keys[$i]}\t$keys[$i+1]\t$genes_to_frequency{$keys[$i+1]}\t$dif\n";
		
		# TODO: now determine some way to pick the empirical threshold in between these

 	}
 	
	# alternative sorting method
	# 	foreach my $key (sort { $genes_to_frequency{$b} <=> $genes_to_frequency{$a} } keys %genes_to_frequency) {
	# 		printf "%s\t%4d\n", $key, $genes_to_frequency{$key};
	#		$dif = $tmp - $genes_to_frequency{$key}; # calculate dif
	# 		$tmp = $genes_to_frequency{$key}; # reset
	# 		
	# 	}
	#     

	return(0.001);

}
__END__
























