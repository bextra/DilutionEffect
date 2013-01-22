#!/usr/bin/perl
# dilution_effect.pl
# by K_Beck
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
my $total_freq;
my %genes_to_frequency; # gene ids and their frequencies
my %genes_to_abundance; # gene ids and their abundance as a proportion
						# of total frequency of all genes

# THRESHOLD CHECKS:
if ($opt_a) {
	die "Usage: $0 -a <threshold 0.0 to 1.0> <file>\n" unless (($opt_a > 0) and ($opt_a <=1));
	# Need one floating point number between 0 and 1 as command line argument 
}

exit; 

my $threshold = $opt_a; #renames threshold for easy calling
print STDERR "The threshold is $threshold\n";

#calls sub to adjust frequency of low abundance transcripts

adjust_abundance() if ($opt_a);

determine_threshold() if ($opt_d);

warn "Other arguments were: @ARGV\n\n";



# # # # # # # # # # #
#
# SUBROUTINES
#
# # # # # # # # # # #

sub adjust_abundance {
	my $header = <>;
	#print "This is the header $header\n"; # checkpoint
	
	#store frequencies and find frequency total
	while (my $line = <>) {
		next if ($line =~ m/^tracking_id/i); #skips first line
		
		#store gene IDs and freq to hash
		my ($gene, $raw_freq) = split(/\s+/, $line);
		#print "$raw_freq\n"; #checkpoint
		$genes_to_frequency{$gene} = $raw_freq;
		
		#takes the total of all frequencies
		$total_freq += $raw_freq;
	}
	printf STDERR "The summative raw frequency of all genes is %.4f\n", $total_freq;
	
	#calculate threshold percent from total and store as threshold index
	#for empirical threshold, determine the desired index and relative threshold
	my $threshold_index = $total_freq * $threshold;
	printf STDERR "The threshold index is %.4f\n", $threshold_index;
	#printf STDOUT "The threshold index is %.4f\n", $threshold_index;
	
	#find the sum of high abundance genes
	my $high_abundance_total = 0;
	my $low_abundance_total = 0;
	
	# my ($high_abundance_total, $low_abundance_total) = (0, 0);

	while (my ($gene, $raw_freq) = each %genes_to_frequency) {
		if ($raw_freq >= $threshold_index) {
			$high_abundance_total += $raw_freq;
		} else {
			$low_abundance_total += $raw_freq;
		}
	}
	
	printf STDERR "The total expression in low abundance is %.4f\n", $low_abundance_total;
	
	if ($high_abundance_total == 0) {
		# will most likely only be pre-puberty set
		print "GeneID\tRaw Expression\n";
		while (my ($gene, $raw_freq) = each %genes_to_frequency) {
			print "$gene $genes_to_frequency{$gene}\n"
		}
		die "No values are above set threshold\n";
	} else { 
	printf STDERR "The total expression in high abundance is %.4f\n", $high_abundance_total;
	}
#calculate dilution proportion and adjustment factor
	my $dilution_proportion = $high_abundance_total / $total_freq;
	my $adj_factor = 1/(1-$dilution_proportion);
	printf STDERR "The dilution proportion is %.6f\n", $dilution_proportion;
	printf STDERR "The adjustment factor is %.6f\n", $adj_factor; 
	#printf STDOUT "The dilution factor is %.6f\n", $dilution_factor;



#apply dilution factor to genes less than threshold index	
	my $dil_adjusted_freq;
	print STDOUT "GeneID\tRawExp\tAdjExp\n";

	my $highgenes;

	while (my ($gene, $raw_freq) = each %genes_to_frequency) {
		if ($raw_freq <= $threshold_index) {
			#adjust the lower abundance genes
			$dil_adjusted_freq = $raw_freq * $adj_factor;
				if ($header =~ m/count/i) {printf STDOUT "$gene\t%d\t%d\n", $raw_freq, $dil_adjusted_freq;
				} else { printf STDOUT "$gene\t%d\t%d\n", $raw_freq, $dil_adjusted_freq;}
		} else {
			printf STDOUT "$gene\t%d\t%d\n", $raw_freq, $raw_freq;
			printf STDERR "$gene\t%d\n", $raw_freq;
			$highgenes++;
			
		}
	}
	print STDERR "There were $highgenes high abundance genes\n";	
}

__END__
