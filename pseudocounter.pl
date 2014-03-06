#!/usr/bin/perl
# pseudocounter.pl
use strict;
use warnings;

# Usage: the script is to pseduocount the output of dilution_effect.pl
# If the expression values are adjusted the file should have 3 columns:
# GeneID	RawExpression	AdjustedExpression
# This script returns the Gene ID and the pseudocounted data.

die "No file specified" unless @ARGV==1;
my ($filename) = @ARGV;

my $header = <>;
my @columns = split(/\s+/, $header);
my $n_columns = scalar@columns;

my $line;
while ($line = <>) {
	my ($gene, $raw_freq, $adj_freq) = split(/\s+/, $line);

	if ($n_columns > 2) {
		$raw_freq = $raw_freq + 1; #pseudocount expression values
		$adj_freq = $adj_freq + 1;
		
		print "$gene\t$raw_freq\t$adj_freq\n";
	} else { # used when the input only has one column of expression values
		$raw_freq = $raw_freq + 1;
		print "$gene\t$raw_freq\n";
	}	
}

exit;

		