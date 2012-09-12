#!/usr/bin/perl
use strict;
use warnings;

# Usage: the script is to pseduocount the two columned files that are split
# on white space by adding 1 to every number in the column. It will return 
# the row ID and the pseudocounted data

die "No file specified" unless @ARGV==1;
my ($filename) = @ARGV;

my $header = <>;
print $header;

my $line;
while ($line = <>) {
	my ($gene, $raw_freq, $adj_freq) = split(/\s+/, $line);	
	if ($filename =~ m/lcounts/i) {
		$raw_freq = $raw_freq + 1;
		# $adj_freq = $adj_freq + 1; # uncomment for three columned file
		print "$gene\t$raw_freq\t$adj_freq\n";
	} 
	else {
		$raw_freq = $raw_freq + 1;
		print "$gene\t$raw_freq\n";
	}	
}

exit;

		