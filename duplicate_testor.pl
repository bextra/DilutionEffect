#!/usr/bin/perl
# duplicate_testor.pl
use strict; use warnings;

# The goal is to check whether some short IDs (e.g. NR_028322) 
# could also exist as part of a longer ID (e.g. NM_001005221)

die "usage: dilution_testor.pl <file.txt>\n" unless @ARGV == 1;

my ($samplefile) = @ARGV; #renames ARGV
my %long_id  = ();
my %short_id = ();

while (my $line = <>) {
	next if ($line =~ m/^tracking_id/);
	chomp($line);

	my($id, $freq, $rest) = split(/\s+/, $line); #splits two columns
	$id =~ s/_dup[0-9]$//; #removes duplicate from end of line	

	if    (length $id == 9)  {$short_id{$id} = 1}
	elsif (length $id == 12) {$long_id{$id}  = 1}
	else  {die "File names are irregularly formatted see lines: $line \n";}
}
#print keys(%long_id), "\n"; test that long_ids are stored properly

#test that the short ids do not exist within the long ids
foreach my $key (keys %short_id) {
	foreach my $key2 (keys %long_id) {
		if ($key2 =~ m/$key/) {print "$key\n";}
		else {die "No short IDs match long IDs\n"};
	}
}
exit;
