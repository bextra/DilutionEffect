#!/usr/bin/perl

###################################################################
# File = AandNotB.pl
# Author = D. Lemay
# Date = 1/23/08
# Purpose = 
# This script reports the (A and Not B) combination of two files, A and B.
#
# This script takes two filenames as arguments.  Each input
# file contains a list of strings, one per line.  The script forms
# the list of strings that occur in the first file, but not the second
# and prints this list to stdout. 
####################################################################

print STDERR "Opening $ARGV[0] and $ARGV[1]\n";
open(FILE1, $ARGV[0]) || die "Can't open $ARGV[0]: $!\n";
open(FILE2, $ARGV[1]) || die "Can't open $ARGV[1]: $!\n";

%my_hash;
while (<FILE1>) {
    #remove ending whitespace
    s/\s+$//;
    $my_hash{$_} = 1;
    #print "FILE1" . " " . $_ . "\n";
}
close(FILE1);

while (<FILE2>) {
    #remove ending whitespace
    s/\s+$//;
    #print "FILE2" . " " . $_ . "\n";
    $my_hash{$_} = 0;
}
close(FILE2);

foreach $key (keys(%my_hash)) {
    if ($my_hash{$key}) {
	print "$key\n";
    }
}
