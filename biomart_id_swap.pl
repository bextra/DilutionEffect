#!/usr/bin/perl
# biomart_id_swap.pl
use strict; use warnings;

# Usage: This will add columns of Biomart ID names and description to a larger
# file by their shared RefSeqIDs

die "Usage: Please specify files to be merged <biomartoutput.txt, DEgenes.txt>" unless @ARGV == 2; 
my ($file1, $file2) = @ARGV;

open (IN1, "<$file1") or die "error reading $file1 for reading\n";
my $header = <IN1>;

# create a hash structure with RefSeqIDs as the keys
my %RefSeq_info; #XXXX comments
while(my $line = <IN1>){
	chomp $line;
	my ($ensembl, $description, $genename, $RefSeqID) = split(/\t/, $line);
#	$ensembl = "NA" if (not defined $ensembl); XXX
	$RefSeq_info{$RefSeqID}{'Ensembl'}         = $ensembl; #hash w ensemble desc genename;
	$RefSeq_info{$RefSeqID}{'Description'}     = $description;
	$RefSeq_info{$RefSeqID}{'Wiki_Gene_Name'}  = $genename;
}
# print scalar keys %RefSeq_info, "\n"; # prints the actual number of keys
# print keys %name_struct, "\n";  prints all RefSeq IDs

close IN1;

warn "Biomart IDs read.\n";
open (IN2, "<$file2") or die "error reading $file2 for reading\n";
my $header2 = <IN2>;
print "RefSeqID\tWikiGeneName\tEnsembleID\tDescription\tbaseMean\tbaseMeanA\tbaseMeanB\tfoldChange\tlog2FoldChange\tpval\tpadj\tadjLog2FoldChange\n";
while (my $line = <IN2>) {
	chomp $line;
	(my @DEstats) = split(/\t/, $line);

	my ($refseq, @rest_of_line) = split(/\t/, $line); #captures first field -RefSeqID	
	my $rest_of_line .= join ("\t", @rest_of_line);
	
	my $output = "$refseq";
	
	if (exists $RefSeq_info{$refseq}){
		$output .= "\t$RefSeq_info{$refseq}{'Wiki_Gene_Name'}\t$RefSeq_info{$refseq}{'Ensembl'}\t$RefSeq_info{$refseq}{'Description'}\t";
	} else{
		$output .= "\tNA\tNA\tNA\t";
	}
	
	$output .= "$rest_of_line";
	print "$output\n";
}
close IN2;

warn "Done.\n";
