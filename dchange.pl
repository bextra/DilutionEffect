#!/usr/bin/perl
# dchange.pl
# K. Beck & I. Korf
# Jan. 10, 2012
use strict; use warnings;

die "usage: $0 <pre lac> <lac> <threshold>\n" unless @ARGV == 3;

my ($f1, $f2, $THRESHOLD) = @ARGV;

# Read counts
# Return reference to a hash GeneID -> count value
my ($c1, $t1) = counts($f1); # pre-pub
my ($c2, $t2) = counts($f2); # lac

# Normalize counts to have the same totals
# The list with the larger count total will be reduced to meet the smaller total
if ($t1 > $t2) {normalize($c1, $t1 / $t2)}
else           {normalize($c2, $t2 / $t1)}


# Prune genes so both expression values are over the threshold of reliability
# Zero or low counts could be a product of silencing or technical error,
# but by just looking at TPM you are unable to tell
my @killed;
foreach my $k (keys %$c1) {
	if ($c1->{$k} < $THRESHOLD and $c2->{$k} < $THRESHOLD) {
		delete $c1->{$k};
		delete $c2->{$k};
		push @killed, $k;
	}
}

print "killed ", scalar @killed, " genes\n";

# Convert to probabilities of the total expression (also pseudo counts here)
my $exp1 = prob($c1, $t1);
my $exp2 = prob($c2, $t2);

my %h; # histogram of change

my @val;

# Represent the change of each gene 
# as the log 2 fold change lac / pre lac
# or absolute change lac - pre lac
## print "Gene\tEarly\tLate\n"; # checkpoint
foreach my $gene (keys %$exp1) {
	#my $change = log($exp2->{$gene} / $exp1->{$gene}) / log(2);
	my $change = ($exp2->{$gene}) - ($exp1->{$gene});
	#$h{int $change}++; # lower resolution hist
	$h{int $change * 10}++; # higher resolution hist
	push @val, $change;
#	print $gene, "\t", $change, "\n";
	printf "$gene\t%.4f\n", $change;
	## print $gene, "\t", $exp1->{$gene}, "\t", $exp2->{$gene}, "\n"; # checkpoint
}

exit;


# Print the histogram of change values
foreach my $val (sort {$a <=> $b} keys %h) {
	#printf "%d\t%d\n", $val, $h{$val}; # lower resolution hist
	printf "%.1f\t%d\n", $val / 10, $h{$val}; # higher resolution hist
}


print scalar keys %$exp1, " total genes\n"; # that remain

@val = sort {$a <=> $b} @val;
my $med =  $val[@val /2];
print "median fold change $med\n";
print 2 ** (-$med), " fold change downward\n";


# TODO:
# Maybe take out the genes that are responsible for some high proportion (i.e. 
# 15% or 30%$) of the total expression and then get the fold change downward
# from the remaining genes



# # # # # # # # # # # # 
#
#     SUBROUTINES
#
# # # # # # # # # # # # 

sub normalize {
	my ($count, $val) = @_;
	foreach my $k (keys %$count) {
		$count->{$k} /= $val;
	}
}


sub prob {
	my ($c) = @_;
	
	my $total;
	my %count;
	# Pseudocount by 1
	foreach my $k (keys %$c) {
		$count{$k} = $c->{$k} + 1;
		$total += $count{$k};
	}
	
	# Convert to proportions of total expression
	my %f;
	foreach my $k (keys %$c) {
		$f{$k} = $count{$k} #/ $total;
	}
	return \%f;
}


sub counts {
	my ($file) = @_;
	
	my %counts;
	my $total;
	open(IN, $file) or die;
	my $header = <IN>;
	while (<IN>) {
		chomp;
		my ($gene, $n) = split;
		$counts{$gene} = $n;
		$total += $n;
	}
		
	return \%counts, $total;
}