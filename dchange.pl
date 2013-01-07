#!/usr/bin/perl
use strict; use warnings;

die "usage: $0 <pre lac> <lac> <threshold>\n" unless @ARGV == 3;

my ($f1, $f2, $THRESHOLD) = @ARGV;

# read counts
my ($c1, $t1) = counts($f1);
my ($c2, $t2) = counts($f2);

# reduce counts in the file with larger expression total
# normalize them to have the same totals
if ($t1 > $t2) {normalize($c1, $t1 / $t2)}
else           {normalize($c2, $t2 / $t1)}

# prune genes based on reliability
# both expression values must be over the threshold of reliability
# this will get rid of the zero (pre lac) and zero (lac) values
# which could be a product of silencing or technical error,
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

# convert to probabilities (also pseudo counts here)
my $exp1 = prob($c1, $t1);
my $exp2 = prob($c2, $t2);

my %h; # histogram of change

my @val;

# represent the change of each gene as the log 2 fold change
# lac / pre lac
foreach my $gene (keys %$exp1) {
	my $change = log($exp2->{$gene} / $exp1->{$gene}) / log(2);
	#$h{int $change}++; # lower resolution hist
	$h{int $change * 10}++; # higher resolution hist
	push @val, $change;
}


# print the histogram of change values
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



####
# Sub
####

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
	# pseudocount by 1
	foreach my $k (keys %$c) {
		$count{$k} = $c->{$k} + 1;
		$total += $count{$k};
	}
	
	# convert to proportions of total expression
	my %f;
	foreach my $k (keys %$c) {
		$f{$k} = $count{$k} / $total;
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