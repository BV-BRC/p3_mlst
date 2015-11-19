#! /usr/bin/perl -w

use strict;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage .=   "Usage: $progname COLNUM [options]\n";
$usage .=   "extracts the given (zero-based) COLNUM (or comma-separated COLNUMs)";
$usage .=   "from a tab-delimited table on STDIN and prints to STDOUT. the order \n";
$usage .=   "of the columns in the input list determines the output order. negative \n";
$usage .=   "COLNUM values are accepted.\n";
$usage .=   "       [-i] invert selection and print columns that are not in the input list.\n";
$usage .=   "\n";

my $invert=0;
my $clist;

while (@ARGV) {
	my $arg=shift;
	if ($arg eq '-h' or $arg eq '-help') {
		die "$usage";
	} elsif ($arg eq '-i' or $arg eq '-invert') {
		$invert=1;
	} else {
		$clist=$arg;
	}
}

die "$usage" unless defined $clist; #&& $clist =~ /^\d+$/;
my @incols = split /,/, $clist;

my %colnums=();
foreach my $colnum (@incols) {
	$colnums{$colnum}=1;
}

while (<>) {
	chomp;
	
	if (/^\s*$/) {
		print "\n";
		next;
	}
	
	my @cols = split /\t/;
	my @prnt = ();

	if ($invert) {
	
		my $i=0;
		foreach my $val (@cols) {
			push @prnt, $val unless defined $colnums{$i};
			$i++;
		}

	} else {
	
		foreach my $colnum (@incols) {
			push @prnt, $cols[$colnum] if defined $cols[$colnum];
		}
	
	}
	
	print join("\t", @prnt) . "\n";

}

exit;


