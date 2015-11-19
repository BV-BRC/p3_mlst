#! /usr/bin/perl -w

use strict;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage .=   "Usage: $progname [options]\n";
$usage .=   "Adds metadata tags to the MLST table on STDIN and prints to STDOUT. \n";
$usage .=   "       [-h HDR]  use HDR as the column header (tag).\n";
$usage .=   "       [-c NUM]  zero-based index of the new column (last).\n";
$usage .=   "\n";


my $colnum;
my $hdr ='tag';


while (@ARGV) {
	my $arg=shift;
	if ($arg eq '-help') {
		die "$usage";
	} elsif ($arg eq '-h' or $arg eq '-head' or $arg eq '-header') {
		defined ($hdr=shift) or die "$usage";
	} elsif ($arg eq '-c' or $arg eq '-col' or $arg eq '-column') {
		defined ($colnum=shift) or die "$usage";
		die "$usage" unless $colnum > -1;
	}
}

my $firstLine=1;
while (<>) {
	chomp;
	
	if (/^\s*$/) {
		print "\n";
		next;
	}
	
	my @cols = split /\t/;
	
	my $val="MLST.$cols[1].$cols[7]";
	$val=$hdr if $firstLine and $cols[0] =~ /^#/;
	$firstLine=0;
	
	if (defined $colnum and $colnum < scalar @cols) {
	
		if ($colnum==0) {
			# insert the new column on the front of the row
			if ($cols[0] =~ /^#/) {
				$cols[0] =~ s/^#//;
				$val = "\#${val}";
			}
			unshift @cols, $val;
		} else {
			# splice the new column into the row
			splice @cols, $colnum, 0, $val;
		}
		
	} else {
	
		# add the new column on the end of the row
		push @cols, $val;
		
	}
	
	print join("\t", @cols) . "\n";
	
}

exit;

