#! /usr/bin/perl -w

use strict;

use Capture::Tiny qw/capture/;

my $indir="ftp://ftp.patricbrc.org/patric2/current_release/fna/";

mkdir "genomes/";

my $ls = capture{ system( "curl -l -s $indir/") };
chomp $ls;

foreach my $line (split(/\n/, $ls)) {
	print STDERR "downloading $line\n";
	`curl "${indir}/$line" -o "genomes/$line"`;
}
