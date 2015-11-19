#! /usr/bin/perl -w

use strict;

use Capture::Tiny qw/capture/;
use VASTR::MLST::LocalMLST;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "Usage: $progname [options]\n";
$usage   .=   "gets a list of genomes and their appropriate MLST database(s) ";
$usage   .=   "and prints to STDOUT. only genomes that have a known MLST database ";
$usage   .=   "are included in the output.\n";
$usage   .=   "     [-i I] Look for genome fna files at url or directory I (http://brcdownloads.patricbrc.org/patric2/genomes/).\n";
$usage   .=   "     [-f F] Limit to genomes that match string F.\n";
$usage   .=   "     [-o O] Write local MLST databases to directory O (~/mlstdb).\n";
$usage   .=   "\n";


#my $loc='/mnt/crick_storage/brcdownloads/patric2/genomes';
my $indir='http://brcdownloads.patricbrc.org/patric2/genomes/';
my $filter;
my $outdir;

while (@ARGV) {
	my $arg=shift;
	if ($arg eq '-i' or $arg eq '-in' or $arg eq '-input') {
		defined ($indir=shift) or die "FATAL: no input url or directory provided with -i.\n$usage";
	} elsif ($arg eq '-f' or $arg eq '-filter') {
		defined ($filter=shift) or die "FATAL: no filter term provided with -f.\n$usage";
	} elsif ($arg eq '-o' or $arg eq '-out' or $arg eq '-output') {
		defined ($outdir=shift) or die "FATAL: no output dir provided with -o.\n$usage";
	}
}


$filter =~ s/\s/_/g if defined $filter;

my $intype='local';
$intype='remote' if $indir =~ /^(http|ftp):\/\//i;

$indir=substr($indir,0,-1) if substr($indir,-1,1) eq '/';

my $mlst = new VASTR::MLST::LocalMLST;
$mlst->setCONN($outdir) if defined $outdir && -d $outdir;

#print "#genome\tfile\tmlst_databases\n";

if ($intype eq 'local') {

	$indir=substr($indir,0,-1) if substr($indir,-1,1) eq '/';
	
	my $ls = capture { system("ls $indir") };
	#print "$ls\n";
	#die;
	
	chomp $ls;
	
	foreach my $entry (split(/\n/, $ls)) {
		next if $entry eq "../" or $entry eq "./";

		# FILTER FOR SPECIFIC GENOMES
		next unless !defined $filter or $entry =~ /$filter/i;
		
		my $fnaFile;
		if (-d "$indir/$entry") {
			# handle top-level fna files				
			$fnaFile = capture{ system("ls ${indir}/${entry}/*.fna") };
			#my $fnaFile = ${base}/${entry}/${entry}.fna;
			chomp $fnaFile;
		} elsif (-f "${indir}/$entry" and $entry =~ /\.fna$/) {
			# handle top-level fna files
			$fnaFile="${indir}/$entry";
			#print STDERR "$fnaFile\n";
		}
		next unless defined $fnaFile && $fnaFile ne '';
		
		my $genomeName=$entry;
		$genomeName =~ s/\.fna$//i;
		$genomeName =~ s/_/ /g;
		
		#print STDERR "$genomeName\n";
		my $dbs = $mlst->matchDatabases("$genomeName");
		unless (scalar @$dbs == 0) {
			my $dblist = join ",", @$dbs;
			print "$genomeName\t$fnaFile\t$dblist\n";
		}

	}
	
} else {

	my $ls = capture{ system( "curl -l -s $indir/") };
	chomp $ls;

	foreach my $line (split(/\n/, $ls)) {
		#print STDERR "$line\n";
		
		my ($fn,$fnaFile);
		if ($line =~ /<a\shref="(.+?)\/">/g) {
			$fn=$1;
			next if ($fn eq "../" or $fn eq "nohup.out" or $fn eq "run.sh");
			
			#my $fna = capture{ system("curl -s ${base}${fn}/${fn}.fna") };
			$fnaFile="${indir}${fn}/${fn}.fna";
		} else {
			$fn=$line;
			$fnaFile="${indir}/$line";
		}

		# FILTER FOR SPECIFIC GENOMES
		next unless !defined $filter or $fn =~ /$filter/i;
		
		my $genomeName = $fn;
		$genomeName =~ s/_/ /g;
		
		my $dbs = $mlst->matchDatabases("$genomeName");
		unless (scalar @$dbs == 0) {
			my $dblist = join ",", @$dbs;
			print "$genomeName\t$fnaFile\t$dblist\n";
		}
	}
	
}


exit;
