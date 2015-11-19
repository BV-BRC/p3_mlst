#! /usr/bin/perl -w

use strict;
use LWP::Simple;
use XML::Simple;
use File::HomeDir qw/home/;
use Capture::Tiny qw/capture/;
use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage .=   "Usage: $progname [options]\n";
$usage .=   "Updates a local store with MLST profiles and sequence data. An update ";
$usage .=   "log is written to STDOUT.\n";
$usage .=   "       [-o OUTDIR] target dir for local store (~/mlstdb/). if OUTDIR ";
$usage .=   "already exists, $progname will assume it contains an old version of ";
$usage .=   "the data and will update only those data that have changed. IMPORTANT: ";
$usage .=   "this may lead to loss of any files in the existing dir! use the -b ";
$usage .=   "flag to backup the existing data store before updating.\n";
$usage .=   "       [-b] create a zipped backup of an existing local store (no).\n";
$usage .=   "\n";

# set up the local dir structure
my $home=File::HomeDir->my_home;
my $outdir=$home."/mlstdb";

my $bak=0;
while (@ARGV) {
	my $arg=shift;
	if ($arg eq '-o' or $arg eq '-out' or $arg eq '-outdir') {
		defined ($outdir=shift) or die "FATAL: -o flag given but no output directory provided.\n$usage";
	} elsif ($arg eq '-b' or $arg eq '-bak' or $arg eq '-backup') {
		$bak=1;
	}
}

my $now = localtime;
print "$now\n";
print "updating local store of MLST data at $outdir.\n";

$outdir=substr($outdir, 0,-1) if substr($outdir,-1,1) eq '/';
# backup the existing store first
if ($bak && -d $outdir) {
	my $t=time;
	my $bakfile = "${outdir}-${t}";
	system("zip -rq $bakfile $outdir/");
	print "$outdir/ backed up to $bakfile.zip\n";
}

mkdir $outdir unless -d $outdir;
my $prodir="$outdir/pro";
mkdir $prodir unless -d $prodir;
my $seqdir="$outdir/seq";
mkdir $seqdir unless -d $seqdir;
my $dbdir ="$outdir/blastdb";
mkdir $dbdir unless -d $dbdir;

# read in metadata for existing local data
# any existing metadata file is saved as a bak file
my %metadata=();
if (-f "$outdir/metadata.txt") {
	open my $fh, "$outdir/metadata.txt";
	while (<$fh>) {
		chomp;
		next if /^#/;
		next if /^\s*$/;
		my ($profile, $date, $count, $dstring) = split /\t/,4;
		$metadata{$profile}={"date"=>$date, "count"=>$count, 'dstring'=>$dstring};
	}
	close $fh;
	system("mv $outdir/metadata.txt $outdir/metadata.bak.txt");
}


# get the metadata xml file from pubmlst.org
my $pubOut = get("http://pubmlst.org/data/dbases.xml");
my $pubXml=XMLin($pubOut);
#print Dumper($pubXml);
#die;
my $a = parseMLSTxml($pubXml);

# get the vp metadata xml file from molvis
#my $vpOut = get("http://molvis.vbi.vt.edu/vp/data/dbases.xml");
#my $vpXml=XMLin($vpOut);
#my $b = parseMLSTxml($vpXml);
#my $total   = $a->[0] + $b->[0];
#my $updated = $a->[1] + $b->[1];

my $total   = $a->[0];
my $updated = $a->[1];


# build the new metadata file
print "updating metadata file.\n";
system("touch $outdir/metadata.txt");
open my $mdfh, ">", "$outdir/metadata.txt";
print $mdfh "#database\tlast_retrieved\tprofiles\tlocus|count:minlen:maxlen\n";
foreach my $pid (keys %metadata) {
	printf $mdfh "$pid\t$metadata{$pid}->{date}\t$metadata{$pid}->{count}\t$metadata{$pid}->{dstring}\n";
}
close $mdfh;
print "done.\n";

print "$updated of $total profiles updated.\n";
#print "end of line.\n\n";
exit;


sub parseMLSTxml {
	my $xml=shift;
	
	my $total  =0;
	my $updated=0;
	
	# for each organism
	foreach my $sp (@{$xml->{"species"}}) {

		my $profileId=$sp->{"content"};
		chomp $profileId;
		$profileId =~ s/\s/_/g;
		$profileId =~ s/\#/_/g;
		$profileId =~ s/\.//g;
		$profileId =~ s/\//-/g;
		
		#print "checking $profileId.\n";
		
		$total++;
		
		# update the metadata file with the new retrieval
		my $newDate=$sp->{"mlst"}->{"database"}->{"retrieved"};
		chomp $newDate;
		my $newCount=$sp->{"mlst"}->{"database"}->{"profiles"}->{"count"};
		chomp $newCount;
		
		if ($newCount==0) {
			system("rm -f $prodir/$profileId.pro.txt");
			system("rm -rf $seqdir/$profileId/");
			next;
		}
		
		# only update if the date and count are different
		unless (doUpdate(\%metadata, $profileId, $newDate, $newCount)) {
			print "$profileId is unchanged.\n";
			next;
		}
		
		print "$profileId has been changed.\n";
		print "  updating profile, loci, and blast data.\n";
	
		# get the profile list and save it locally
		# tab-delim file, header is:
		# ST	locus1..locusn	[clonal_complex|mlst_clade]
		if (defined $sp->{'mlst'} and defined $sp->{'mlst'}->{'database'} and defined $sp->{'mlst'}->{'database'}->{'profiles'} and defined $sp->{'mlst'}->{'database'}->{'profiles'}->{'url'}) {
			my $u = $sp->{'mlst'}->{'database'}->{'profiles'}->{'url'};
			chomp $u;
			system("curl -s -S $u -o $prodir/$profileId.pro.txt");
		} else {
			print STDERR "ERROR: unable to get profile for $profileId.\n";
			next;
		}
		
		unless (defined $sp->{"mlst"}->{"database"}->{"loci"} and defined $sp->{"mlst"}->{"database"}->{"loci"}->{"locus"}) {
			print STDERR "ERROR: No locus information available for $profileId!\n";
			next;
		}
		
		my $count = scalar @{$sp->{"mlst"}->{"database"}->{"loci"}->{"locus"}};
		if ($count<1) {
			print STDERR "ERROR: No locus information available for $profileId!\n";
			next;
		}
		
		print "  updating $count loci.\n";
		
		# for each locus in this profile
		mkdir "$seqdir/$profileId";
		my @locusStats=();
		
		foreach my $locus (@{$sp->{"mlst"}->{"database"}->{"loci"}->{"locus"}}) {
	
			# download the latest seq file (fa format; header is locus_id)
			my $locusName=$locus->{"content"};
			chomp $locusName;
			my $locusUrl=$locus->{"url"};
			chomp $locusUrl;
			system("curl -s -S $locusUrl -o $seqdir/$profileId/$locusName.fna");
			
			open my $testFH, "<", "$seqdir/$profileId/$locusName.fna";
			unless (defined $testFH) {
				print STDERR "ERROR: $profileId $locusName: No sequence file!\n";
				next;
			}
			my $test = do { local $/; <$testFH> };
			if ($test =~ /Invalid\slocus/) {
				print STDERR "ERROR: $profileId $locusName: Invalid locus selected!\n";
				next;
			}
			if ($test eq '') {
				print STDERR "ERROR: $profileId $locusName: No content!\n";
				next;
			}

			# calc the min and max lengths and total count for each seq at this locus
			my $stats = calcLocusStats("$seqdir/$profileId/$locusName.fna");
			push @locusStats, "$locusName|$stats";
			
			# make a blast db
			my $null = capture{ system("makeblastdb -dbtype nucl -in $seqdir/$profileId/$locusName.fna -out $dbdir/$profileId/$locusName") };
		}
		
		# update the internal metadata hash
		$metadata{$profileId} = {"date"=>$newDate, "count"=>$newCount, 'dstring'=>join("\t", @locusStats)};
		
		$updated++;
	}
	
	return [$total, $updated];
	
}

sub doUpdate {
	my ($mdHash, $id, $date, $count) = @_;
	return 1 unless defined $mdHash->{$id};
	return 1 unless $mdHash->{$id}->{"count"} == $count;
	return 1 unless $mdHash->{$id}->{"date"} eq $date;
	return 0;
}


sub calcLocusStats {
	my $seqFna=shift;
	open my $fh, "$seqFna" or return "0:0:0";
	
	my $min=1000000;
	my $max=0;
	my $count=0;
	my $seq='';
	while (<$fh>) {
		chomp;
		if (/^>/) {
			if ($seq ne '') {
				my $seqlen = ($seq =~ tr/ACTGactg//);
				$min=$seqlen unless $min < $seqlen;
				$max=$seqlen unless $max > $seqlen;
			}
			$seq='';
			$count++;
		} else {
			$seq .= $_;
		}
	}
	close $fh;

	if ($seq ne '') {
		my $seqlen = ($seq =~ tr/ACTGactg//);
		$min=$seqlen unless $min < $seqlen;
		$max=$seqlen unless $max > $seqlen;
	}

	return "$count:$min:$max";
}

