#! /usr/bin/perl -w 

use strict;


use VASTR::MLST::LocalMLST;
use VASTR::MLST::PubMLST;
use Data::Dumper;

my $seq="TTTGATACCGTTGCCGAAGGTTTGGGCGAAATTCGTGATTTATTGCGCCGTTATCATCAT
GTCAGCCATGAGTTGGAAAATGGTTCGAGTGAGGCTTTGTTGAAAGAACTCAACGAATTG
CAACTTGAAATCGAAGCGAAGGACGGCTGGAAACTGGATGCGGCAGTCAAGCAGACTTTG
GGGGAACTCGGTTTGCCGGAAAATGAAAAAATCGGCAACCTTTCCGGCGGTCAGAAAAAG
GGCGTCGCCTTGGCTCAGGCTTGGGTGCAAAAGCCCGACGTATTGCTGCTGGACGAGCCG
TCCAACCATTTGGATATCGACGCGATTATTTGGCTGGAAAATCTGCTCAAAGCGTTTGAA
GGCAGCTTGGTTGTGATTACCCACGACCGCCGTTTTTTGGACAATATCGCCACGCGGATT
GTCGAACTCGATC";
$seq =~ s/\n//g;

#my $api = new VASTR::MLST::LocalMLST;
my $api = new VASTR::MLST::PubMLST;

my $dblist = $api->getAvailableDatabases;
die "Unable to access api." unless scalar keys %$dblist > 0;
print Dumper($dblist);
exit;

# get list of available MLST databases
my $dbs = $api->getAvailableDatabases();
#die "Unable to access api." unless scalar @$dbs > 0;
print Dumper($dbs);
exit;

# get list of available MLST databases
my $dbs = $api->matchDatabases("Streptococcus");
die "Unable to access api." unless scalar @$dbs > 0;
print Dumper($dbs);

=cut
my $loci = $api->getGeneLoci($dbs->[0]);
print Dumper($loci);
exit;

my $alleles = $api->matchToAllele($dbs->[0],"abcZ",$seq);
print Dumper($alleles);
exit;

my $profile = $api->getProfile($dbs->[0], 1);
print Dumper($profile);
exit;
=cut

my $loci = $api->getLoci($dbs->[0]);
print Dumper(keys %$loci);

my $profile = $api->getProfile($dbs->[0], 1);
print Dumper($profile);
#delete $profile->{"adk"};

my $sts = $api->getSTs($dbs->[0], $profile);
print Dumper($sts);
exit;

