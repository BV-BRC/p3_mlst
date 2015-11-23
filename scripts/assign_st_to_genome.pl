#! /usr/bin/perl -w

use strict;
use VASTR::MLST::LocalMLST;
use VASTR::Tools;
use Capture::Tiny qw/capture/;
use File::HomeDir qw/home/;
use Getopt::Long::Descriptive;

use Data::Dumper;

my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;


my($opt, $usage) = describe_options("%c %o",
				    ["retrieves the sequence type (ST) for the fna string on STDIN and writes it to STDOUT."],
				    ["mlst-dir|m=s", "location of local MLST databases"],
				    ["sensitivity|s=s", "sensitivity of response:",{ default => 3 }],
				    ["0: print only if all alleles are full-length and an exact ST is found"],
				    ["1: print if a single ST match is found, even if all alleles are not full-length."],
				    ["2: print if all alleles are full-length and any closely matching STs are found."],
				    ["3: all results, even if no ST is found (default)."],
				    ["db-list|d=s", "comma-separated list of MLST database(s) to query. tries to use the fasta header from STDIN if absent.", { default => "" }],
				    ["extend|e", "extend local MLSTdb with any new alleles and STs."],
				    ["alignments|a", "print all locus alignments to file (default is to print only if there is a mismatch)."],
				    ["threads|T=i", "use N threads to run blast.", { default => 1 }],
				    ["org=s", "Organism name"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 0;

my $dbdir = $opt->mlst_dir;
my $sens   = $opt->sensitivity;
my $dblist =$opt->db_list;
my $addNew =$opt->extend;
my $threads=$opt->threads;

my $org = $opt->org;
my $printAllAlignments = $opt->alignments;

my $mlst = new VASTR::MLST::LocalMLST;
$mlst->setCONN($dbdir) if (defined $dbdir && -d $dbdir);

# set up tmp dir and sequence file
my $ext="fna";
my $query="query";

my $tmp = File::Temp::tempdir(CLEANUP => 1);

# get the MLST databases to query
my @dbs=split /,/, $dblist;
#printf "%s\n", join ",", @dbs;
#exit;

chomp $org if defined $org;
my $orgSlice=$org;
# if necessary, identify the MLST database(s) using the fasta header of the input
unless (scalar @dbs > 0) {
    if (defined $org) {
	my $dbsLoc = $mlst->matchDatabases($org);
	if (scalar @$dbsLoc == 0) {
	    my ($g, $s) = split /\s/, $org;
			$orgSlice = "$g $s";
	    $dbsLoc = $mlst->matchDatabases($orgSlice);
	    if (scalar @$dbsLoc == 0) {
		$orgSlice = "$g";
		$dbsLoc = $mlst->matchDatabases($orgSlice);
	    }
	}
	@dbs=@$dbsLoc;
    }
}
if (@dbs == 0)
{
    print STDERR "No dbs for $org\n";
    exit 0;
}

# blast requires the fna sequence be written to file, so write one to the tmp dir
# also read the organism/genome name from the fasta header
my $gid;
my %stypes=('contig'=>0,'sequence'=>0,'plasmid'=>0,'genome'=>0);
open my $fnaFH, ">", "$tmp/$query.$ext" or die "unable to open fna file for writing: $!";
while (<>) {
	if (/^>/) {
		my @hdr = split /\s{3}/, $_;
		$gid=$hdr[0];
#		$org = $hdr[scalar(@hdr)-1];
#		$org =~ s/\[|\]//g;

		chomp $gid;
		if (@hdr > 1)
		{
		    if ($hdr[1] =~ /contig/i) {
			$stypes{'contig'}=$stypes{'contig'}+1;
		    } elsif ($hdr[1] =~ /plasmid/i) {
			$stypes{'plasmid'}=$stypes{'plasmid'}+1;
		    } elsif ($hdr[1] =~ /chromosome/i) {
			$stypes{'genome'}=$stypes{'genome'}+1;
		    } elsif ($hdr[1] =~ /complete\sgenome/i) {
			$stypes{'genome'}=$stypes{'genome'}+1;
		    } else {
			$stypes{'sequence'}=$stypes{'sequence'}+1;
		    }
		} else {
		    $stypes{'sequence'}=$stypes{'sequence'}+1;
		}
	}
	print $fnaFH "$_";
}
print $fnaFH "\n";
close $fnaFH;

$gid = substr $gid, 1 if defined $gid;

#print Dumper(\@dbs);
#exit;
die "no valid organism, genome, or MLST database defined." unless defined $org or scalar @dbs > 0;


my $stdout;
my $stderr;

# create a blast database of the input sequences
($stdout, $stderr) = capture { system( "makeblastdb -in $tmp/$query.$ext -dbtype nucl -out $tmp/$query" ) };
#print STDERR "$stdout\n";
#print STDERR "$stderr\n";

my %allMatches=();
foreach my $db (@dbs) {
	
	$allMatches{$db}={};

	# get the loci for the current db, keyed by gene symbol:
	# min_length, max_length, can_vary, count
	my $loci = $mlst->getLoci($db);
	
	my %matches=();

	# extract the gene sequences for each locus from the input sequence
	# use the seqs to build a profile
	foreach my $locusName (keys %$loci) {
		
		# store details about the genome sequence we id as the appropriate locus
		my %query = (
			'seq_id'  => '',
			'seq_len' => 0,
			'length'  => 0,
			'start'   => 0,
			'end'     => 0,
			'seq'     => ''
		);
		
		# use the sequence of allele 1 at each locus to extract from genome
		my $locusSeq = $mlst->getAllele($db, $locusName, 1);
#print STDERR "$db: $locusName: $locusSeq\n";
#next;

		# find closest sequence match to MLST allele 1 in the input genome (query)
		($stdout, $stderr) = capture { system( "echo \">mlst|$db|locus|$locusName|   ${locusName}_1   [$org]\n$locusSeq\n\" | blastn -task megablast -evalue .000001 -max_target_seqs 1 -outfmt 11 -db $tmp/$query -query - -out $tmp/$query-$locusName.asn -num_threads $threads" ) };
		# print STDERR "XXX $stderr\n'$locusSeq'\n";
		my $result = capture { system( "blast_formatter -archive $tmp/$query-$locusName.asn -outfmt \"6 sseqid sstart send qstart qend qlen length sseq\"" ) };
		chomp $result;
		
		# unable to find anything resembling an allele for this locus in the input sequence
		# fail silently
		if (!defined $result or $result eq "") {
			$matches{$locusName} = {'query'=>{}, 'matches'=>[{'database'=>$db,'locus'=>$locusName,'allele'=>-1}]};
			next;
		}
		
		my ($seq_id, $seq_start, $seq_end, $a1_start, $a1_end, $a1_len, $match_len, $inseq) = split /\t/, $result;
		
		# blast result returns an empty sequence for unknown reason
		# fail silently
		if (!defined $inseq or $inseq eq "") {
			$matches{$locusName} = {'query'=>{}, 'matches'=>[{'database'=>$db,'locus'=>$locusName,'allele'=>-1}]};
			next;
		}
		
#print "$locusName\t$result\n";
		
		# extract the genome/scaffold/contig sequence from the fna file
		my $genomeseq='';

		my $foundSeq=0;
		open my $seqFH, "$tmp/$query.$ext";
		while (<$seqFH>) {
			chomp;
			if (/^>/) {
				last if $foundSeq;
				my @hdr = split /\s{3}/, $_;
				$hdr[0] =~ s/>//;
				if ($hdr[0] eq $seq_id) {
					$foundSeq=1;
					next;
				}
			} else {
				$genomeseq .= "$_" if $foundSeq;
			}
		}
		close $seqFH;

		%query = (
			'seq_id'     => $seq_id,
			'seq_len'    => length $genomeseq,
			'length'     => $match_len,
			'start'      => $seq_start,
			'start_adj'  => 0,
			'end'        => $seq_end,
			'end_adj'    => 0,
			'target_len' => $a1_len,
			'seq'        => $inseq,
			'full'       => 0
		);
		$query{'full'}=1 if $a1_start==1 && $a1_end==$a1_len;
		
		#
		# HANDLE CASES WHERE INPUT ALLELE LENGTH AND LOCUS ALLELE LENGTH ARE NOT EQUAL!
		#
		
		# adjust start of query, if necessary
		my $leaderseq ='';
		if ($a1_start > 1) {
			my $diffInLen5 = $a1_start-1;
			my $new_start=1;
			
			if ($query{'start'} > $query{'end'}) {
				$new_start = $query{'start'}+$diffInLen5;
				if ($new_start > $query{'seq_len'}) {
					$new_start=$query{'seq_len'};
					$diffInLen5 = $query{'seq_len'}-$query{'start'};
				}
			} else {
				$new_start = $query{'start'}-$diffInLen5;
				if ($new_start<1) {
					$new_start=1;
					$diffInLen5=$query{'start'}-1;
				}
			}
			
			unless ($genomeseq eq '') {
				if ($query{'start'} > $query{'end'}) {
					$leaderseq = substr($genomeseq, $query{'start'}, $diffInLen5);
					$leaderseq = VASTR::Tools::rc($leaderseq);
				} else {
					$leaderseq = substr($genomeseq, $new_start-1, $diffInLen5);
				}
				$query{'start'}=$new_start;
				$query{'length'}=$query{'length'}+length($leaderseq);
				$query{'start_adj'}=$diffInLen5;
			}
			
		}
		
		# adjust end of query, if necessary
		my $trailerseq='';
		if ($a1_end < $a1_len) {
			my $diffInLen3 = $a1_len-$a1_end;
			
			my $new_end=1;
			if ($query{'start'} > $query{'end'}) {
				$new_end=$query{'end'}-$diffInLen3;
				if ($new_end<1) {
					$new_end=1;
					$diffInLen3=$query{'end'}-1;
				}
			} else {
				$new_end=$query{'end'}+$diffInLen3;
				if ($new_end>$query{'seq_len'}) {
					$new_end=$query{'seq_len'};
					$diffInLen3 = $new_end-$query{'end'};
				}
			}

			unless ($genomeseq eq '') {
				if ($query{'start'} > $query{'end'}) {
					$trailerseq = substr($genomeseq, $new_end-1, $diffInLen3);
					$trailerseq = VASTR::Tools::rc($trailerseq);
				} else {
					$trailerseq = substr($genomeseq, $query{'end'}-1, $diffInLen3);
				}
				$query{'end'}=$new_end;
				$query{'length'}=$query{'length'}+length($trailerseq);
				$query{'end_adj'}=$diffInLen3;
			}
			
		}
		
		$query{'seq'} = $leaderseq . $query{'seq'} . $trailerseq;
		$query{'full'}=1 if $query{'length'}==$a1_len;
		
#print Dumper(\%query);

		$matches{$locusName}={'query'=>\%query};
=cut
		match query allele to the closest allele(s) in the MLST db. returns hash array:
			{ 'database'      => $db,
		    'locus'         => $locusName,
		    'allele'        => allele id in profile,
		    'diffs'         => number of diffs between query and match
		  }
		
		if there are differences between the query and the match, also includes:
		    'sequence'      => sequence of matching allele,
		    'alignment'     => blast-formatted pairwise alignment between match and query,
		    'substitutions' => array of substitution definitions
		
		an empty array indicates no close match to query (eval < 1e-6)
=cut
		my $arr = $mlst->matchToAllele($db, $locusName, $query{'seq'});
		
		# add a new allele to the database if the sequence is full but does not match an existing allele exactly
		if ($addNew) {
			if (scalar @$arr == 1) {
				$matches{$locusName}->{'matches'}=$arr;
			} else {
				my $newaid = $mlst->addAllele($db, $locusName, $query{'seq'});
				my $newArr = [{'database'=>$db,'locus'=>$locusName,'allele'=>$newaid,'diffs'=>0}];
				$matches{$locusName}->{'matches'} = [{'database'=>$db,'locus'=>$locusName,'allele'=>$newArr}];
			}
		} else {
			if (scalar @$arr == 0) {
				$matches{$locusName}->{'matches'} = [{'database'=>$db,'locus'=>$locusName,'allele'=>-1}];
			} else {
				$matches{$locusName}->{'matches'}=$arr;
			}
		}

#print Dumper(\%matches);
#last;

	}	# end of loop to find a match for each locus
	
	
	# build an MLST profile for the query sequence
	my @lociNames = sort keys %$loci;
	my $stDeterminants={};	# exact match loci that will be used to query the ST database
	#my @proArr=();	# all allele ids, including inexact matches and missing alleles (-1)
	my %rawProfile=();
	
	# syntax note:
	# loci with no exact or close match in the MLSTdb are represented by a single 
	# array with allele id of -1

	foreach my $locusName (@lociNames) {
	
		if (!defined $matches{$locusName} || scalar(@{$matches{$locusName}->{'matches'}}) == 0 || $matches{$locusName}->{'matches'}->[0]->{'allele'} == -1) {

			# no match at all for this locus
			$rawProfile{$locusName}=-1;
						
		} elsif (scalar(@{$matches{$locusName}->{'matches'}})==1 && $matches{$locusName}->{'query'}->{'full'} && $matches{$locusName}->{'matches'}->[0]->{'diffs'}==0) {
			
			# single exact match for this full-length locus
			# we use these to query for STs

			$rawProfile{$locusName} = $matches{$locusName}->{'matches'}->[0]->{'allele'};
			$stDeterminants->{$locusName} = $matches{$locusName}->{'matches'}->[0]->{'allele'};

		} else {
		
			# multiple exact or close matches for this locus, full-length or not
			
			my @subProf=();
			foreach my $mhash (@{$matches{$locusName}->{'matches'}}) {
				my $prostr=$mhash->{'allele'};
				$prostr .= '*' if $mhash->{'diffs'} > 0;
				$prostr .= '^' unless $matches{$locusName}->{'query'}->{'full'};
				push @subProf, $prostr;
			}
			$rawProfile{$locusName} = join('|',@subProf);
			
		}
		
	}
	
	# find ST for query's profile in MLST database using the stDeterminants hash
	# returns hash of profile hashes keyed by ST
	my $sts={};
	$sts=$mlst->getSTs($db, $stDeterminants) if scalar keys %$stDeterminants > 0;
	
	# add a new ST if all alleles exist but no ST is defined
	if ($addNew && scalar keys %$stDeterminants == scalar @lociNames && scalar keys %$sts == 0) {
		my $newST = $mlst->addST($db, $stDeterminants);
		$sts = {$newST => $stDeterminants};
	}
	
#print Dumper($stDeterminants);
#print Dumper($sts);
	
	my @proArr=();
	foreach my $ln (sort keys %rawProfile) {
		push @proArr, $rawProfile{$ln};
	}
	
	$allMatches{$db} = {
		'seq_id'    => $gid,
		'organism'  => $org,
		'loci'      => \@lociNames,
		'profile'   => \@proArr,
		'STs'       => $sts,
		'seqtypes'  => \%stypes,
		'matchData' => \%matches
	};
	
}

cleanup();
printMatches(\%allMatches);

exit;



sub printMatches {

	my $matches=shift;
	
	foreach my $db (keys %$matches) {
		
		my $line='';

		my $m=$matches->{$db};
		my $mData=$m->{'matchData'};
		
#warn "$sens";
		
		# handle sensitivity of return value
		my $prostr = join ',', @{$m->{'profile'}};
#warn "$prostr";
		next if $sens == 0 && (scalar keys %{$m->{'STs'}} != 1 || $prostr =~ /\*|\^|\|/);
		next if $sens == 1 && scalar keys %{$m->{'STs'}} != 1;
		next if $sens == 2 && $prostr =~ /\^/;
		
		$line .= "$m->{organism}" if defined $m->{'organism'};
		$line .= "\t";
		$line .= "$db\t";
		
		my @sids;
		foreach my $locusName (sort keys %$mData) {
			my $sid = $mData->{$locusName}->{'query'}->{'seq_id'};
			$sid=-1 unless defined $sid;
			push @sids, $sid;
		}
		$line .= join ', ', @sids;
		$line .= "\t";

		my $typeStr='';
		foreach my $stype (keys %{$m->{'seqtypes'}}) {
			my $count = $m->{'seqtypes'}->{$stype};
			$typeStr .= "$count $stype, " unless $count==0;
		}
		$typeStr = substr($typeStr, 0, -2) unless $typeStr eq '';
		$line .= "$typeStr\t";

		my @coords;
		foreach my $locusName (sort keys %$mData) {
			my $s = $mData->{$locusName}->{'query'}->{'start'};
			$s=0 unless defined $s;
			my $e = $mData->{$locusName}->{'query'}->{'end'};
			$e=0 unless defined $e;
			push @coords, "$s-$e";
		}
		$line .= join ', ', @coords;
		$line .= "\t";

		my $lociStr = join ', ', @{$m->{'loci'}} if defined $m->{'loci'};
		$line .= "$lociStr" if defined $lociStr;
		$line .= "\t";
		
		my $profileStr = join ', ', @{$m->{'profile'}} if defined $m->{'profile'};
		$line .= "$profileStr" if defined $profileStr;
		$line .= "\t";
		
		my $stStr='';
		if (defined $m->{'STs'} && scalar keys %{$m->{'STs'}} > 0) {
			$stStr = join ', ', keys %{$m->{'STs'}};
		} elsif (defined $m->{'profile'}) {
			my $pre=0;
			my $post=0;
			my $tot = scalar @{$m->{'profile'}};
			foreach my $id (@{$m->{'profile'}}) {
				if ($id eq '-1') {
					$pre++;
				} elsif ($id =~ /\*|\^|\|/) {
					$post++;
				}
			}
			$stStr = "-$pre";
			$stStr .= ".$post" unless $post == 0;
			$stStr = '0*' if $pre+$post<$tot;
			$stStr = '0' if $pre+$post==0;
		}
		$line .= "$stStr";
		$line .= "\n";
		
		# print the alignment for any inexact matches
		for (my $i=0; $i<scalar(@{$m->{'profile'}}); $i++) {
			my $profileId = $m->{'profile'}->[$i];
			next if $profileId eq '-1';
			
			my $locusName = $m->{'loci'}->[$i];
			
			# profile id for this locus contains one or more inexact matches
			# write the best pairwise alignment to the tmp directory
			printAlignment($mData, $m->{'seq_id'}, $m->{'organism'}, $db, $locusName, $printAllAlignments) if $profileId =~ /\*|\||\^/ or $printAllAlignments;
		}

=cut
		my $exactMatchStr  ='';
		my $inexactMatchStr='';
		my $absentSeqStr   ='';
		
		for (my $i=0; $i<scalar(@{$m->{'profile'}}); $i++) {
			my $profileId = $m->{'profile'}->[$i];
			my $locusName = $m->{'loci'}->[$i];
			
		}
		
		foreach my $locusName (@{$m->{'loci'}}) {
			my $mArr = $m->{'matches'}->{$locusName};
			
			# handle missing sequence
			if ($mArr->[0]->{'allele'} == -1) {
				$absentSeqStr .= "|$locusName";
				next;
			}
			
			# loop through matches and partition into exact or inexact
			foreach my $m (@$mArr) {
				if ($m->{'diffcount'} == 0) {
					$exactMatchStr .= "|$locusName";
				} else {
					$inexactMatchStr .= "|$locusName:$m->{btop}";
					# write the pairwise alignment to the tmp directory
					printAlignment($tmp, $matches->{$db}, $db, $locusName, $m);
				}
			}
			
		}
		
		$exactMatchStr   = substr $exactMatchStr, 1 unless $exactMatchStr eq '';
		$inexactMatchStr = substr $inexactMatchStr, 1 unless $inexactMatchStr eq '';
		$absentSeqStr    = substr $absentSeqStr, 1 unless $absentSeqStr eq '';
		
		$exactMatchStr   =~ s/\|/, /g;
		$inexactMatchStr =~ s/\|/, /g;
		$absentSeqStr    =~ s/\|/, /g;
		
		print "$matches->{$db}->{genome_id}\t$matches->{$db}->{genome}\t$db\t$lociStr\t$profileStr\t$stStr\t";
		print "$exactMatchStr\t$inexactMatchStr\t$absentSeqStr\n";
=cut
		
		print $line;
		
	}
}



sub printAlignment {
	my ($matchData,$gid,$organism,$db,$locusName,$printAll)=@_;
	
#print Dumper($matchData);
#exit;
	
	$printAll=0 unless defined $printAll;
	
	my $outdir = "alignments/";
	mkdir $outdir unless -d $outdir;
	my $path=$organism;
	$path =~ s/\s+/_/g;
	$path =~ s/\///g;
	$path =~ s/://g;
	$outdir .= $path;
	mkdir($outdir) unless -d $outdir;
	
	my $gfn = "${locusName}.${db}.align";
	$gfn =~ s/\///g;
	$gfn =~ s/://g;
	
	# write the pairwise alignments to file in the tmp directory
	open my $outFH, ">", "$outdir/$gfn.txt" or warn "ERROR: something went wrong opening $outdir/$gfn.txt : $!";
	foreach my $m (@{$matchData->{$locusName}->{'matches'}}) {
		next unless $m->{'diffs'} > 0 or $printAll;
		my $sid   = $matchData->{$locusName}->{'query'}->{'seq_id'};
		$sid=-1 unless defined $sid;
		my $s = $matchData->{$locusName}->{'query'}->{'start'};
		$s=0 unless defined $s;
		my $e = $matchData->{$locusName}->{'query'}->{'end'};
		$e=0 unless defined $e;
		print $outFH ">$sid   $locusName:$db:$s-$e ($m->{diffs} diffs)   [$organism]\n";
		print $outFH "$m->{alignment}\n\n";
	}
	close $outFH;
	
}




sub cleanup {
	system( "rm -r $tmp/*.asn" );
	system( "rm -r $tmp/*.nhr" );
	system( "rm -r $tmp/*.nin" );
	system( "rm -r $tmp/*.nsq" );
	system( "rm -r $tmp/$query.$ext" );
}



