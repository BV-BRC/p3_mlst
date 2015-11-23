package VASTR::MLST::LocalMLST;

use strict;
use warnings;

use base 'VASTR::MLST::MLST';

use Capture::Tiny qw/capture/;
use File::HomeDir qw/home/;
use POSIX qw/strftime/;

#use Data::Dumper;


sub new {
	my $class=shift;

	my $home=File::HomeDir->my_home;

	my $tmpdir=shift;
	$tmpdir="$home/mlst-tmp" unless defined $tmpdir && -d $tmpdir;
	
	my $self = {
		CONN   => "$home/mlstdb",
		TMPDIR => $tmpdir
		#CONN   => "/Users/tim/mlstdb",
		#TMPDIR => "/Users/tim/mlst-tmp"
	};
	
	bless($self,$class);
	return $self;
}


sub setCONN {
	my ($self,$loc)=@_;
	$loc = substr($loc, 0, -1) if substr($loc,-1,1) eq "/";
	if (-d $loc) {
		if (-f "$loc/metadata.txt") {
			if (-d "$loc/blastdb/" && -d "$loc/pro/" && -d "$loc/seq/") {
				$self->{'CONN'} = $loc;
			} else {
				warn "unable to set new connection! $loc is missing one or more required subdirs (blastdb/, pro/, or seq/).";
			}
		} else {
			warn "unable to set new connection! $loc/metadata.txt does not exist.";
		}
	} else {
		warn "unable to set new connection! $loc/ is not a dir.";
	}
}


sub getAvailableDatabases {
	my $self=shift;
	
	my $dbs={};
	
	open my $fh, "$self->{CONN}/metadata.txt" or return $self->handleError("getAvailableDatabases", 1, "$self->{CONN}/metadata.txt: $!");
	while (<$fh>) {
		chomp;
		next if /^#/;
		next if /^\s*$/;
		my ($id, $date, $count) = split /\t/;
		my $description = $id;
		$description =~ s/_/ /g;
		$description =~ s/(\d)/\#$1/g;
		$description =~ s/spp/spp\./g;
		$dbs->{$id}=$description;
	}
	close $fh;
	
	$dbs=$self->_fixDbQuery($dbs);
	$self->{DBS}=$dbs;
	return $dbs;
}


sub _fixDbQuery {
	my ($self,$dbs)=@_;
	my %fixedDbs=();
	foreach my $databaseId (keys %$dbs) {
		my $genome = $dbs->{$databaseId};
		if ($genome =~ /^(\w+\s\w+)\sand\s(\w+\s\w+)$/) {
			$genome="$1;$2";
		} elsif ($genome =~ /^(\w+)\s(\w+)/) {
			if ($2 eq "spp") {
				$genome=$1;
			} else {
				$genome="$1 $2";
			}
		}
		$fixedDbs{$databaseId}=$genome;
	}
	return \%fixedDbs;
}


=cut
defined in base module
sub matchDatabases {}
=cut


sub getLoci {
	my ($self,$databaseId) = @_;
	
	my %loci=();
	
	open my $fh, "$self->{CONN}/metadata.txt" or return $self->handleError("getLoci", 2, "$self->{CONN}/metadata.txt: $!");
	while (<$fh>) {
		chomp;
		next if /^#/;
		next if /^\s*$/;
		my @cols = split/\t/;
		next unless shift @cols eq $databaseId;
		shift @cols;	# last retrieved
		shift @cols;	# profile count
		foreach my $locusStr (@cols) {
			my ($name, $data) = split /\|/, $locusStr;
			my ($count, $min, $max) = split /:/, $data;
			$loci{$name}={"min_length"=>$min, "max_length"=>$max, "can_vary"=>($min != $max), "count"=>$count};
		}
		last;
	}
	close $fh;

	return \%loci;
}


sub getAllele {
	my ($self,$databaseId,$locus,$alleleId)=@_;
	
	my $seq="";
	my $capture=0;
	my $file = "$self->{CONN}/seq/$databaseId/$locus.fna";
	# print STDERR "file=$file\n";
	open my $fh, $file or return $self->handleError("getAllele", 5, "$self->{CONN}/seq/$databaseId/$locus.fna: $!");
	while (<$fh>) {
	    chomp;
	    my $line=$_;
	    if ($line =~ /^>/) {
		last unless $seq eq "";
		if ($line eq ">${locus}_${alleleId}"
		    or $line eq ">${locus}-${alleleId}"
		    # the following was from when we had mismatches in the locus
		    # tag in the metadata. Fixed by fixing the update.
		    # or $line =~ /^>[A-Za-z]+_$locus[-_]$alleleId/
		   )
		{
		    $capture=1;
		}
	    } elsif ($capture) {
		$seq .= $line;
	    }
	}
	close $fh;
	return $seq;
}


sub findMatches {
	my ($self,$databaseId,$locus,$seq,$max)=@_;
	
	$max=1 unless defined $max;
	my @matches=();
	
	my $asnfile="$self->{TMPDIR}/${databaseId}-${locus}.asn";
	system("echo \"$seq\" | blastn -task megablast -evalue 0.001 -max_target_seqs $max -outfmt 11 -db $self->{CONN}/blastdb/$databaseId/$locus -query - -out $asnfile");
		
	my $blastResult = capture{ system("blast_formatter -archive $asnfile -outfmt \"6 sseqid mismatch gaps length sseq\"") };
	system("rm $asnfile");
	
	for my $result (split /\n/, $blastResult) {
		my ($aid, $mismatches, $gaps, $length, $aseq) = split /\t/, $result;

		# extract the id from the fasta header
		$aid =~ s/^.+?_(\d+)$/$1/g;

		my $match={
			"database"	 => $databaseId,
			"locus"			 => $locus,
			"allele"		 => $aid,
			"sequence"	 => $aseq,
			"mismatches" => $mismatches,
			"gaps"			 => $gaps,
			"length"		 => $length
		};
		push @matches, $match;
	}
	
	
	return \@matches;
}


sub matchToAllele {
	my ($self,$databaseId,$locus,$seq)=@_;
	
	my @alleles=();	# the return object

	my $asnfile="$self->{TMPDIR}/${databaseId}-${locus}.asn";
	system("echo \"$seq\" | blastn -task megablast -evalue 0.000001 -max_target_seqs 25 -outfmt 11 -db $self->{CONN}/blastdb/$databaseId/$locus -query - -out $asnfile");
	
	my $blastResult = capture{ system("blast_formatter -archive $asnfile -outfmt \"6 sseqid mismatch gaps length sseq btop sstart send qstart qend bitscore slen\"") };
	return \@alleles unless defined $blastResult && $blastResult ne '';
	
	# get the alignments and store them keyed by allele id (locus_allele)
	my $blastAlign = capture{ system("blast_formatter -archive $asnfile -num_alignments 25") };
	my @raw_alignments = split /\n{3}/, $blastAlign;
	my %alignments=();
	foreach my $align (@raw_alignments) {
		my $key = $1 if $align =~ />\s*([A-Za-z0-9_\-]+)/i;
		next unless defined $key;
		next if defined $alignments{$key};
		$alignments{$key}=$align;
	}
	#$blastAlign =~ s/\n/~/g;
	#$blastAlign=$1 if $blastAlign =~ /~{2}(>.+?)~{3}/i;
	#$blastAlign =~ s/~/\n/g;
	system("rm $asnfile");
	
	
	# process the blast results
	my @results=();
	for my $result (split /\n/, $blastResult) {
		my @a = split /\t/, $result;
		my $diffcount = $a[1]+$a[2]+($a[11]-$a[3]);
		push @results, {
			'sseqid'   => $a[0],
			#'mismatch' => $a[1],
			#'gaps'     => $a[2],
			'length'   => $a[3],
			'sseq'  	 => $a[4],
			#'qseq' 	  => $a[6].
			'btop'  	 => $a[5],
			'diffs' 	 => $diffcount,
			'sstart'   => $a[6],
			'send'     => $a[7],
			'qstart'   => $a[8],
			'qend'     => $a[9],
			'bitscore' => $a[10]
		};
	}
	return \@alleles unless scalar @results > 0;
	
	
	# extract the best match(es) from the blast results
	# want to keep equivalent best matches (by num of diffs)
	# initially set to high number to ensure capturing at least the first match
	my $mindiffs = length $seq;
	my $topscore=0;
	#for my $result (sort {$a->{'diffs'} <=> $b->{'diffs'}} @results) {
	for my $result (sort {$b->{'bitscore'} <=> $a->{'bitscore'}} @results) {
		
		last if $result->{'bitscore'} < $topscore;
		$topscore=$result->{'bitscore'};
		
		last if $mindiffs < $result->{'diffs'};
		$mindiffs=$result->{'diffs'};

		# extract the allele id from the fasta header
		my $aid = $result->{'sseqid'};
		if ($aid =~ /^\w+?_(\d+)$/) {
			$aid=$1;
		} elsif ($aid =~ /^\w+?-(\d+)$/) {
			$aid=$1;
		}
		
		my %allele=(
			"database"  => $databaseId,
			"locus"		  => $locus,
			"allele"	  => $aid,
			#'btop'      => $result->{'btop'},
			'bitscore'  => $result->{'bitscore'},
			"diffs"		  => $result->{'diffs'}
		);
		
		if ($result->{'diffs'} > -1) {
		
			#$allele{'soffset'}=$result->{'sstart'}-$result->{'qstart'};
			#$allele{'eoffset'}=$result->{'send'}-$result->{'qend'}+1;
			
			#$allele{'sequence'}=$result->{'sseq'};
			
			$allele{'alignment'}=$alignments{$result->{'sseqid'}};
		
			my @subs;
			
			# account for seqs too short on the 5' end
			
			
			my $pos=1;
			while ($result->{'btop'} =~ /(\d+)([ACTGN-]+)/gi) {
				$pos += $1;
				$pos++ if substr($2,1,1) eq '-';
				push @subs, {"position_in_query"=>$pos, "allele_nuc"=>substr($2,1,1), "query_nuc"=>substr($2,0,1)};
			}
			$allele{"substitutions"} = \@subs if (@subs);
		
		}
		
		push @alleles, \%allele;
		
	}
	
	return \@alleles;
}


sub getSTs {
	my ($self, $databaseId, $profile) = @_;
	# { locus1=>id1, locus2=>id2, locus3=>id3 }
	
	my %matches=();

	my @loci = keys %$profile;
	my %colNums=();
	my $hdrline=1;

	open my $fh, "$self->{CONN}/pro/$databaseId.pro.txt" or return $self->handleError("getSTs", 7, "$self->{CONN}/pro/$databaseId.pro.txt: $!");
	while (<$fh>) {
		chomp;
		my @cols = split /\t/;
		
		my $currentST = shift @cols;
		
		if ($hdrline) {
			# store the column that corresponds to each locus, using the header line
			for (my $i=0; $i<scalar(@cols); $i++) {
				my $locusName=$cols[$i];
				next if $locusName eq 'clonal_complex' or $locusName eq 'mlst_clade';
				$colNums{$locusName}=$i;
			}
			$hdrline=0;
			next;
		}
		
		# test if any of the allele ids fail to match the corresponding allele id in the given profile 
		my $doesMatch=1;
		foreach my $locus (@loci) {
			my $colNum=$colNums{$locus};
			unless ($cols[$colNum] == $profile->{$locus}) {
				$doesMatch=0;
				last;
			}
		}
		
		# store this ST iff all loci in the current profile match the given input profile
		if ($doesMatch) {
		
			my %prohash=();
			foreach my $locus (sort {$colNums{$a} <=> $colNums{$b}} keys %colNums) {
				$prohash{$locus}=$cols[$colNums{$locus}];
			}
			$matches{$currentST} = \%prohash;
		}
		
	}
	close $fh;
	
	return \%matches;
}


sub getProfile {
	my ($self, $databaseId, $st) = @_;
	
	my %profile=();

	my %colNums=();
	my $hdrline=1;
	open my $fh, "$self->{CONN}/pro/$databaseId.pro.txt" or return $self->handleError("getProfile", 9, "$self->{CONN}/pro/$databaseId.pro.txt: $!");
	while (<$fh>) {
		chomp;
		my @cols = split /\t/;
		
		my $currentST = shift @cols;
		
		if ($hdrline) {
			# store the column that corresponds to each locus, using the header line
			for (my $i=0; $i<scalar(@cols); $i++) {
				my $locusName=$cols[$i];
				next if $locusName eq 'clonal_complex' or $locusName eq 'mlst_clade';
				$colNums{$locusName}=$i;
			}
			$hdrline=0;
			next;
		} elsif ($currentST == $st) {
			# get the profile for the ST that matches the query ST
			foreach my $locus (keys %colNums) {
				my $colNum = $colNums{$locus};
				$profile{$locus} = $cols[$colNum];
			}
			last;
		}
	
	}
	return \%profile;
}


sub addAllele {
	my ($self, $databaseId, $locusName, $seq) = @_;
	
	my $seqlen = length $seq;
	
	# get the next available locus id
	my $locusHdrs = capture { system("grep '^>' $self->{CONN}/seq/$databaseId/$locusName.fna") };
	chomp $locusHdrs;
	my @a = split /\n/, $locusHdrs;
	my $lastHdr = pop @a;
	my $newId = pop @{split /_|-/, pop(@a)};
	$newId++;
	
	# add the locus sequence to the file
	$seq =~ s/(.{80})/$1\n/gi;
	$seq = substr($seq,0,-1) if substr($seq, -1, 1) eq "\n";
	system("echo \">${locusName}_${newId}\n$seq\n\" >> $self->{CONN}/seq/$databaseId/$locusName.fna");

	# make the new blastdb
	system("rm $self->{CONN}/blastdb/$databaseId/$locusName.n*");
	my $null = capture{ system("makeblastdb -dbtype nucl -in $self->{CONN}/seq/$databaseId/$locusName.fna -out $self->{CONN}/blastdb/$databaseId/$locusName") };
	
	# update the metadata file
	my $now = strftime "%Y-%b-%d", localtime;

	open my $fh, "$self->{CONN}/metadata.txt" or return $self->handleError("addAllele", 8, "$self->{CONN}/metadata.txt: $!");
	while (<$fh>) {
		chomp;
		my $line=$_;
		
		my @cols = split /\t/, $line;
		
		unless (scalar @cols > 0 && $cols[0] eq $databaseId) {
			print $fh "$line\n";
			next;
		}
		
		$cols[1] = $now;
		for (my $i=3; $i<scalar(@cols); $i++) {
			next unless $cols[$i] =~ /$locusName\|(\d+):(\d+):(\d+)/i;
			my $count =$1;
			my $minlen=$2;
			my $maxlen=$3;
			$count++;
			$minlen=$seqlen if $seqlen < $minlen;
			$maxlen=$seqlen if $seqlen > $maxlen;
			$cols[$i]="$locusName|$count:$minlen:$maxlen";
			last;
		}
		
		print $fh join "\t", @cols;
		print $fh "\n";
	}
	close $fh;
	
	return $newId;
}


sub addST {
	my ($self, $databaseId, $profile) = @_;
	
	my $newST=0;
	my @locusNames = sort keys %$profile;
	
	# update the profile file
	open my $pfh, "$self->{CONN}/pro/$databaseId.pro.txt" or return $self->handleError("addST", 9, "$self->{CONN}/pro/$databaseId.pro.txt: $!");
	while (<$pfh>) {
		chomp;
		my @cols = split /\t/;
		$newST=$cols[0];
	}
	close $pfh;
	$newST++;
	
	my $prostr="$newST";
	foreach my $locus (@locusNames) {
		$prostr .= "\t$profile->{$locus}";
	}
	system("echo \"$prostr\n\" >> $self->{CONN}/pro/$databaseId.pro.txt");

	# update the metadata file
	my $now = strftime "%Y-%b-%d", localtime;

	open my $fh, "$self->{CONN}/metadata.txt" or return $self->handleError("addST", 10, "$self->{CONN}/metadata.txt: $!");
	while (<$fh>) {
		chomp;
		my $line=$_;
		
		my @cols = split /\t/, $line;
		
		unless (scalar @cols > 0 && $cols[0] eq $databaseId) {
			print $fh "$line\n";
			next;
		}
		
		$cols[1] = $now;
		$cols[2] = $cols[2]+1;
		
		print $fh join "\t", @cols;
		print $fh "\n";
	}
	close $fh;
	
	return $newST;
	
}

1;
