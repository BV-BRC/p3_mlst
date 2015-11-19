package VASTR::MLST::PubMLST;

use base 'VASTR::MLST::MLST';
use strict;
use warnings;

use SOAP::Lite;


sub new {
	my $class=shift;
	my $self = new VASTR::MLST::MLST;

	my $soap = SOAP::Lite
		-> uri('http://pubmlst.org/MLST')
		-> proxy('http://pubmlst.org/cgi-bin/mlstdbnet/mlstFetch.pl');
	$self->{CONN}=$soap;
	
	bless($self,$class);
	return $self;
}

sub getAvailableDatabases {
	my $self=shift;
	my $soapResponse = $self->{CONN}->getDatabaseList();
	
	my $dbs={};
	unless ($soapResponse->fault) {
	
		for my $t ($soapResponse->valueof('//database')) {
			$dbs->{$t->{'name'}}=$t->{'description'};
		}
		$dbs=$self->_fixDbQuery($dbs);
		$self->{DBS}=$dbs;
		return $dbs;
		
	} else {
		return $self->handleError("getAvailableDatabases", 1, $soapResponse->faultstring);
	}
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
	
	my $soapResponse = $self->{CONN}->getLocusList($databaseId);
	
	if ($soapResponse->fault) {
		return $self->handleError("getLoci", 2, $soapResponse->faultstring);
	}
	
	my @errors=();
	my %loci=();

	for my $locus ($soapResponse->valueof('//locus')) {

		$loci{$locus} = {"length"=>0, "variable"=>0, "count"=>0};

		my $soapResponse2 = $self->{CONN}->getLocusLength($databaseId, $locus);
		if ($soapResponse2->fault) {
			push @errors, $self->handleError("getLoci", 3, $soapResponse2->faultstring);
		} else {
			$loci{$locus}->{"length"} = $soapResponse2->valueof('//length');
			$loci{$locus}->{"variable"}=1 if ($soapResponse2->valueof('//variable') eq "true");
		}

		my $soapResponse3 = $self->{CONN}->getAlleleCount($databaseId, $locus);
		if ($soapResponse3->fault) {
			push @errors, $self->handleError("getLoci", 4, $soapResponse3->faultstring);
		} else {
			$loci{$locus}->{"count"} = $soapResponse3->result();
		}

	}
	
	if (scalar @errors > 0) {
		$loci{"errorCode"}=\@errors;
	}
	
	return \%loci;
	
}


sub getAllele {
	my ($self,$databaseId,$locus,$alleleId)=@_;
	
	my $seq="";
	my $soapResponse = $self->{CONN}->getAlleleSequence($databaseId, $locus, $alleleId);

	unless ($soapResponse->fault) {
		$seq .= $soapResponse->result();
		return $seq;
	} else {
		return $self->handleError("getAllele", 5, $soapResponse->faultstring);
	}
}


sub findMatches {
	my ($self,$databaseId,$locus,$seq,$max)=@_;
	
	$max=1 unless defined $max;
	my @matches=();
	my $soapResponse = $self->{CONN}->locusBlast($databaseId,$locus,$seq,$max);

	if ($soapResponse->fault) {
		return $self->handleError("findMatches", 6, $soapResponse->faultstring);
	} else {
		for my $t ($soapResponse->valueof('//blastMatch')) {
			my $match={
				"database"	 => $databaseId,
				"locus"			 => $locus,
				"allele"		 => $t->{'id'},
				"sequence"	 => '',
				"mismatches" => $t->{'mismatches'},
				"gaps"			 => $t->{'gaps'},
				"length"		 => $t->{'alignment'}
			};
			if ($match->{'mismatches'}==0 && $match->{'gaps'}==0) {
				$match->{'sequence'} = $seq;
			} else {
				$match->{'sequence'} = $self->getAllele($databaseId,$locus,$match->{'allele'});
			}
			push @matches, $match;
		}
	}
	return \@matches;
}


sub matchToAllele {
	my ($self,$databaseId,$locus,$seq)=@_;
	
	my %allele=(
		"database" => $databaseId,
		"locus"		 => $locus,
		"allele"	 => -1,
		"sequence" => '',
		"diffs"		 => 0
	);
	my $soapResponse = $self->{CONN}->locusQuery($databaseId,$locus,$seq);

	if ($soapResponse->fault) {
		return $self->handleError("matchToAllele", 7, $soapResponse->faultstring);
	} else {
		my @subs;
		for my $t ($soapResponse->valueof('//substitution')) {
			push @subs, {"position"=>$t->{'position'}, "allele_nuc"=>$t->{'nucleotide1'}, "seq_nuc"=>$t->{'nucleotide2'}};
		}
		$allele{"allele"} = $soapResponse->result();
		if (@subs) {
			$allele{"diffs"} = scalar @subs;
			$allele{"substitutions"} = \@subs;
		}
	}

	return [\%allele];
}


sub getSTs {
	my ($self, $databaseId, $profile) = @_;
	
	my @profileElements=();
	foreach my $locus (keys %$profile) {
		push @profileElements, SOAP::Data->name('alleleNumber' => \SOAP::Data->value(
			SOAP::Data->name('locus' => $locus),
			SOAP::Data->name('id' => $profile->{$locus})));
	}
	my $soapResponse = $self->{CONN}->getSTs($databaseId, @profileElements);
	
	my @sts=();
	unless ($soapResponse->fault) {
		for my $st ($soapResponse->valueof('//ST')) {
			my $prohash = $self->getProfile($databaseId, $st);
			push @sts, {$st=>$prohash};
		}
		return \@sts;
	} else {
		return $self->handleError("getSTs", 8, $soapResponse->faultstring);
	}

}


sub getProfile {
	my ($self, $databaseId, $st) = @_;

	my $soapResponse = $self->{CONN}->getProfile($databaseId, $st);
	unless ($soapResponse->fault) {
		my %profile=();
		for my $t ($soapResponse->valueof('//alleleNumber')) {
			$profile{$t->{'locus'}} = $t->{'id'};
		}
		return \%profile;
	} else {
		return $self->handleError("getProfile", 9, $soapResponse->faultstring);
	}
}


1;
