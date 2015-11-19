package VASTR::MLST::MLST;

use strict;
use warnings;

use File::HomeDir qw/home/;

sub new {
	my $class=shift;
	
	my $tmpdir=shift;
	unless (defined $tmpdir && -d $tmpdir) {
		my $home=File::HomeDir->my_home;
		$tmpdir="$home/mlst-tmp";
	}
	
	my $self = {
		CONN    => undef,
		TMPDIR  => $tmpdir
	};
	
	bless($self,$class);
	return $self;
}


=cut
defined in sub-modules
=cut
sub getAvailableDatabases {
	return 0;
}

=cut
sub getDatabasesByGenome2 {
	my ($self,$genomeName) = @_;
	
	unless (defined $self->{DBS}) {
		my $dbs = $self->getAvailableDatabases;
		return $dbs if (defined $dbs->{"errorCode"});
	}
	
	my $dbsToQuery=[];
	foreach my $dbName (keys %{$self->{DBS}}) {
		my @descriptions = ($self->{DBS}->{$dbName});
		if ($descriptions[0] =~ /^(\w+\s\w+)\sand\s(\w+\s\w+)$/) {
			$descriptions[0]=$1;
			$descriptions[1]=$2;
		} elsif ($descriptions[0] =~ /^(\w+)\sspp\./ ) {
			$descriptions[0]=$1;
		} else {
			my @a = split /\s/, $descriptions[0], 3;
			pop @a unless scalar @a < 3;
			$descriptions[0]=join(" ", @a);
		}
		for (my $i=0; $i<scalar(@descriptions); $i++) {
			my @g = split /\s/, $genomeName;
			my @d = split /\s/, $descriptions[$i];
			if (scalar @d == 1 && (lc $descriptions[$i] eq lc $g[0])) {
				push(@$dbsToQuery, $dbName);
			} elsif (scalar @d == 2 && (lc $descriptions[$i] eq lc($g[0])." ".lc($g[1]))) {
				push(@$dbsToQuery, $dbName);
			}
		}
	}
	return $dbsToQuery;
}
=cut

sub matchDatabases {
	my ($self,$inputString) = @_;
	
	my $dbsToQuery=[];

	unless (defined $self->{DBS}) {
		my $dbs = $self->getAvailableDatabases;
		return $dbs if (defined $dbs->{"errorCode"});
	}
	
	my @searchStringElements = split /\s/, $inputString;

	# limit input to 'genus species'
	while (scalar @searchStringElements > 2) {
		pop @searchStringElements;
	}
	my $searchString = join " ", @searchStringElements;
	$searchString =~ s/\s/\\s/gi;
	
	# first look for a stringent 'genus species' match
	#warn "$searchString\n";
	foreach my $dbName (keys %{$self->{DBS}}) {
		my $description = $self->{DBS}->{$dbName};
		if ($description =~ /$searchString/i || $dbName =~ /$searchString/i) {
			push(@$dbsToQuery, $dbName);
			#warn "  $dbName\n";
		}
	}
	
	# if no 'genus species' match, look for genus alone
	if (scalar @$dbsToQuery == 0) {

		# limit search to dbs of format 'GENUS_spp'
		foreach my $dbName (keys %{$self->{DBS}}) {
			push(@$dbsToQuery, $dbName) if $dbName =~ /(.+?)_spp/ && lc $searchStringElements[0] eq lc $1;
		}
		
	}
	return $dbsToQuery;
}


=cut
defined in sub-modules
=cut
sub getLoci {
	return 0;
}


=cut
defined in sub-modules
=cut
sub getAllele {
	return 0;
}


=cut
defined in sub-modules
=cut
sub findMatches {
	return 0;
}


=cut
defined in sub-modules
=cut
sub matchToAllele {
	return 0;
}


=cut
defined in sub-modules
=cut
sub getSTs {
	return 0;
}


=cut
defined in sub-modules
=cut
sub getProfile {
	return 0;
}



sub handleError {
	my ($self, $errSource, $errorCode, $errStr) = @_;
	
	my $errObj={"errorCode"=>$errorCode, "thrownBy"=>$errSource, "description"=>$errStr};
	return $errObj;
}


1;

