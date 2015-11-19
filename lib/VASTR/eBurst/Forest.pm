package VASTR::eBurst::Complex;

use strict;
use warnings;

use Data::Dumper;

sub new {
	my ($class,$a)=@_;
	
	my $self = {
		ALIST         => [],
		IX2ID         => [],
		NODES         => {},
		EDGES         => {},
		ALLOW_DUP_IDS => 1
	};
	$self->{ALLOW_DUP_IDS}=$a if defined $a;
	
	bless($self,$class);
	return $self;
}


sub addNode {
	my ($self,$nid,$size)=@_;
	$size=1 unless defined $size;
	
	my $index;
	if (!$self->{ALLOW_DUP_IDS} and $self->hasNode($nid)) {
		# don't add a new node if already present and dups not allowed
		$index=$self->getIndex($nid)->[0];
	} else {
		$index=scalar @{$self->{ALIST}};
		$self->{NODES}->{$nid} = [] unless $self->hasNode($nid);
		push @{$self->{NODES}->{$nid}}, {'index'=>$index, 'size'=>$size};
		$self->{ALIST}->[$index]=[];
		$self->{IX2ID}->[$index]=$nid;
	}

	return $index;
		
}


sub hasNode {
	my ($self, $nid)=@_;
	my $count=0;
	$count=scalar(@{$self->{NODES}->{$nid}}) if defined $self->{NODES}->{$nid};
	return $count;
}


sub getNode {
	my ($self, $nid)=@_;
	if ($self->hasNode($nid)) {
		return $self->{NODES}->{$nid};
	} else {
		return [];
	}
}


sub getNodes {
	my $self=shift;
	my @nids = keys %{$self->{NODES}};
	return \@nids;
}


sub setNodeProp {
	my ($self, $nid, $propHash)=@_;
	foreach my $node (@{$self->getNode($nid)}) {
		foreach my $k (%$propHash) {
			$node->{$k} = $propHash->{$k};
		}
	}
}


sub getParent {
	my ($self, $nid)=@_;

	return undef unless $self->hasNode($nid);
	
	my @pnids   =();
	my $numArr  =$self->getIndex($nid);
	my %numHash =map {$_ => 1} @$numArr;
	
	for (my $snum=0; $snum<scalar(@{$self->{ALIST}}); $snum++) {
		foreach my $tnum (@{$self->{ALIST}->[$snum]}) {
			push @pnids, $self->{IX2ID}->[$snum] if defined $numHash{$tnum};
		}
	}
	
	return \@pnids;
}



sub addEdge {
	my ($self, $sid, $tid, $len)=@_;
	$len=1 unless defined $len;
	
	return if $self->hasEdge($sid, $tid);
	
	my $snums = $self->getIndex($sid);
	my $tnums = $self->getIndex($tid);
	
	foreach my $snum (@$snums) {
		foreach my $tnum (@$tnums) {

			push @{$self->{ALIST}->[$snum]}, $tnum;
			$self->{EDGES}->{$snum}={} unless defined $self->{EDGES}->{$snum};
			$self->{EDGES}->{$snum}->{$tnum} = {
				'source'=> $sid,
				'target'=> $tid,
				'uid'   => "$snum-$tnum",
				'length'=> $len
			};

			push @{$self->{ALIST}->[$tnum]}, $snum;
			$self->{EDGES}->{$tnum}={} unless defined $self->{EDGES}->{$tnum};
			$self->{EDGES}->{$tnum}->{$snum} = {
				'source'=> $tid,
				'target'=> $sid,
				'uid'   => "$tnum-$snum",
				'length'=> $len
			};
		}
	}
}


sub hasEdge {
	my ($self, $sid, $tid)=@_;

	my $snums = $self->getIndex($sid);
	my $tnums = $self->getIndex($tid);
	
	foreach my $snum (@$snums) {
		foreach my $tnum (@$tnums) {
			foreach my $e (@{$self->{ALIST}->[$snum]}) {
				return 1 if $e == $tnum;
			}
		}
	}
	return 0;
}


sub getEdge {
	my ($self, $sid, $tid)=@_;
	
	my @edges=();
	
	my $snums = $self->getIndex($sid);
	my $tnums = $self->getIndex($tid);
	my %thash = map {$_ => 1} @$tnums;
	
	foreach my $snum (@$snums) {
		for (my $i=0; $i<scalar(@{$self->{ALIST}->[$snum]}); $i++) {
			my $tnum=$self->{ALIST}->[$snum]->[$i];
			if (defined $thash{$tnum}) {
				if (defined $self->{EDGES}->{$snum} and defined $self->{EDGES}->{$snum}->{$tnum}) {
					push(@edges, $self->{EDGES}->{$snum}->{$tnum});
				}
				if (defined $self->{EDGES}->{$tnum} and defined $self->{EDGES}->{$tnum}->{$snum}) {
					push(@edges, $self->{EDGES}->{$tnum}->{$snum});
				}
			}
		}
	}
	return \@edges;
}


sub setEdgeProp {
	my ($self, $sid, $tid, $propHash)=@_;
	foreach my $edge (@{$self->getEdge($sid, $tid)}) {
		foreach my $k (keys %$propHash) {
			$edge->{$k} = $propHash->{$k};
		}
	}
}

sub getEdgeProp {
	my ($self, $sid, $tid, $key)=@_;
	
	my %props=();
	foreach my $edge (@{$self->getEdge($sid, $tid)}) {
		$props{$edge->{'uid'}}=$edge->{$key};
	}
	return \%props;
}

sub _getEdgePropByIndex {
	my ($self, $snum, $tnum, $key)=@_;
	
	my $val=0;
	if (defined $self->{EDGES}->{$snum} and defined $self->{EDGES}->{$snum}->{$tnum} and defined $self->{EDGES}->{$snum}->{$tnum}->{$key}) {
		$val=$self->{EDGES}->{$snum}->{$tnum}->{$key};
	} elsif (defined $self->{EDGES}->{$tnum} and defined $self->{EDGES}->{$tnum}->{$snum} and defined $self->{EDGES}->{$tnum}->{$snum}->{$key}) {
		$val=$self->{EDGES}->{$tnum}->{$snum}->{$key};
	}
	return $val;
}


sub BFS {
	my ($self, $sid, $tid)=@_;
	my %discovered=();

	$discovered{$sid}=1;
		
	my @list=([]);
	push @{$list[0]}, $sid;
	my $i=0;
	
	while (scalar @{$list[$i]} > 0) {

		$list[$i+1]=[];

		for (my $j=0; $j<scalar(@{$list[$i]}); $j++) {
			my $u=$list[$i]->[$j];

			foreach my $v (@{$self->{ALIST}->[$u]}) {
				unless (defined $discovered{$v}) {
					$discovered{$v}=1;
					push @{$list[$i+1]}, $v;
				}
			}
		}
		$i++;
	}
	
	return \@list;
}


sub getPath {
	my ($self, $snum, $tnum, $discovered)=@_;
	$discovered={} unless defined $discovered;
	$discovered->{$snum} = 1;

	my @path=($snum);
	foreach my $child (@{$self->{ALIST}->[$snum]}) {
		unless (defined $discovered->{$child}) {
			$discovered->{$child}=1;
			if ($child==$tnum) {
				push @path, $child;
				last;
			} else {
				my $subpath=$self->getPath($child, $tnum, $discovered);
				my $found=0;
				foreach my $subnum (@$subpath) {
					$discovered->{$subnum}=1;
					$found=1 if $subnum == $tnum;
				}
				if ($found) {
					@path = (@path, @$subpath);
					last;
				}
			}
		}
	}
	
#print STDERR "$snum -> $tnum\n";
#print STDERR join (",", @path) . "\n";

	return \@path;

}


sub getDistance {
	my ($self, $sid, $tid)=@_;
	
	my @d=();
	
	my $snums = $self->getIndex($sid);
	my $tnums = $self->getIndex($tid);
	my %thash = map {$_ => 1} @$tnums;
	foreach my $snum (@$snums) {
		foreach my $tnum (@$tnums) {
			my $pathd=0;
			my $path=$self->getPath($snum, $tnum);
#print STDERR Dumper($path);
			my $six;
			my $tix= shift @$path;
			while (scalar @$path > 0) {
				$six= $tix;
				$tix= shift @$path;
				my $edged = $self->_getEdgePropByIndex($six, $tix, 'length');
#print STDERR "$s -> $t : $edged\n";
				$pathd += $self->_getEdgePropByIndex($six, $tix, 'length');
			}
			push @d, $pathd;
		}
	}
	
	return \@d;
}


sub printAdjMatrix {
	my $self=shift;
	
	# print target ids across the top
	print "\t" . join("\t", @{$self->{IX2ID}}) . "\n";
	
	for (my $i=0; $i<scalar(@{$self->{IX2ID}}); $i++) {
		my $sid=$self->{IX2ID}->[$i];
		print "$sid";
		for (my $j=0; $j<scalar(@{$self->{IX2ID}}); $j++) {
			my $tid=$self->{IX2ID}->[$j];
			my $edge=0;
			$edge=$self->_getEdgePropByIndex($i, $j, 'length') if $self->hasEdge($sid, $tid);
			print "\t$edge";
		}
		print "\n";
	}
}



sub getIndex {
	my ($self, $nid)=@_;
	my @a=();
	for (my $i=0; $i<scalar(@{$self->{IX2ID}}); $i++) {
		push @a, $i if "$self->{IX2ID}->[$i]" eq "$nid";
	}
	return \@a;
}



1;
