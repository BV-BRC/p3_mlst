#! /usr/bin/perl -w

use strict;
use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage .=   "Usage: $progname [options]\n";
$usage .=   "Calculates an eBurst forest from the MLST data on STDIN and writes \n";
$usage .=   "adjacency lists to STDOUT.\n";
$usage .=   "       [-d D] use D as threshold distance for including edges (1).\n";
$usage .=   "       [-m M] also print the distance matrix to file M (none).\n";
$usage .=   "       [-s]   include singleton nodes in output, connected to nearest \n";
$usage .=   "              neighbor (false). singletons ignore the threshold setting.\n";
$usage .=   "\n";

my $matrixf;
my $threshold =1;
my $singletons=0;

while (@ARGV) {
	my $arg=shift;
	if ($arg eq '-h' or $arg eq '-help') {
		die "$usage";
	} elsif ($arg eq '-d' or $arg eq '-dist' or $arg eq '-distance') {
		defined ($threshold=shift) or die "FATAL: -d flag found but no threshold provided.\n$usage";
	} elsif ($arg eq '-m' or $arg eq '-matrix') {
		defined ($matrixf=shift) or die "FATAL: -m flag found but no output file name provided.\n$usage";
	} elsif ($arg eq '-s' or $arg eq '-singles' or $arg eq '-singletons') {
		$singletons=1;
	}
}

my %profiles=();
my %dm      =();
my $loci    =0;

# read in the MLST profile data
my $rmCC=-1;
while (<>) {
	chomp;
	my @cols=split /\t/, "$_", -1;

	if ($rmCC==-1) {
		if (lc $cols[-1] eq 'clonal_complex' or lc $cols[-1] eq 'mlst_clade' or lc $cols[-1] eq 'species') {
			$rmCC=1;
		} else {
			$rmCC=0;
		}
		next;
	}
	pop @cols if $rmCC;
	
	my $st=shift @cols;
	$profiles{$st}=\@cols;
	$dm{$st}={};
	
	$loci=scalar(@cols) unless $loci;
}


my @alists=();
my @stIds = sort {$a<=>$b} keys %profiles;

# init the adjacency lists. index 0 holds SLV,DLV,TLV totals
for (my $i=0; $i<scalar(@stIds); $i++) {
	$alists[$stIds[$i]]=[[0,0,0]];
	for (my $j=1; $j<=$loci; $j++) {
		push @{$alists[$stIds[$i]]}, [];
	}
}

# calculate the distances
for (my $i=0; $i<scalar(@stIds); $i++) {
	my $st1=$stIds[$i];

	for (my $j=$i+1; $j<scalar(@stIds); $j++) {
		my $st2=$stIds[$j];

		my $d=0;
		for (my $k=0; $k<scalar(@{$profiles{$st1}}); $k++) {
			$d++ unless $profiles{$st1}->[$k] == $profiles{$st2}->[$k];
			#last if $d>$threshold;
		}
		$dm{$st1}->{$st2}=$d if defined $matrixf;
		
		next if $d==0;
		push @{$alists[$st1]->[$d]}, $st2;
		#push @{$alists[$st2]->[$d]}, $st1;
		if ($d<4) {
			$alists[$st1]->[0]->[$d-1]=$alists[$st1]->[0]->[$d-1]+1;
			$alists[$st2]->[0]->[$d-1]=$alists[$st2]->[0]->[$d-1]+1;
		}
	}
}


# print the adjacency lists
for (my $i=0; $i<scalar(@alists); $i++) {
	next unless defined $alists[$i];
	print "$i";
	for (my $j=0; $j<=$threshold; $j++) {
		print "\t" . join(",", @{$alists[$i]->[$j]});
	}
	if ($singletons and $alists[$i]->[0]->[0] == 0) {
		# assumes SLV graph only; refactor to use threshold instead
		for (my $j=$threshold+1; $j<scalar(@{$alists[$i]}-1); $j++) {
			print "\t" . join(",", @{$alists[$i]->[$j]});
			last if scalar @{$alists[$i]->[$j]} > 0;
		}
	}
	print "\n";
}


# print the full distance matrix if appropriate
if (defined $matrixf) {
	open my $dmfh, ">", "$matrixf" or warn "Unable to print distance matrix: $!";
	if (defined $dmfh) {
		#print $dmfh join("\t", @stIds) . "\n";
		for (my $i=0; $i<scalar(@stIds); $i++) {
			my $st1=$stIds[$i];
			#print $dmfh "$st1";
			my @arr=();
			for (my $j=$i+1; $j<scalar(@stIds); $j++) {
				my $st2=$stIds[$j];
				#print $dmfh "\t";
				if (defined $dm{$st1}->{$st2}) {
					push @arr, "$dm{$st1}->{$st2}";
				} elsif (defined $dm{$st2}->{$st1}) {
					push @arr, "$dm{$st2}->{$st1}";
				} else {
					push @arr, "0";
				}
			}
			print $dmfh join("\t", @arr) . "\n" unless scalar @arr == 0;
		}
		close $dmfh;
	}
}


exit;


=cut
{
		var T={};	// hash of edge ids
		
		var set2node={};	// maps set id to hash of node ids
		var node2set={};	// maps node id to its set id
		
		var E=this.networkModel.data.edges;
		E.sort(function(a,b) {
			return (a.weight < b.weight) ? 1 : ((b.weight < a.weight) ? -1 : 0);
		});
		
		var setIncr=1;
		for (var i=0; i<E.length; i++) {
			var u=E[i].source;
			var v=E[i].target;
			
			if (!node2set.hasOwnProperty(u)) {
				// u has not been examined yet
				if (node2set.hasOwnProperty(v)) {
					// v has been examined, so add u to v's set
					node2set[u]=node2set[v];
				} else {
					// neither u nor v have been examined, add them both to the next empty set
					node2set[u]=setIncr;
					node2set[v]=setIncr;
					setIncr++;
				}
			}
			
			if (!node2set.hasOwnProperty(v)) {
				// v has not been examined yet
				if (node2set.hasOwnProperty(u)) {
					// u has been examined, so add v to u's set
					node2set[v]=node2set[u];
				} else {
					// neither v nor u have been examined, add them both to the next empty set
					node2set[u]=setIncr;
					node2set[v]=setIncr;
					setIncr++;
				}
			}

			if (node2set[u] != node2set[v]) {
				// u and v are in different sets, so add the (u,v) to T
				T[E[i]]=1;
				// merge u's set and v's set into one
				var uSet=node2set[u];
				var vSet=node2set[v];
				for (var nodeId in set2node[vSet]) {
					if (set2node[vSet].hasOwnProperty(nodeId)) {
						node2set[nodeId]=uSet;
					}
				}
				delete set2node[vSet];
				
			}
		}
		
	}
=cut


