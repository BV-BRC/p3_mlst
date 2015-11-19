package VASTR::Tools;

use strict;


sub rc {
	my $s=shift;
	chomp $s;
	my $revseq = reverse $s;
	$revseq =~ tr/actgACTG/tgacTGAC/;
	return $revseq;
}



1;
