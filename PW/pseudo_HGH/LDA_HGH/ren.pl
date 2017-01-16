#!/usr/local/bin/perl
#

open (DATA,"atomic.data") || die("File not found...\n");
while(<DATA>) {
	chop;
	($n,$o) = split /\t/;
	$name{$n} = $o;
}
close DATA;
#
#
opendir(DIR,"./");
@allfiles = grep /hgh$/, readdir DIR;
closedir DIR;
#
#
foreach $k1 (@allfiles) {
	($k = $k1) =~ tr/A-Z/a-z/;
	($n,$r) = split /\./,$k;
	$r =~ s/hgh/sc/;
	$n1 = @name{$n}.$n."\.".$r."\.hgh";
	link $k1,$n1;
#	print "$k1  $n1\n";
}
