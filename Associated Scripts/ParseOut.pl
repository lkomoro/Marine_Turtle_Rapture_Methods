#!/usr/bin/perl

$file = $ARGV[0];

open(FILE, "<$file") or die;

$line1 = "";
$line2 = "";

while (<FILE>) {
	$line = $_; chomp($line);
	if ($line1 eq "") {
		$line1 = $line;
	} elsif ($line2 eq "") {
		$line2 = $line;
	}
	
	if ($line1 ne "" && $line2 ne "") {
		@tabs1 = split(/\t/,$line1);
		@tabs2 = split(/\t/,$line2);
		if (($tabs1[3] + 9) == $tabs2[1] && $tabs1[0] eq $tabs2[0]) {
			print "$line1\t$line2\n";
		}
		$line1 = $line2;
		$line2 = "";
	}

}
close FILE;


