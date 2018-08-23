#!/usr/bin/perl

$file = $ARGV[0];

open(FILE, "<$file") or die;

while (<FILE>) {
        $name = $_; chomp($line);
	$seq = <FILE>; chomp($seq);

	substr($name,0,1) = "";
	($name) = split(/\s/, $name);

			$reg = substr($seq, 380, 128);
$seq3 = reverse($reg); $seq3 =~ tr/ACGTacgt/TGCAtgca/;
				print ">" . $name . $pos . "\n" . $seq3 . "\n";

			$x++;
}
close FILE;
