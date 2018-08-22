#!/usr/bin/perl

$file = $ARGV[0];

open(FILE, "<$file") or die;

while (<FILE>) {
        $name = $_; chomp($line);
	$seq = <FILE>; chomp($seq);

	substr($name,0,1) = "";
	($name) = split(/\s/, $name);

			$reg = substr($seq, 500, 128);
				print ">" . $name . $pos . "\n" . $reg . "\n";

			$x++;
}
close FILE;
