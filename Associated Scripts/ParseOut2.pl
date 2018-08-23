#!/usr/bin/perl

$file = $ARGV[0];
$file2 = $ARGV[1];

open(FILE, "<$file") or die;

while (<FILE>) {
	$line = $_; chomp($line);
	
	@tabs = split(/\t/,$line);
	$scaff = $tabs[0];
	$pos = $tabs[3] + 1;

	push(@{$hash{$scaff}}, $pos);

}
close FILE;

#print scalar(@{$hash{scaffold28}});

open(FILE, "<$file2") or die;

while (<FILE>) {
        $name = $_; chomp($line);
	$seq = <FILE>; chomp($seq);

	substr($name,0,1) = "";
	($name) = split(/\s/, $name);
#	print "$name\n";

	if ($hash{$name} ne "") {
		$x = 0;
		while ($x < scalar(@{$hash{$name}})) {
			$pos = $hash{$name}[$x];
			$start = $pos - 501;
			$reg = substr($seq, $start, 1008);
			$site = substr($seq, ($pos-1), 8);

			if ($site eq "CCTGCAGG" && length($reg) == 1008) {
				print ">" . $name . "_" . $pos . "\n" . $reg . "\n";
			} else {
				#print "$site\t" . length($reg) . "\n";
			}
			$x++; 
		}
	}

}
close FILE;



