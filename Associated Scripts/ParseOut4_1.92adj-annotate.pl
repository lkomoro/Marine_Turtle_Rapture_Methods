#!/usr/bin/perl

$file = $ARGV[0];
$file2 = $ARGV[1];

open(FILE, "<$file") or die;

while (<FILE>) {
	$line = $_; chomp($line);

	@tabs = split(/\t/,$line);
	$scaff = $tabs[0]; $pos = $tabs[1]; $maj = $tabs[2]; $min = $tabs[3]; $maf = $tabs[4]; $nind = $tabs[5];
	if ($maf >= 0.1 && $maf <= 0.4) {
		$key = $scaff . "_" . $pos;
		$hash{$key} = "$pos,$maj,$min,$maf,$nind";
	} #this takes the pac.mafs.txt file and identifies all the places where there is a potential SNP (as we define of maf>0.1 and <0.4)
#memorizes the scaffold_SNPposition as the key here, as well as the supporting information in the hash line but not sure, need to look up more to understand that part
#but overall, in the end have a list of scaffold_SNPposition to use below
}
close FILE;


open(FILE, "<$file2") or die;

while (<FILE>) {
        $name = $_; chomp($name);
	$seq = <FILE>; chomp($seq);

	substr($name,0,1) = "";
	($scaff, $pos) = split(/\_/, $name);#defines scaff and pos as what is in the all.out4.txt where each sequence has a header that is "scaffoldX_pos"

	$stringl = "";
	$countl = 0;
	$x = $pos-84; # here 'pos' is now redefined, and is now the first position of the cut site, in reference to where it is on its respective scaffold
#so this means that "x" here is the start of the 'good coverage' area we want to look at on the left side pairs that straddle the cut site
	while ($x < $pos) { #so while the x is less than the start of the cut site (meaning that we have a left side sequence)
		$key = $scaff . "_" . $x; #is this from the key above to call in our SNP positions that meet the criteria of being in the region we want,
		# or is this now being redefined? Clarify
		if ($hash{$key} ne "") { # i think this is just saying do it except for when the hash key is not equal to nothing (?)
			$countl++; #'++ placed after a variable increments the variable by one after returning the value'
			#not sure about these lines -is looping through the memorized list from above for when it meets the criteria?
			$stringl = $stringl . $hash{$key} . ";";
		}
		$x++;
	}

	$stringr = "";
	$countr = 0;
        $x = $pos+8; #since pos is first position of the cut site, x is the start of the 'good coverage' for the right side of the pairs that straddle the cut site
        while ($x < ($pos+92)) { #so while the x is less than the end of the 'good coverage' section on the right side
		$key = $scaff . "_" . $x; #save the information for the scaffold ID and where the 'good section' starts for the right sequences
                if ($hash{$key} ne "") {#see above...
			$countr++;
			$stringr = $stringr . $hash{$key} . ";";
		}
                $x++;
        }

	chop($stringl);  chop($stringr);
#so outcome of this section are the lists of SNPs that fit into our regions of interest?

	$pass = "";
	if ($countl == 1 && $countr == 0) {#if there is one SNP on the left and none on the right, then...
		$seq2 = substr($seq, 500-120, 128);#then seq2 is a substring of the sequence starting at position 500 (start of cut site) minus 120, going for 128 bases
		#this creates sequence of 120 bp plus cut site 8base sequence at end
		$seq3 = reverse($seq2); $seq3 =~ tr/ACGTacgt/TGCAtgca/;#seq3 is the reverse complement of of seq2-which we want for the left side ones
		#note that now it starts with the 8base cut site instead of ending with it-these will get taken out before ordering baits
		print ">$name" . "_R " . "$stringl\n";#print the name and "R" (to designate it's a reverse complement sequence) and
		#the info we saved above with the SNP position,major,minor alleles, MAF, and #individuals (with a return to start a new line)
		print "$seq3\n";#print the 128 base sequence on the next line
		$pass = $stringl; #go to next one?
	} elsif ($countl == 0 && $countr == 1) {#if there is one SNP on the right and none on the left, then...
		$seq2 = substr($seq, 500, 128);#sequence starting at cut site, going for 128 bases to right
		print ">$name" . "_F " . "$stringr\n"; #print name...see above
                print "$seq2\n"; #print sequence...see above
		$pass = $stringr;#go to next one?
	} elsif ($countl == 1 && $countr == 1) {#if there is one SNP on either side, choose the right one, follows that of above
		$seq2 = substr($seq, 500, 128);
                print ">$name" . "_F " . "$stringr\n";
                print "$seq2\n";
		$pass = $stringr;
	}

	if ($pass ne "") { #if it doesn't meet any of the criteria above, skip it
		($site) = split(/,/, $pass);
		#print "$scaff\t$site\n"; #for when we want to print all the identifiers for the possible SNPs
	}


}
close FILE;
