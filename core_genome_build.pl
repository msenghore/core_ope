#!/usr/bin/perl

$fin = shift;
$ref = shift;
$relax = shift;

$edge = 3;
$intergap = 100;

open (FIN, "<$fin") or die;
while (<FIN>) {
	if($_ =~ /^>(\S+)/) {
		$name[@name] = $1;
	} else {
		chomp $_;
		$seq{$name[@name-1]} .= $_;
	}
}
close FIN;
if($ref eq "") {$ref = $name[0];}

local %gaps;
for(local $j=0; $j<@name; ++$j) {
	local @blank; local %g;
	for(local $m=0; $m<length($seq{$name[0]}); $m+=1000) {
		if(substr($seq{$name[$j]}, $m, 1000) =~ /[\.-]/) {
			for(local $i=$m; $i<length($seq{$name[0]}) && $i < $m + 1000; ++$i) {
				if(substr($seq{$name[$j]}, $i, 1) eq "-") {
					if(@blank > 0 && $i == $blank[@blank-1][1]+1) {
						$blank[@blank-1][1] = $i;
					} else {
						@{$blank[@blank]} = ($i, $i);
					}
				} elsif (substr($seq{$name[$j]}, $i, 1) eq ".") {
					$g{$i} ++;
					substr($seq{$name[$j]}, $i, 1, "-");
				} elsif (substr($seq{$name[$j]}, $i, 1) eq "n" || substr($seq{$name[$j]}, $i, 1) eq "N") {
					substr($seq{$name[$j]}, $i, 1, "-");
				}
			}
		}
	}
	for(local $i=0; $i<@blank; ++$i) {
		for(local $m = $blank[$i][0] - $edge; $m<=$blank[$i][1] + $edge; ++$m) {
			$g{$m} ++;
		}
		if($blank[$i][1] - $blank[$i][0] < 10) {
			splice(@blank, $i, 1);
			$i --;
		} else {
			if($blank[$i][0]-$edge >= 0) {
				substr($seq{$name[$j]}, $blank[$i][0]-$edge, $edge, "-" x $edge);
			}
			if($blank[$i][1]+$edge < length($seq{$name[$j]})) {
				substr($seq{$name[$j]}, $blank[$i][1]+1, $edge, "-" x $edge);
			}
		}
	}
	
	for(local $i=1; $i<@blank; ++$i) {
		if($blank[$i][0] - $blank[$i-1][1] <= $intergap) {
			for(local $m = $blank[$i-1][1]+1; $m < $blank[$i][0]; ++$m) {
				$g{$m} ++;
			}
			substr($seq{$name[$j]}, $blank[$i-1][1]+1, $blank[$i][0]-$blank[$i-1][1]-1, "-" x ($blank[$i][0]-$blank[$i-1][1]-1));
		}
	}
	foreach $site (keys %g) {
		$gaps{$site} ++;
	}
}

local @region;
for(local $i=0; $i<length($seq{$name[0]}); ++$i) {
	if($gaps{$i} > $relax) {
		if(@region > 0 && $region[@region-1][1] == $i-1) {
			$region[@region-1][1] = $i;
		} else {
			@{$region[@region]} = ($i, $i);
		}
	}
}

print ">$ref\n$seq{$ref}\n";

for(local $j=0; $j<@name; ++$j) {
	if($name[$j] eq $ref) {
		next;
	}
	for(local $i=0; $i<@region; ++$i) {
		substr($seq{$name[$j]}, $region[$i][0], $region[$i][1]-$region[$i][0]+1, "-" x ($region[$i][1]-$region[$i][0]+1));
	}
	print ">$name[$j]\n$seq{$name[$j]}\n";
}
