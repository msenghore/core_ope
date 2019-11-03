#! /usr/bin/perl

$fin = shift;
$region = shift;
$edge = shift;
if($edge eq "") {$edge = 10;}

open (REG, "<$region") or die;
while (<REG>) {
	@part = split /\s+/, $_;
	@{$repeat[@repeat]} = ($part[1], $part[2]);
}
close REG;

open (FIN, "<$fin") or die;
while (<FIN>) {
	if ($_ =~ /^>(\S+)/) {
		$name = $1;
		$names[@names] = $name;
	} else {
		chomp $_;
		$seq{$name} .= $_;
	}
}

print ">$names[0]\n$seq{$names[0]}\n";

for(local $i=1; $i<@names; ++$i) {
	foreach $r (@repeat) {
		for(local $j= $$r[0]-$edge-1 > 0 ? $$r[0]-$edge-1 : 0; $j < $$r[1]+$edge && $j < length($seq{$names[$i]}); ++$j) {
			substr($seq{$names[$i]}, $j, 1, "-");
		}
	}
	print ">$names[$i]\n$seq{$names[$i]}\n";
}
