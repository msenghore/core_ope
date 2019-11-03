#! /usr/bin/perl

$list = shift;
open (LIST, "<$list") or die;
local $head = <LIST>; chomp $head;
@names = split /\t/, $head;
while (<LIST>) {
	chomp $_;
	local @part = split /\t/, $_;
	for(local $i=1; $i<@names; ++$i) {
		$seq[$i] .= $part[$i];
	}
}
close LIST;

for(local $i=1; $i<@names; ++$i) {
	$names[$i] =~ s/:\d+-\d+$//;
	print ">$names[$i]\n$seq[$i]\n";
}
