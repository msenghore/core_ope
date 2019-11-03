#! /usr/bin/perl

$reference = shift;
$query = shift;
$out = shift;

`nucmer --mum -p $query $reference $query`;
`delta-filter -1 $query.delta > $query.delta.2`;
`perl /usr/local/bioinf/core_ope/strict_filter.pl $query.delta.2 > $query.delta.3`;
`perl /usr/local/bioinf/core_ope/aln_extract.pl $query.delta.3 > $out`;
`rm $query.delta $query.delta.?`;
