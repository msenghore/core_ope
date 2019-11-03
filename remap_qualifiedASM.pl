#! /usr/bin/perl

$fin = shift;
$f1 = shift;
$f2 = shift;
$fout = shift;
$times = shift;

if($times eq "") {
	$times = 999;
}
$bbuild = "bowtie2-build";
$bowtie2 = "bowtie2";
$samtools = "samtools";
$bcftools = "bcftools";

$fin2 = $fin;
local $prev_pair = 0;
for(local $ite=0; $ite<$times; ++$ite) {
	if($ite == $times-1 && $fin =~ /Scaffold.fa/) {
		next;
	}
	local @name; local %seq;
	open (FIN, "<$fin2") or die;
	while (<FIN>) {
		if($_ =~ /^>(\S+)/) {
			$name[@name] = $1;
		} else {
			chomp $_;
			$seq{$name[@name-1]} .= $_;
		}
	}
	close FIN;
	`$bbuild $fin2 $fin2`;
	`$bowtie2 -p 2 -x $fin2 -1 $f1 -2 $f2 -S $fin.cor.sam`;
	
	`rm -f $fin2.*.bt2`;
	`$samtools faidx $fin2`;
	`$samtools view -bt $fin2.fai -o $fin.cor.bam $fin.cor.sam`;
	`rm -f $fin.cor.sam $fin2.fai`;
	`$samtools sort $fin.cor.bam $fin.cor.sort`;
	`rm -f $fin.cor.bam`;
	# `$samtools mpileup -C50 -gf $fin2 $fin.cor.sort.bam > $fin.cor.bcf`;
	# `rm -f $fin.cor.sort.bam`;

	if($ite == $times-1) {
		`$samtools mpileup -I -C50 -gf $fin2 $fin.cor.sort.bam > $fin.cor.bcf`;
		foreach $n (@name) {
			$seq{$n} =~ s/[a-zA-Z]/N/g;
		}
		open (SNP, "$bcftools view -cp 1 $fin.cor.bcf |") or die;
		while (<SNP>) {
			chomp $_;
			if (substr($_, 0, 1) eq "#") {next;}
			local @part = split /\t/, $_;
			if (substr($part[7], 0, 5) eq "INDEL") {
				next;
			}
			if ($part[7] =~ /AF1=([e\d\.]+)/) {
				if ($1 >= 0.2 && $1 <= 0.8) {
	                        	next;
	                	}
	        	}
		        if($part[7] =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) {
	        	        if($1 + $2 + $3 + $4 < 3) {
	                	        next;
		                } elsif (($1 + $2) > 0.8 * ($1 + $2 + $3 + $4) && $1 > 0.7 * ($1+$3) && $2 > 0.7 * ($2+$4)) {
	        	                substr($seq{$part[0]}, $part[1]-1, 1, $part[3]);
	                	} elsif (($3 + $4) > 0.8 * ($1 + $2 + $3 + $4) && $3 > 0.7 * ($1+$3) && $4 > 0.7 * ($2+$4)) {
                	        	substr($seq{$part[0]}, $part[1]-1, 1, substr($part[4],0,1));
		                }
	        	}

		}
		
		open (FOUT, ">$fout") or die;
		foreach $n (@name) {
			print FOUT ">$n\n";
			for(local $i=0; $i<length($seq{$n}); $i+=100) {
				print FOUT substr($seq{$n}, $i, 100)."\n";
			}
		}
		close FOUT;
	} else {
		`$samtools mpileup -C50 -gf $fin2 $fin.cor.sort.bam > $fin.cor.bcf`;
		open (SNP, "$bcftools view -vp 1 $fin.cor.bcf |") or die;
		local @indel; local $change = 0;
		while (<SNP>) {
			chomp $_;
			local @part = split /\t/, $_;
			
			if($part[7] =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) {
				if($1 + $2 < $3 + $4) {
					if(substr($part[7], 0, 5) ne "INDEL") {
						substr($seq{$part[0]}, $part[1]-1, 1, substr($part[4],0,1));
					} else {
						@alter = split /,/, $part[4];
						@{$indel[@indel]} = ($part[0], $part[1], $part[3], $alter[0], $1+$2, $3+$4);
					}
				}
			}
		}
		for(local $i=0; $i<@indel; ++$i) {
			local $end = $indel[$i][1] + length($indel[$i][2]) -1;
			for(local $j=$i+1; $j<@indel; ++$j) {
				if($end >= $indel[$i][1]) {
					if($indel[$j][3] =~ /N/) {
						splice (@indel, $j, 1); $j --; next;
					} elsif ($indel[$i][3] =~ /N/) {
						splice(@indel, $i, 1); $i --; last;
					} elsif($indel[$j][5]-$indel[$j][4] > $indel[$i][5]-$indel[$i][4]) {
						splice(@indel, $i, 1); $i --; last;
					} else {
						splice (@indel, $j, 1); $j --; next;
					}
				}
			}
		}
		
		for(local $i= @indel-1; $i>=0; $i --) {
			substr($seq{$indel[$i][0]}, $indel[$i][1]-1, length($indel[$i][2]), $indel[$i][3]);
		}
		
		open (FOUT, ">$fin.$ite.fasta") or die;
		foreach $n (@name) {
			print FOUT ">$n\n";
			for(local $i=0; $i<length($seq{$n}); $i+=100) {
				print FOUT substr($seq{$n}, $i, 100)."\n";
			}
		}
		close FOUT;
		`rm -f $fin2.fai`;
		if($fin ne $fin2) {
			`rm -f $fin2`;
		}
		$fin2 = "$fin.$ite.fasta";
		
		if (@indel < 1) {
			$ite = $times -2;
		}
	}
	# `rm -f $fin.cor.bcf`;
}
