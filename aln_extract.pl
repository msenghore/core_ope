#! /usr/bin/perl

$delta = shift;
$snp = shift;

if($snp ne "") {
	open (SNP, "<$snp") or die; $n =0;
	while (<SNP>) {
		$n ++;
		@part =split /\t/, $_;
		$snp_list{"$part[17]:$part[1]"} = $part[3];
	}
	if($n == 0) {
		exit;
	}
}

open (DELTA, "<$delta") or die;

while (<DELTA>) {
	if($_ =~ /^>(\S+) (\S+) (\d+) (\d+)/) {
		chomp $_;
		$ref = $1; $qry = $2; $l_r = $3;
		
		if($seq{$ref} eq "") {
			for(local $i=0; $i<$l_r; ++$i) {
				$seq{$ref} .= "-";
			}
		}
		
		@align_file = `show-aligns -q $delta \"$ref\" \"$qry\"`;
		for(local $i=0; $i< @align_file; ++$i) {
			$line = $align_file[$i];
			if($line =~ /^-- BEGIN alignment \[ (..) (\d+) - (\d+) \| (..) (\d+) - (\d+) \]/) {
				$ref_seq = ""; $qry_seq = ""; $ref_d = $1; $ref_s = $2; $ref_e = $3; $qry_d = $4; $qry_s = $5; $qry_e = $6;
				for(;$i < @align_file; ++$i) {
					$i+=3;
					chomp $align_file[$i];
					if($align_file[$i] =~ /^$/) {
						last;
					}
					$ref_line = $align_file[$i]; $ref_line =~ /(\d+)\s+(\S+)/;
					$ref_seq .= $2;
					$qry_line = $align_file[$i+1]; $qry_line =~ /(\d+)\s+(\S+)/;
					$qry_seq .= $2;
				}
				local @s;
				local $ref_i = $ref_s - $ref_d; local $qry_i = $qry_s - $qry_d;
				for(local $i=0; $i<length($ref_seq); ++$i) {
					local $r_base = substr($ref_seq, $i, 1);
					local $q_base = substr($qry_seq, $i, 1);
					if($r_base ne ".") {
						$ref_i += $ref_d;
					}
					if($q_base ne ".") {
						$qry_i += $qry_d;
					}
					if($q_base eq "." || $r_base eq ".") {
					@{$s[@s]} = ($ref, $qry, $ref_i, $qry_i, $qry_d, $r_base, $q_base, $i);
					}
				}
				local @div_reg;
				for(local $m=0; $m<@s; ++$m) {
					local $tot_dis = 0;
					local $prev_r = $s[$m][2]; local $prev_q = $s[$m][3];
					for(local $n=$m+1; $n<@s; ++$n) {
						if($n > $m + 3) {
							last;
						}
						if($s[$n][0] eq $s[$m][0] && $s[$n][1] eq $s[$m][1]) {
							$tot_dis += ($s[$n][2]-$prev_r < $s[$n][4]*($s[$n][3] - $prev_q)) ? ($s[$n][2]-$prev_r) : $s[$n][4]*($s[$n][3]-$prev_q);
							if($tot_dis >15) {
								last;
							}
							$prev_r = $s[$n][2]; $prev_q = $s[$n][3];
						}
					}
					if($tot_dis <= 15 && $n < @s) {
						if(@div_reg > 0 && $div_reg[$#div_reg][0] eq $s[$m][0] && $div_reg[$#div_reg][2] <= $m && $div_reg[$#div_reg][3] >= $m) {
							$div_reg[$#div_reg][3] = $m+3;
						} else {
							@{$div_reg[@div_reg]} = ($s[$m][0], $s[$m][1], $m, $m+3);
						}
					}
				}
				for(local $i=@div_reg-1; $i >= 0; $i --) {
					local $r_base = substr($ref_seq, $s[$div_reg[$i][2]][7], $s[$div_reg[$i][3]][7] - $s[$div_reg[$i][2]][7] + 1);
					local $q_base = substr($qry_seq, $s[$div_reg[$i][2]][7], $s[$div_reg[$i][3]][7] - $s[$div_reg[$i][2]][7] + 1);
					$r_base =~ s/\.//g; $q_base =~ s/\.//g; $r_len = length($r_base); $q_len = length($q_base);
					for(local $j=0; $j<$q_len; ++$j) {
						$r_base .= ".";
					}
					for(local $j=0; $j<$r_len; ++$j) {
						$q_base = ".".$q_base;
					}
					substr($ref_seq, $s[$div_reg[$i][2]][7], $s[$div_reg[$i][3]][7] - $s[$div_reg[$i][2]][7] + 1, $r_base);
					substr($qry_seq, $s[$div_reg[$i][2]][7], $s[$div_reg[$i][3]][7] - $s[$div_reg[$i][2]][7] + 1, $q_base);
				}
				
				for(local $i=0; $i<length($ref_seq); ++$i) {
					local $r_base = substr($ref_seq, $i, 1);
					if($r_base eq ".") {
						substr($ref_seq, $i, 1, ""); substr($qry_seq, $i, 1, ""); $i --; next;
					}
					local $q_base = substr($qry_seq, $i, 1);
					if($r_base ne $q_base && $q_base ne "." && $q_base ne "n" && $r_base ne "n") {
						$site = "$ref:".($i+$ref_s);
						if($snp_list{"$site"} ne "") {
							substr($qry_seq, $i, 1, $snp_list{"$site"});
						} else {
							if($snp ne "") {
								exit;
							}
						}
					}
				}
				if(length($ref_seq) != $ref_e - $ref_s +1) {
				}
				substr($seq{$ref}, $ref_s-1, length($ref_seq), $qry_seq);
			}
		}
	}
}
close DELTA;

@name = keys %seq;
@name = sort @name;
foreach $n (@name) {
	print ">${n}_$delta\n$seq{$n}\n";
}
