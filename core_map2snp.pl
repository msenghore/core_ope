#! /usr/bin/perl

$map = shift;
%translate = (
	"ttt" => "Phe", "ttc" => "Phe", "tta" => "Leu", "ttg" => "Leu",
	"tct" => "Ser", "tcc" => "Ser", "tca" => "Ser", "tcg" => "Ser",
	"tat" => "Tyr", "tac" => "Tyr", "taa" => "***", "tag" => "***",
	"tgt" => "Cys", "tgc" => "Cys", "tga" => "***", "tgg" => "Trp",
	
	"ctt" => "Leu", "ctc" => "Leu", "cta" => "Leu", "ctg" => "Leu",
	"cct" => "Pro", "ccc" => "Pro", "cca" => "Pro", "ccg" => "Pro",
	"cat" => "His", "cac" => "His", "caa" => "Gln", "cag" => "Gln",
	"cgt" => "Arg", "cgc" => "Arg", "cga" => "Arg", "cgg" => "Arg",
	
	"att" => "Ile", "atc" => "Ile", "ata" => "Ile", "atg" => "Met",
	"act" => "Thr", "acc" => "Thr", "aca" => "Thr", "acg" => "Thr",
	"aat" => "Asn", "aac" => "Asn", "aaa" => "Lys", "aag" => "Lys",
	"agt" => "Ser", "agc" => "Ser", "aga" => "Arg", "agg" => "Arg",
	
	"gtt" => "Val", "gtc" => "Val", "gta" => "Val", "gtg" => "Val",
	"gct" => "Ala", "gcc" => "Ala", "gca" => "Ala", "gcg" => "Ala",
	"gat" => "Asp", "gac" => "Asp", "gaa" => "Glu", "gag" => "Glu",
	"ggt" => "Gly", "ggc" => "Gly", "gga" => "Gly", "ggg" => "Gly"
	);


foreach $fanno (@ARGV) {
	$anno_file = $fanno;

	open (ANN, "<$anno_file") or die;
	while (<ANN>) {
		local @part = split /\t/, $_;
		if($part[2] eq "source") {
				next;
		}
		@{$annotation[@annotation]} = (0, $part[0], $part[2], $part[3], $part[4], $part[6]);

		for(local $i=int(($part[3]-100)/1000); $i<= int(($part[4]+100)/1000); ++$i) {
			${$site_in{$i}}[@{$site_in{$i}}] = scalar(@annotation)-1;
		}
	}
	close (ANN);
}


open (MAP, "<$map") or die;
while (<MAP>) {
	if( $_ =~ /^#/ ) {
		print $_;
		next;
	} elsif($_ =~ /^=/) {
		doit();

	} elsif( $_ =~ /^> *(.+?):(\d+)-(\d+)/ ) {
		$n[@n] = $1; $s[@s] = $2; $e[@e] = $3;
		if($s[$#s] > $e[$#e]) {
			$dir[@dir] = "-";
		} else {
			$dir[@dir] = "+";
		}
		if($_ =~ /^> *(.+?):(\d+)-(\d+) ([+-])/) {
			$dir[@dir-1] = $4;
			if($4 eq "-") {
				local $tmp = $s[$#s]; $s[$#s] = $e[$#e]; $e[$#e] = $tmp;
			}
		}
	} elsif ($_ =~ /^>(\S+)/) {
		$n[@n] = $1; $s[@s] = 1; $e[@e] = 1000000; $dir[@dir] = "+";
	}else {
		chomp $_;
		$seq[@n-1] .= $_;
	}
}

if(@n > 0) {
	doit();
}

sub doit {
		foreach $s1 (@seq) {
			$s1 =~ tr/A-Z/a-z/;
		}
		print "#Site";
		local $all_tag_ref;
		for(local $i=0; $i<@s; ++$i) {
			print "\t$n[$i]";
			$all_tag_ref .= "-";
		}
		print "\n";

		$realj = 0;

		for(local $j=0; $j<length($seq[0]); ++$j) {
			if(substr($seq[0], $j, 1) =~ /[a-zA-Z-]/) {
				$realj ++;
			}
			
			local $mutate = 0; local $prev = ""; local @this_base;
			local $all_tag = $all_tag_ref;
			for(local $i=0;$i<@s; ++$i) {
				$this_base[$i] = substr($seq[$i], $j, 1);
				if($this_base[$i] eq ".") {
					$this_base[$i] = "-";
				}
				if($prev eq "" && $this_base[$i] ne "-" && $this_base[$i] ne "n") {
					$prev = $this_base[$i];
				} elsif ($this_base[$i] ne "-" && $this_base[$i] ne "n" && $this_base[$i] ne "N") {
					if($prev ne $this_base[$i] || $this_base[$i] eq "-") {
						$mutate = 1;
					}
				}
			}

			if($mutate >=1) {
				local %present;
				local $tag = 65;
				for(local $ii=0; $ii<@s; ++$ii) {
					if($ii > 0) {
						print "\t".$this_base[$ii];
					} else {
						print "$realj\t".$this_base[$ii];
					}
					if($present{$ii} == 1 || $this_base[$ii] eq "-") {
						next;
					} else {
						substr($all_tag, $ii, 1, chr($tag));
					}

					for(local $jj=$ii+1; $jj<@s; ++$jj) {
						if($this_base[$ii] eq $this_base[$jj]) {
							substr($all_tag, $jj, 1, chr($tag));
							$present{$jj} = 1;
						}
					}
					$tag ++;
				}
				print "\t$all_tag";

				foreach (@{$site_in{int($realj/1000)}}) {
					local $a = $annotation[$_];
					$target_j = $realj;
					if($target_j >= $$a[3] && $target_j <= $$a[4]) {
						print "\t$$a[1]($$a[2]:$$a[5])\t";
					} elsif ($target_j >= $$a[3]-100 && $target_j < $$a[3] && $$a[5] eq "+") {
						print "\t$$a[1]($$a[2]:$$a[5])\t";
					} elsif ($target_j <= $$a[4]+100 && $target_j > $$a[3] && $$a[5] eq "-") {
						print "\t$$a[1]($$a[2]:$$a[5])\t";
					}else {
						next;
					}
					if($$a[2] eq "CDS" &&$target_j >= $$a[3] && $target_j <= $$a[4] ) {
						local $start_dist; local $start_site;
						if($$a[5] eq "+") {
							$start_dist = $realj - $$a[3];
							$start_site = $start_dist%3;
						} else {
							$start_dist = $$a[4] - $realj;
							$start_site = 2-$start_dist%3;
						}
						if($e[$$a[0]] < $s[$$a[0]]) {
							$start_site = 2 - $start_site;
						}
						print "$start_dist:$start_site\t";
						local %used; local @toaa;
						for(local $ii=0; $ii<@s; ++$ii) {
							local $tmp =substr($seq[$ii], $j-$start_site,3);
							if( ($$a[5] eq "-" && $e[$$a[0]] > $s[$$a[0]]) || ($$a[5] eq "+" && $e[$$a[0]] < $s[$$a[0]])) {
								$tmp = reverse $tmp;
								$tmp =~ tr/atcg/tagc/;
							}
							if($used{$tmp} < 1) {
								print "$tmp,";
								$toaa[@toaa] = $tmp;
								$used{$tmp} = 1;
							} 
						}
						print "\t";
						local $type = "S"; local $prev = "";
						foreach (@toaa) {
							local $tra = $translate{$_};
							if($tra eq "") {
								$tra = "IUB";
								if($_ =~ /-/) {
									$type = "FS";
								}
							}
							print $tra,",";
							if($prev ne "" && $tra ne $prev && $type ne "FS") {
								$type = "NS";
							}
							$prev = $tra;
						}
						print "\t$type";
					} else {
						local $start_dist;
						if($$a[5] eq "+") {
							$start_dist = $realj - $$a[3];
						} else {
							$start_dist = $$a[4] - $realj;
						}
						print "$start_dist:-\t";
						local %used; local $indel;
						for(local $ii=0; $ii < @s; ++$ii) {
							local $tmp = substr($seq[$ii], $j, 1);
							if($used{$tmp} < 1) {
								print "$tmp,";
								$used{$tmp} = 1;
								if($tmp eq "-") {
									$indel = 1;
								}
							} 
						}
						if($indel == 1) {
							print "\t-\tFS";
						} else {
							print "\t-\tNC";
							}
					}
				}
				print "\n";
			}
		}
		undef @n; undef @s; undef @e; undef @seq;
}
