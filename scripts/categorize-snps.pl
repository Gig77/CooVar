#!/usr/bin/perl

use strict;

print "[categorize-snps.pl] Start executing script on ";
system("date");

my %ref_cdna=();
my %tar_cdna=();

sub classify_grantham
{
	my $grantham_score = shift;
	my $category="";

	if($grantham_score>150)
	{
		$category = "RADICAL";
	}
        elsif($grantham_score<=150 && $grantham_score>=101)
        {
		$category = "MODERATELY_RADICAL";
        }
        elsif($grantham_score<=100 && $grantham_score>=51)
        {
		$category = "MODERATELY_CONSERVATIVE";
        }
        elsif($grantham_score<=50 && $grantham_score>=1)
        {
		$category = "CONSERVATIVE";
        }
        else
        {
                print "[categorize-snps.pl] WARNING: Grantham score out of boundary: $grantham_score\n";
        }
	return $category;
}


open(D,$ARGV[0]) || die "[categorize-snps.pl] $!: $ARGV[0]\n";		#reference cDNA with start and end position for exons
open(E,$ARGV[1]) || die "[categorize-snps.pl] $!: $ARGV[1]\n";         #target cDNA with start and end position for exons
open(F,$ARGV[2]) || die "[categorize-snps.pl] $!: $ARGV[2]\n";		#Grantham Matrix
open(G,$ARGV[3]) || die "[categorize-snps.pl] $!: $ARGV[3]\n";		#all snps: will be used to categorized those as intronic/intergenic
open(H,$ARGV[4]) || die "[categorize-snps.pl] $!: $ARGV[4]\n";		#splice junction coordinates
open(I,$ARGV[5]) || die "[categorize-snps.pl] $!: $ARGV[5]\n";		#gff3_trans file for annotating SNPs
open(FREQ,">$ARGV[6]") || die "[categorize-snps.pl] $!: $ARGV[6]\n";  # snps frequency table
my $out_dir = $ARGV[7] or die "[categorize-snps.pl] output directory not specified\n";

my %translate = (
        # . - stop
        'taa'=>'*','tag'=>'*','tga'=>'*',
        # a - alanine
        'gct'=>'A','gcc'=>'A','gca'=>'A','gcg'=>'A',
        # c - cysteine
        'tgt'=>'C','tgc'=>'C',
        # d - aspartic acid
        'gat'=>'D','gac'=>'D',
        # e - glutamic acid
        'gaa'=>'E','gag'=>'E',
        # f - phenylalanine
        'ttt'=>'F','ttc'=>'F',
        # g - glycine
        'ggt'=>'G','ggc'=>'G','gga'=>'G','ggg'=>'G',
        # h - histidine
        'cat'=>'H','cac'=>'H',
        # i - isoleucine
        'att'=>'I','atc'=>'I','ata'=>'I',
        # k - lysine
        'aaa'=>'K','aag'=>'K',
        # l - leucine
        'ctt'=>'L','ctc'=>'L','cta'=>'L','ctg'=>'L',
        'tta'=>'L','ttg'=>'L',
        # m - methionine
        'atg'=>'M',
        # n - asparagine
        'aat'=>'N','aac'=>'N',
        # p - proline
        'cct'=>'P','ccc'=>'P','cca'=>'P','ccg'=>'P',
        # q - glutamine
        'caa'=>'Q','cag'=>'Q',
        # r - arginine
        'cgt'=>'R','cgc'=>'R','cga'=>'R','cgg'=>'R',
        'aga'=>'R','agg'=>'R',
        # s - serine
        'tct'=>'S','tcc'=>'S','tca'=>'S','tcg'=>'S',
        'agt'=>'S','agc'=>'S',
        # t - threonine
        'act'=>'T','acc'=>'T','aca'=>'T','acg'=>'T',
        # v - valine
        'gtt'=>'V','gtc'=>'V','gta'=>'V','gtg'=>'V',
        # w - tryptophan
        'tgg'=>'W',
        # y - tyrosine
        'tat'=>'Y','tac'=>'Y',
);

my %target_seq=();
my $target_id;
my %target_chrom=();

while(<E>)
{
        chomp($_);
        if($_=~/^\>(\S+)\t(\S+)\t\S+\t\S+\t(\S+)/)
        {
			$target_id=$1;
			$target_chrom{$target_id}=$2;
        }
        else
        {
            $target_seq{$target_id} .=$_ . "\n";
        }
}
close(E);

my %ref_seq=();
my $ref_id;

while(<D>)
{
        chomp($_);
        if($_=~/^\>(\S+)/)
        {
                $ref_id=$1;
        }
        else
        {
                $ref_seq{$ref_id} .=$_ . "\n";
        }
}
close(D);

my %grantham=();
my @amin;
my $count_rows = 0;
while(<F>)
{
	chomp($_);
	if($_=~/^A/)
	{
		@amin=split(//,$_);
	}
	else
	{
		last if ($_=~/\/\//);
		$_=~s/^\s+//;
		$_=~s/\s+$//;
		$_=~s/\.//g;
		my @scores = split(/\s+/,$_);
		for(my $i=0;$i<@scores;$i++)
		{
			my $coord1 = $amin[$count_rows] . "_" . $amin[$i];
			my $coord2 = $amin[$i] . "_" . $amin[$count_rows];
			$grantham{$coord1} = $scores[$i];
			$grantham{$coord2} = $scores[$i];
		}
		$count_rows++;
	}
}
close(F);

my %all_snps = ();
while(<G>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	$all_snps{"$line[0]\:$line[1]"}="$line[2] $line[3]";
}
close(G);

my %sj=();

while(<H>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	$line[1]=~/(\S+)\:(\d+)\.\.(\d+)/;
	my $chrom = $1;
	my $start = $2;
	my $end   = $3;
	$sj{"$chrom\:$start"}.=$line[0] . "\;$line[2]\;$line[4] ";
	$sj{"$chrom\:$end"}.=$line[0] . "\;$line[2]\;$line[4] ";
}
close(H);

open(CAT,">$out_dir/intermediate-files/categorized_snp_coords.list") || die "[categorize-snps.pl] $!: $out_dir/intermediate-files/categorized_snp_coords.list\n";
open(GVF,">$out_dir/categorized-gvs.gvf") || die "[categorize-snps.pl] $!: $out_dir/categorized-gvs.gvf\n";
print GVF "\#\#gff-version 3\n";
print GVF "\#\#gvf-version 1\.05\n";

my %seen_snps=();

my %snp2gvf_variant=();
my %snp2gvf_note=();
my $freqs = 0;

for my $key (keys %target_seq)
{
	my $synonymous=0;
	my $mis_sense=0;
	my $non_sense=0;
	my $extension=0;

	my @exons_tar=split(/\n/,$target_seq{$key});
	my @exons_ref=split(/\n/,$ref_seq{$key});

	my @snps;

	for (my $i=0;$i<@exons_tar;$i++)
	{
		$exons_tar[$i]=~/\((\d+)\-(\d+)\)\s+(\S+)/;
		my $start=$1;
		my $end=$2;
		my $subseq=$3;
		$tar_cdna{$key} .=$subseq;

		my @aux_subseq=split(//,$subseq);
		for(my $j=0;$j<@aux_subseq;$j++)
		{
			next if ($aux_subseq[$j]!~/[ACTG]/);
			my $coord;
			if($start<$end)
			{
				$coord=$start+$j;
			}
			else
			{
				$coord=$start-$j;
			}
			my $snp="$target_chrom{$key}\:$coord";
			push @snps,$snp;
		}
	}

        for (my $i=0;$i<@exons_ref;$i++)
        {
                $exons_ref[$i]=~/\((\d+)\-(\d+)\)\s+(\S+)/;
                my $start=$1;
                my $end=$2;
                $ref_cdna{$key} .=$3;
        }

	my @seq_ref=split(//,$ref_cdna{$key});
	my @seq_tar=split(//,$tar_cdna{$key});

	if(scalar(@seq_ref) != scalar(@seq_tar))
	{
		print "[categorize-snps.pl] WARNING: Inconsistent amount of DNA for $key\n";
	}

	my $count_codons=0;

	for(my $i=0;$i<@seq_ref;$i+=3)
	{
		my $n_subs=0;
		$count_codons++;
	        my @codon_bias;
        	$codon_bias[0] = 0;
        	$codon_bias[1] = 0;
        	$codon_bias[2] = 0;

		for(my $j=0;$j<=2;$j++)
		{
			next if ($seq_tar[$i+$j]!~/[ACTG]/);
			$codon_bias[$j]++;
			$n_subs++;
		}

		my $codon_ref= $seq_ref[$i] . "" .  $seq_ref[$i+1] . "" . $seq_ref[$i+2];
		my $codon_tar= $seq_tar[$i] . "" .  $seq_tar[$i+1] . "" . $seq_tar[$i+2];
		
		my $amino_ref=$translate{$codon_ref};
		my $amino_tar=$translate{lc($codon_tar)};

		if($codon_tar=~/[ACTG]/)
		{
			if($amino_ref eq $amino_tar)
			{
				$synonymous+=$n_subs;
				for(my $k=0;$k<@codon_bias;$k++)
				{
					next if ($codon_bias[$k] == 0);
					my $snp=shift @snps;

					$snp2gvf_variant{"$snp\_\_synonymous\_codon"}.= $key . "\,";
					$snp2gvf_note{"$snp\_\_synonymous\_codon"}.="$codon_ref\>$codon_tar\_$amino_ref\>$amino_tar\_"; 

					my $codon_pos = $k+1;
					$snp2gvf_note{"$snp\_\_synonymous\_codon"}.="aa$count_codons\_";
					$snp2gvf_note{"$snp\_\_synonymous\_codon"}.="codon_loc$codon_pos\,";

					print CAT "$snp\t$key\tsynonymous\t$codon_ref\>$codon_tar\t$amino_ref\>$amino_tar\t";
					print CAT "\t\t", $k+1,"\n";
					$seen_snps{$snp}++;
				}
			}
			elsif(($amino_ref ne $amino_tar) &&  ($amino_tar eq '*'))
			{
				$non_sense+=$n_subs;
				for(my $k=0;$k<@codon_bias;$k++)
				{
					next if ($codon_bias[$k] == 0);
					my $snp=shift @snps;

					$snp2gvf_variant{"$snp\_\_stop\_gained"}.= $key . "\,";
                    $snp2gvf_note{"$snp\_\_stop\_gained"}.="$codon_ref\>$codon_tar\_$amino_ref\>$amino_tar\_";
					my $codon_pos = $k+1;
                                        $snp2gvf_note{"$snp\_\_stop\_gained"}.="aa$count_codons\_";
                                        $snp2gvf_note{"$snp\_\_stop\_gained"}.="codon_loc$codon_pos\,";

					print CAT "$snp\t$key\tnon_sense\t$codon_ref\>$codon_tar\t$amino_ref\>$amino_tar\t";
					print CAT "\t\t", $k+1,"\n";
					$seen_snps{$snp}++;
				}
			}
			elsif(($amino_ref ne $amino_tar) &&  ($amino_ref eq '*'))
			{
				$extension+=$n_subs;
                                for(my $k=0;$k<@codon_bias;$k++)
                                {
					next if ($codon_bias[$k] == 0);
                                        my $snp=shift @snps;

					$snp2gvf_variant{"$snp\_\_stop\_lost"}.= $key . "\,";
					$snp2gvf_note{"$snp\_\_stop\_lost"}.="$codon_ref\>$codon_tar\_$amino_ref\>$amino_tar\_";
                    my $codon_pos = $k+1;
                    $snp2gvf_note{"$snp\_\_stop\_lost"}.="aa$count_codons\_";
                    $snp2gvf_note{"$snp\_\_stop\_lost"}.="codon_loc$codon_pos\,";

                    print CAT "$snp\t$key\textension\t$codon_ref\>$codon_tar\t$amino_ref\>$amino_tar\t"; 
					print CAT "\t\t", $k+1,"\n";
					$seen_snps{$snp}++;
                                }

			}
			else
			{
				$mis_sense+=$n_subs;
                for(my $k=0;$k<@codon_bias;$k++)
                {
					next if ($codon_bias[$k] == 0);
                    my $snp=shift @snps;

					my $coord_grantham = $amino_ref . "_" . $amino_tar;
                    print CAT "$snp\t$key\tmis_sense\t$codon_ref\>$codon_tar\t$amino_ref\>$amino_tar\t"; 
					print CAT "$count_codons\t";
					my $grantham_category = classify_grantham($grantham{$coord_grantham});
					print CAT "$grantham_category\($grantham{$coord_grantham}\)\t";
					print CAT $k+1,"\n";

					if($grantham_category eq 'CONSERVATIVE' || $grantham_category eq 'MODERATELY_CONSERVATIVE')
					{
						$snp2gvf_variant{"$snp\_\_conservative\_missense\_codon"}.= $key . "\,";
						$snp2gvf_note{"$snp\_\_conservative\_missense\_codon"}.="$codon_ref\>$codon_tar\_$amino_ref\>$amino_tar\_";
                        my $codon_pos = $k+1;
                        $snp2gvf_note{"$snp\_\_conservative\_missense\_codon"}.="aa$count_codons\_";
                        $snp2gvf_note{"$snp\_\_conservative\_missense\_codon"}.="codon_loc$codon_pos\_";
						$snp2gvf_note{"$snp\_\_conservative\_missense\_codon"}.="$grantham_category\($grantham{$coord_grantham}\)\,";
					}
					else
					{
						$snp2gvf_variant{"$snp\_\_non\_conservative\_missense\_codon"}.= $key . "\,";
                        $snp2gvf_note{"$snp\_\_non\_conservative\_missense\_codon"}.="$codon_ref\>$codon_tar\_$amino_ref\>$amino_tar\_";
                        my $codon_pos = $k+1;
                        $snp2gvf_note{"$snp\_\_non\_conservative\_missense\_codon"}.="aa$count_codons\_";
                        $snp2gvf_note{"$snp\_\_non\_conservative\_missense\_codon"}.="codon_loc$codon_pos\_";
                        $snp2gvf_note{"$snp\_\_non\_conservative\_missense\_codon"}.="$grantham_category\($grantham{$coord_grantham}\)\,";

					}

					$seen_snps{$snp}++;
                                }
			}
		}
	}

	print FREQ "$key\t", length($ref_cdna{$key}), "\t$synonymous\t$mis_sense\t$non_sense\t$extension\t";
	print FREQ 'Total=' , $synonymous + $mis_sense + $non_sense + $extension , "\n";
	$freqs ++;
}
close(FREQ);
print "[categorize-snps.pl] $freqs lines written to $ARGV[6]\n";

#find SNPs that impact splice junctions
print "[categorize-snps.pl] Determining SNPs impacting splice junctions\n";
for my $key (keys %all_snps)
{
	#next if (defined $seen_snps{$key});
	if(defined $sj{$key})
	{
		my @trans=split(/\s+/,$sj{$key});
		for (my $i=0;$i<@trans;$i++)
		{
			$trans[$i]=~/(\S+)\;(\S+)\;(\S+)/;
			my $trans_name = $1;
			my $sj_seq = $2;
			my $type_sj = $3;
			if($type_sj eq 'donor')
			{
				$snp2gvf_variant{"$key\_\_splice\_donor\_variant"}.= $trans_name . "\,";
			}
			else
			{
				$snp2gvf_variant{"$key\_\_splice\_acceptor\_variant"}.= $trans_name . "\,";
			}
			print CAT "$key\t$trans_name\tsplice_junction\t\t\t\t\t$type_sj\n";
		}
		$seen_snps{$key}++;
	}
}

for my $key (keys %all_snps)
{
	next if (defined $seen_snps{$key});
	$snp2gvf_variant{"$key\_\_silent\_mutation"}++;
	print CAT "$key\t\tintronic\/intergenic\t\t\t\t\t\n";
}
close(CAT);

my $snp_id=1;

my %transcript2snp_id=();

my $gvfs = 0;
for my $key (keys %snp2gvf_variant)
{
	my @snp_type = split(/\_\_/,$key);
	$snp_type[0]=~/(\S+)\:(\d+)/;
	my $chrom = $1;
	my $coord = $2;

	$all_snps{$snp_type[0]}=~/(\S+)\s+(\S+)/;
	my $ref = $1;
	my $tar = $2;

	print GVF "$chrom\tCooVar\tSNV\t$coord\t$coord\t\.\t\+\t\.\t";
	print GVF "ID\=snp\_$snp_id\;Variant\_seq=$tar\;Reference\_seq\=$ref\;";
	print GVF "Variant\_type\=$snp_type[1];";  # 2011-11-21 | CF | additional tag for easy track grouping in GBrowse
	print GVF "Variant\_effect\=$snp_type[1]";

	if($snp_type[1] ne 'silent_mutation')
	{
		chop($snp2gvf_variant{$key});	#removes comma
		print GVF " 0 mRNA $snp2gvf_variant{$key}\;";

        my @trans = split(/\,/,$snp2gvf_variant{$key});
        for(my $i=0;$i<@trans;$i++)
        { 
        	$transcript2snp_id{$trans[$i]}.="snp\_$snp_id\($snp_type[1]\),";
        }
	}

	if(defined $snp2gvf_note{$key})
	{
		chop($snp2gvf_note{$key});	#removes comma
		print GVF "Note\=$snp2gvf_note{$key}\n";
		$gvfs ++;
	}
	else
	{
		print GVF "\n";
		$gvfs ++;
	}
	$snp_id++;
}
close(GVF);

print "[categorize-snps.pl] $gvfs GVs written to $out_dir/categorized-gvs.gvf\n";

#now we annotate the SNPs IDs in the gff3 file for transcripts

open(TMP_GFF3,">$out_dir/intermediate-files/transcripts_snps_applied.gff3.tmp") || die "[categorize-snps.pl] $!: $out_dir/intermediate-files/transcripts_snps_applied.gff3.tmp\n";
while(<I>)
{
	chomp($_);
	if($_=~/^\#/ || $_=~/\tCooVar\tCDS\t/)
	{
		print TMP_GFF3 $_,"\n";
		next;
	}
	else
	{
		$_=~/\tID\=([^;\s]+)/;
        print "[categorize-snps.pl] ERROR: could not parse ID from following line:\n$_\n" if (!defined $1);
		my $trans = $1;
		print TMP_GFF3 $_,"\;Note\=$transcript2snp_id{$trans}\n";
	}
}
close(I);
#system("rm transcripts.gff3.tmp")

print "[categorize-snps.pl] Done at ";
system("date");
