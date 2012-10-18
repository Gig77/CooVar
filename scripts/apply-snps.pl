#!/usr/bin/perl

use strict;

print "[apply_SNPs.pl] Start executing script on ";
system("date");

my $reference_pfa = $ARGV[0];
my $snp_list = $ARGV[1];
my $out_dir = $ARGV[2] or die "[apply_SNPs.pl] output directory not specified\n";

open(D,$reference_pfa) || die "[apply_SNPs.pl] $!: $reference_pfa\n";

open(E,$snp_list) || die "[apply_SNPs.pl] $!: $snp_list\n";

my %snps = ();

my %trans_exons= ();
my %trans_start= ();
my %trans_end  = ();
my %trans_chrom= ();
my %trans_strand=();

my $id;

while(<D>)
{
        chomp($_);
        if($_=~/^\>(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
        {
                $id=$1;
                my $chrom=$2;
                my $start = $3;
                my $end = $4;
                my $strand = $5;
                $trans_chrom{$id} = $chrom;
                $trans_start{$id} = $start;
                $trans_end{$id} = $end;
                $trans_strand{$id} = $strand;
        }
        else
        {
                $trans_exons{$id}.=$_ . "\n";
        }
}
close(D);

while(<E>)
{
        chomp($_);
        my @line = split(/\t/,$_);
        my $chrom = $line[0];
        my $coord = $line[1];
        my $subs  = $line[3];
        $snps{"$chrom\:$coord"} = uc($subs);
}
close(E);

open(SNP_PFA,">$out_dir/intermediate-files/target_cdna_snps.exons") or die ("[apply_SNPs.pl] ERROR writing to file $out_dir/intermediate-files/target_cdna_snps.exons\n");
open(SNP_FA,">$out_dir/intermediate-files/target_cdna_snps.fasta") or die ("[apply_SNPs.pl] ERROR writing to file $out_dir/intermediate-files/target_cdna_snps.fasta\n");

my ($snp_pfa_lines, $snp_fa_lines) = (0, 0);
my $c;
for my $key (keys %trans_exons)
{
	my $chrom = $trans_chrom{$key};
	my $strand = $trans_strand{$key};
	
	print SNP_PFA "\>$key\t$chrom\t$trans_start{$key}\t$trans_end{$key}\t$strand\n";
	print SNP_FA "\>$key\t$chrom\t$trans_start{$key}\t$trans_end{$key}\t$strand\n";

	my @exons = split(/\n/,$trans_exons{$key});
	for(my $j=0;$j<@exons;$j++)
    {
      	my ($exon_prefix,$exon_start,$exon_end, $exon_seq) = $exons[$j] =~ /^(\((\d+)\-(\d+)\))\s+(\S+)/;
		print SNP_PFA "$exon_prefix ";

		my $abs_coord=$exon_start;
		for (my $i=0;$i<length($exon_seq);$i++)
		{
			my $coord = "$chrom\:$abs_coord";
			if(defined $snps{$coord})
			{
				my $snp_seq = $snps{$coord};
				if($strand eq '-')
				{
					$snp_seq =~tr/[ACTG]/[TGAC]/;
				}
                print SNP_PFA $snp_seq;
                print SNP_FA $snp_seq;
			}
			else
			{
				$c = substr($exon_seq, $i, 1);
				print SNP_PFA $c;
				print SNP_FA $c;
			}
			if($strand eq '-')
			{
				$abs_coord--;
			}
			else
			{
				$abs_coord++;
			}
		}
		print SNP_PFA "\n";
		$snp_pfa_lines ++;
	}
	print SNP_FA "\n";
	$snp_fa_lines ++;
}
close(SNP_PFA);
close(SNP_FA);

print "[apply_SNPs.pl] $snp_pfa_lines lines written to $out_dir/intermediate-files/target_cdna_snps.exons\n";
print "[apply_SNPs.pl] $snp_fa_lines lines written to $out_dir/intermediate-files/target_cdna_snps.fasta\n";

print "[apply_SNPs.pl] Done at ";
system("date");
