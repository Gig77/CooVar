#!/usr/bin/perl

use strict;

my $reference_pfa = $ARGV[0];
my $snp_list = $ARGV[1];

open(D,$reference_pfa) || die "$!\n";

open(E,$snp_list) || die "$!\n";

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

open(SNP_PFA,">target_cDNA_SNPs.exons");
open(SNP_FA,">target_cDNA_SNPs.fasta");

for my $key (keys %trans_exons)
{
	print SNP_PFA "\>$key\t$trans_chrom{$key}\t$trans_start{$key}\t$trans_end{$key}\t$trans_strand{$key}\n";
	print SNP_FA "\>$key\t$trans_chrom{$key}\t$trans_start{$key}\t$trans_end{$key}\t$trans_strand{$key}\n";

	my @exons = split(/\n/,$trans_exons{$key});
	for(my $j=0;$j<@exons;$j++)
        {
        	my ($exon_prefix,$exon_start,$exon_end, $exon_seq);
        	$exons[$j]=~/^(\((\d+)\-(\d+)\))\s+(\S+)/;
                $exon_prefix=$1;
                $exon_start = $2;
                $exon_end = $3;
                $exon_seq = $4;

		print SNP_PFA "$exon_prefix ";

		my $abs_coord=$exon_start;
		my @nts= split(//,$exon_seq);

		for (my $i=0;$i<@nts;$i++)
		{
			my $coord = "$trans_chrom{$key}\:$abs_coord";
			if(defined $snps{$coord})
			{
				my $snp_seq = $snps{$coord};
				if($trans_strand{$key} eq '-')
				{
					$snp_seq =~tr/[ACTG]/[TGAC]/;
				}
                                print SNP_PFA $snp_seq;
                                print SNP_FA $snp_seq;
			}
			else
			{
				print SNP_PFA $nts[$i];
				print SNP_FA $nts[$i];
			}
			if($trans_strand{$key} eq '-')
			{
				$abs_coord--;
			}
			else
			{
				$abs_coord++;
			}
		}
		print SNP_PFA "\n";
	}
	print SNP_FA "\n";
}
close(SNP_PFA);
close(SNP_FA);
