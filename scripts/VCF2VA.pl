#!/usr/bin/perl

use strict;
use File::Basename;

print "[VCF2VA.pl] Start executing script on ";
system("date");

#this script takes as input a VCF file and generates the necessary input files for SNPs, insertions and deletions

my $out_dir = $ARGV[1] or die "[VCF2VA.pl] output directory not specified\n";
open(VCF,$ARGV[0]) || die "$!\n";
#my $filter_non_pass=$ARGV[1];

my $snp_file = "$out_dir/".basename($ARGV[0]) . '.VA_SNP';
my $ins_file = "$out_dir/".basename($ARGV[0]) . '.VA_INS';
my $del_file = "$out_dir/".basename($ARGV[0]) . '.VA_DEL';

open(SNP,">$snp_file");
open(INS,">$ins_file");
open(DEL,">$del_file");

my ($snps, $inss, $dels) = (0, 0, 0);
while(<VCF>)
{
	chomp($_);
	next if ($_=~/^\#/);
	my @line = split(/\t/,$_);
	#next if($line[6] eq 'PASS' && $filter_non_pass ne '0');
	my $chrom = $line[0];
	my $start = $line[1];
	my $ref = $line[3];
	my @alt = split(/\,/,$line[4]);
	my $tar = $alt[0];	#arbitrarily takes the first alternative
	next if ($tar eq '.');	#monomorphic, no change

	my $len_ref = length($ref);
	my $len_tar = length($tar);

	#SNP
	if($len_ref == 1 && $len_tar == 1)
	{
		print SNP "$chrom\t$start\t$ref\t$tar\n";
		$snps ++;	
	}
	#SNP shown in an indel manner
	elsif($len_ref == $len_tar && $len_ref > 1)
	{
		print "[VCF2VA.pl] WARNING: Check $chrom\t$start\n";
		$ref=~/^\S(\S+)$/;
		my $original = $1;
		$tar=~/^\S(\S+)$/;
		my $subs = $1;
		my @to_original = split(//,$original);
		my @to_subs = split(//,$subs);
		my $aux_start = $start+1;
		for(my $i=0;$i<@to_subs;$i++)
		{
			next if ($to_original[$i] eq $to_subs[$i]);	#within the substituted block could be same base pairs
			print SNP "$chrom\t$aux_start\t$to_original[$i]\t$to_subs[$i]\n";
			$aux_start++;
			$snps ++;			
		}
	}
	#DEL
	elsif($len_ref > $len_tar)
	{
		print DEL "$chrom\t" , $start + $len_tar  , "\t" , $start + $len_ref - 1 , "\n";
		$dels ++;	
	}
	#INS
	elsif($len_ref < $len_tar)
	{
		$ref=~/^(\S)(\S*)$/;
		my $ref_base = $1;
		my $ref_rest = $2;
		$tar=~/^$ref_base(\S+)$ref_rest$/;
		my $ins = $1;
		print INS "$chrom\t" , $start  , "\t" ,  $start + 1 , "\t$ins\n";
		$inss ++;	
	}
	else
	{
		print "[VCF2VA.pl] WARNING: Unrecognized line: $_\n";
	}
}
close(VCF);
close(SNP);
close(INS);
close(DEL);

print "[VCF2VA.pl] $snps SNPs written to $snp_file...\n";
print "[VCF2VA.pl] $inss insertions written to $ins_file...\n";
print "[VCF2VA.pl] $dels deletions written to $del_file...\n";

print "[VCF2VA.pl] Done at ";
system("date");
