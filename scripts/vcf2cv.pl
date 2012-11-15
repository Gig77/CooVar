#!/usr/bin/perl

use strict;
use File::Basename;

print "[vcf2cv.pl] Start executing script on ";
print localtime()."\n";

#this script takes as input a VCF file and generates the necessary input files for SNPs, insertions and deletions

print "[vcf2cv.pl] Parsing VCF input file on ";
print localtime()."\n";

my $out_dir = $ARGV[1] or die "[vcf2cv.pl] output directory not specified\n";
open(VCF,$ARGV[0]) || die "$!\n";
#my $filter_non_pass=$ARGV[1];

my $snp_file = "$out_dir/".basename($ARGV[0]) . '.snp';
my $ins_file = "$out_dir/".basename($ARGV[0]) . '.ins';
my $del_file = "$out_dir/".basename($ARGV[0]) . '.del';

open(SNP,">$snp_file");
open(INS,">$ins_file");
open(DEL,">$del_file");

my ($snps, $inss, $dels, $mult_alt) = (0, 0, 0, 0);
while(<VCF>)
{
	chomp($_);
	next if ($_=~/^\#/);
	my @line = split(/\t/,$_);
	#next if($line[6] eq 'PASS' && $filter_non_pass ne '0');
	my $chrom = $line[0];
	my $start = $line[1];
	my $ref = uc($line[3]);

	my @alt = split(/\,/,$line[4]);	
	$mult_alt ++ if (@alt > 1);
	my $tar = uc($alt[0]);

	# sanity check
	if (!$tar or !$ref)
	{
		print "[vcf2cv.pl]   WARNING: Invalid VCF entry: $_\n";
		next;
	};
	
	# ignore monomorphic variants (no change)
	if ($tar eq '.' or $ref eq $tar)
	{
		print "[vcf2cv.pl]   WARNING: Ignoring monomorphic variant: $_\n";
		next;
	};

	my $len_ref = length($ref);
	my $len_tar = length($tar);

	my ($p, $sr, $st) = (0, $len_ref, $len_tar);
	while ($p < $len_ref && $p < $len_tar && substr($ref, $p, 1) eq substr($tar, $p, 1)) { $p ++; } # skip matching prefix
	while ($sr > $p && $st > $p && substr($ref, $sr-1, 1) eq substr($tar, $st-1, 1)) { $sr--, $st--; } # skip matching suffix
		
	# SNPs = different bases b/w ref and target
	while($p < $sr and $p < $st) 
	{
		my ($br, $bt) = (substr($ref, $p, 1), substr($tar, $p, 1));
		if ($br ne $bt)
		{
			print SNP "$chrom\t".($start+$p)."\t$br\t$bt\n";
			$snps ++;											
		}
		$p ++;
	}
		
	# deletion = remaining sequence in ref
	if ($p < $sr)
	{
		print DEL "$chrom\t" , $start + $p  , "\t" , $start+$sr-1 , "\n";
		$dels ++;	
	}
		
	# insertion = remaining sequence in target
	if ($p < $st)
	{
		my $ins = substr($tar, $p, $st-$p);
		print INS "$chrom\t" , $start + $p - 1 , "\t" ,  $start + $p , "\t$ins\n";
		$inss ++;		
	}
}
close(VCF);
close(SNP);
close(INS);
close(DEL);

print "[vcf2cv.pl]   WARNING: $mult_alt variants had multiple alternative alleles; only first alternative allele will be considered\n"
	if ($mult_alt > 0);

print "[vcf2cv.pl]   $snps SNPs written to $snp_file...\n";
print "[vcf2cv.pl]   $inss insertions written to $ins_file...\n";
print "[vcf2cv.pl]   $dels deletions written to $del_file...\n";

print "[vcf2cv.pl] Done at ";
print localtime()."\n";
