#!/usr/bin/perl

use strict;
use File::Basename;

print "[TAB2CV.pl] Start executing script on ";
system("date");

my $out_dir = $ARGV[1] or die "[TAB2CV.pl] output directory not specified\n";

open(TAB,$ARGV[0]) || die "$!\n";

my $snp_file = "$out_dir/".basename($ARGV[0]) . '.CV_SNP';
my $ins_file = "$out_dir/".basename($ARGV[0]) . '.CV_INS';
my $del_file = "$out_dir/".basename($ARGV[0]) . '.CV_DEL';
 
open(SNP,">$snp_file") or die "[TAB2CV.pl] could not write to $snp_file!\n";
open(INS,">$ins_file") or die "[TAB2CV.pl] could not write to $ins_file!\n";
open(DEL,">$del_file") or die "[TAB2CV.pl] could not write to $del_file!\n";

my ($snps, $ins, $del) = (0, 0, 0);
while(<TAB>)
{
	chomp($_);
	$_=~/^(\S+)\s+(.+)$/;
	my $type = $1;
	my $rest = $2;

	if($type=~/SNP/i)
	{
		print SNP $rest,"\n";
		$snps ++;
	}
	elsif($type=~/INS/i)
	{
		print INS $rest,"\n";
		$ins ++;
	}
	elsif($type=~/DEL/i)
	{
		print DEL $rest,"\n";
		$del ++;
	}
	else
	{
		print "[TAB2CV.pl] WARNING: not SNP/INS/DEL\:" , $_,"\n";
	}
}
close(TAB);
close(SNP);
close(INS);
close(DEL);

print "[TAB2CV.pl] $snps SNPs written to $snp_file...\n";
print "[TAB2CV.pl] $ins insertions written to $ins_file...\n";
print "[TAB2CV.pl] $del deletions written to $del_file...\n";

print "[TAB2CV.pl] Done at ";
system("date");
