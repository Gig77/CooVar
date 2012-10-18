#!/usr/bin/perl

use strict;
use File::Basename;

print "[tab2cv.pl] Start executing script on ";
system("date");

my $out_dir = $ARGV[1] or die "[tab2cv.pl] output directory not specified\n";

open(TAB,$ARGV[0]) || die "$!\n";

my $snp_file = "$out_dir/".basename($ARGV[0]) . '.snp';
my $ins_file = "$out_dir/".basename($ARGV[0]) . '.ins';
my $del_file = "$out_dir/".basename($ARGV[0]) . '.del';
 
open(SNP,">$snp_file") or die "[tab2cv.pl] could not write to $snp_file!\n";
open(INS,">$ins_file") or die "[tab2cv.pl] could not write to $ins_file!\n";
open(DEL,">$del_file") or die "[tab2cv.pl] could not write to $del_file!\n";

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
		print "[tab2cv.pl] WARNING: not SNP/INS/DEL\:" , $_,"\n";
	}
}
close(TAB);
close(SNP);
close(INS);
close(DEL);

print "[tab2cv.pl] $snps SNPs written to $snp_file...\n";
print "[tab2cv.pl] $ins insertions written to $ins_file...\n";
print "[tab2cv.pl] $del deletions written to $del_file...\n";

print "[tab2cv.pl] Done at ";
system("date");
