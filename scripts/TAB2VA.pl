#!/usr/bin/perl

use strict;

open(TAB,$ARGV[0]) || die "$!\n";

my $snp_file = $ARGV[0] . '.VA_SNP';
my $ins_file = $ARGV[0] . '.VA_INS';
my $del_file = $ARGV[0] . '.VA_DEL';
 
open(SNP,">$snp_file");
open(INS,">$ins_file");
open(DEL,">$del_file");

while(<TAB>)
{
	chomp($_);
	$_=~/^(\S+)\s+(.+)$/;
	my $type = $1;
	my $rest = $2;

	if($type=~/SNP/i)
	{
		print SNP $rest,"\n";
	}
	elsif($type=~/INS/i)
	{
		print INS $rest,"\n";
	}
	elsif($type=~/DEL/i)
	{
		print DEL $rest,"\n";
	}
	else
	{
		print STDERR "Not SNP/INS/DEL\:" , $_,"\n";
	}
}
close(TAB);
close(SNP);
close(INS);
close(DEL);
