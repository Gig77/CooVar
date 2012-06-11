#!/usr/bin/perl

use strict;
use File::Basename;

print "[revise-GV-files.pl] Start executing script on ";
system("date");

my $out_dir = $ARGV[3];

open(SNP,$ARGV[0]) || die "[revise-GV-files.pl] $!: $ARGV[0]\n";	
open(INS,$ARGV[1]) || die "[revise-GV-files.pl] $!: $ARGV[1]\n";
open(DEL,$ARGV[2]) || die "[revise-GV-files.pl] $!: $ARGV[2]\n";

my $output_snp = "$out_dir/VA_SNPs/kept_" . basename($ARGV[0]);
my $excluded_snp = "$out_dir/VA_SNPs/excluded_" . basename($ARGV[0]);

open(KEPT_SNP,">$output_snp");
open(EXC_SNP,">$excluded_snp");

my $output_ins = "$out_dir/VA_Insertions/kept_" . basename($ARGV[1]);
my $excluded_ins = "$out_dir/VA_Insertions/excluded_" . basename($ARGV[1]);

open(KEPT_INS,">$output_ins");
open(EXC_INS,">$excluded_ins");

my $output_del = "$out_dir/VA_Deletions/kept_" . basename($ARGV[2]);
my $excluded_del = "$out_dir/VA_Deletions/excluded_" . basename($ARGV[2]);

open(KEPT_DEL,">$output_del");
open(EXC_DEL,">$excluded_del");

my $kept = 0;
my $excluded = 0;
while(<SNP>)
{
	chomp($_);
	$_=~s/^\s+//;
	$_=~s/\s+$//;
	my @line = split(/\s+/,$_);
	if ((scalar(@line)!=4) || ($line[1]!~/^\d+$/) || ($line[2]!~/[ACTGactg]/) || ($line[2]!~/\w/)
	|| ($line[3]!~/[ACTGactg]/) || ($line[3]!~/\w/))
	{
		print EXC_SNP $_,"\n";
		$excluded ++;
	}
	else
	{
		$_=~s/\s+/\t/g;
		print KEPT_SNP $_,"\n";
		$kept ++;
	}
}
close(SNP);
close(KEPT_SNP);
close(EXC_SNP);

print "[revise-GV-files.pl] $kept kept SNPs written to $output_snp\n";
print "[revise-GV-files.pl] $excluded excluded SNPs written to $excluded_snp\n";


$kept = 0;
$excluded = 0;
while(<INS>)
{
        chomp($_);
        $_=~s/^\s+//;
        $_=~s/\s+$//;
        my @line = split(/\s+/,$_);
        if ((scalar(@line)!=4) || ($line[1]!~/^\d+$/) || ($line[2]!~/^\d+$/) || ($line[3]!~/[ACTGactg]/))
        {
            print EXC_INS $_,"\n";
			$excluded ++;
        }
        else
        { 
			$_=~s/\s+/\t/g;
        	print KEPT_INS $_,"\n";
			$kept ++;
       	}
}
close(INS);
close(KEPT_INS);
close(EXC_INS);

print "[revise-GV-files.pl] $kept kept insertions written to $output_ins\n";
print "[revise-GV-files.pl] $excluded excluded insertions written to $excluded_ins\n";

$kept = 0;
$excluded = 0;
while(<DEL>)
{
        chomp($_);
        $_=~s/^\s+//;
        $_=~s/\s+$//;
        my @line = split(/\s+/,$_);
        if ((scalar(@line)!=3) || ($line[1]!~/^\d+$/) || ($line[2]!~/^\d+$/) || ($line[1]>$line[2]))
        {
        	print EXC_DEL $_,"\n";
			$excluded ++;
        }
		else
		{
			$_=~s/\s+/\t/g;
            print KEPT_DEL $_,"\n";
			$kept ++;
        }
}
close(DEL);
close(KEPT_DEL);
close(EXC_DEL);

print "[revise-GV-files.pl] $kept kept deletions written to $output_del\n";
print "[revise-GV-files.pl] $excluded excluded deletions written to $excluded_del\n";

print "[revise-GV-files.pl] Done at ";
system("date");
