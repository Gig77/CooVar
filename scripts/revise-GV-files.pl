#!/usr/bin/perl

use strict;

open(SNP,$ARGV[0]) || die "SNPs:$!\n";	
open(INS,$ARGV[1]) || die "Ins:$!\n";
open(DEL,$ARGV[2]) || die "Del:$!\n";

my $output_snp = "kept_" . $ARGV[0];
my $excluded_snp = "excluded_" . $ARGV[0];

open(KEPT_SNP,">$output_snp");
open(EXC_SNP,">$excluded_snp");

my $output_ins = "kept_" . $ARGV[1];
my $excluded_ins = "excluded_" . $ARGV[1];

open(KEPT_INS,">$output_ins");
open(EXC_INS,">$excluded_ins");

my $output_del = "kept_" . $ARGV[2];
my $excluded_del = "excluded_" . $ARGV[2];

open(KEPT_DEL,">$output_del");
open(EXC_DEL,">$excluded_del");

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
	}
	else
	{
		$_=~s/\s+/\t/g;
		print KEPT_SNP $_,"\n";
	}
}
close(SNP);
close(KEPT_SNP);
close(EXC_SNP);

while(<INS>)
{
        chomp($_);
        $_=~s/^\s+//;
        $_=~s/\s+$//;
        my @line = split(/\s+/,$_);
        if ((scalar(@line)!=4) || ($line[1]!~/^\d+$/) || ($line[2]!~/^\d+$/) || ($line[3]!~/[ACTGactg]/))
        {
                print EXC_INS $_,"\n";
        }
        else
        { 
		$_=~s/\s+/\t/g;
                print KEPT_INS $_,"\n";
       	}
}
close(INS);
close(KEPT_INS);
close(EXC_INS);

while(<DEL>)
{
        chomp($_);
        $_=~s/^\s+//;
        $_=~s/\s+$//;
        my @line = split(/\s+/,$_);
        if ((scalar(@line)!=3) || ($line[1]!~/^\d+$/) || ($line[2]!~/^\d+$/) || ($line[1]>$line[2]))
        {
                print EXC_DEL $_,"\n";
        }
	else
	{
		$_=~s/\s+/\t/g;
                print KEPT_DEL $_,"\n";
        }
}
close(DEL);
close(KEPT_DEL);
close(EXC_DEL);
