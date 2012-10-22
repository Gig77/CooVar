#!/usr/bin/perl

use strict;
use File::Basename;

print "[revise-gv-files.pl] Start executing script on ";
system("date");

my $out_dir = $ARGV[3];

# reading contig information
my %contigs;
print "[revise-gv-files.pl] Reading contig information from $out_dir/intermediate-files/contigs.summary on ";
system("date");
open(CHR, "$out_dir/intermediate-files/contigs.summary")
	or die("[revise-gv-files.pl]   ERROR: Could not read contig information.\n");
while(<CHR>)
{
	chomp;
	my ($contig, $len) = split("\t");
	$contigs{$contig} = $len;
}
close(CHR);

open(SNP,$ARGV[0]) || die "[revise-gv-files.pl] $!: $ARGV[0]\n";	
open(INS,$ARGV[1]) || die "[revise-gv-files.pl] $!: $ARGV[1]\n";
open(DEL,$ARGV[2]) || die "[revise-gv-files.pl] $!: $ARGV[2]\n";

my $output_snp = "$out_dir/snps/kept_" . basename($ARGV[0]);
my $excluded_snp = "$out_dir/snps/excluded_" . basename($ARGV[0]);

open(KEPT_SNP,">$output_snp");
open(EXC_SNP,">$excluded_snp");

my $output_ins = "$out_dir/insertions/kept_" . basename($ARGV[1]);
my $excluded_ins = "$out_dir/insertions/excluded_" . basename($ARGV[1]);

open(KEPT_INS,">$output_ins");
open(EXC_INS,">$excluded_ins");

my $output_del = "$out_dir/deletions/kept_" . basename($ARGV[2]);
my $excluded_del = "$out_dir/deletions/excluded_" . basename($ARGV[2]);

open(KEPT_DEL,">$output_del");
open(EXC_DEL,">$excluded_del");

# revising SNPs
print "[revise-gv-files.pl] Revising SNPs... on ";
system("date");

my $kept = 0;
my $excluded = 0;
while(<SNP>)
{
	chomp($_);
	$_=~s/^\s+//;
	$_=~s/\s+$//;
	my @line = split(/\s+/,$_);
	if ((scalar(@line)!=4) || ($line[1]!~/^\d+$/))
	{
		print EXC_SNP $_,"\terror parsing line\n";
		$excluded ++;
	}
	elsif (!exists $contigs{$line[0]})
	{
		print EXC_SNP $_,"\tinvalid contig name\n";
		$excluded ++;
	}
	elsif ($line[1] > $contigs{$line[0]})
	{
		print EXC_SNP $_,"\tposition exceeds contig length\n";
		$excluded ++;		
	}
	elsif (($line[2]!~/[ACTGactg]/) || ($line[2]!~/\w/))
	{
		print EXC_SNP $_,"\tinvalid reference allele\n";
		$excluded ++;
	}
	if (($line[3]!~/[ACTGactg]/) || ($line[3]!~/\w/))
	{
		print EXC_SNP $_,"\tinvalid alternative allele\n";
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

if ($excluded > 0)
{
	print "[revise-gv-files.pl]   **************\n";
	print "[revise-gv-files.pl]   *** WARNING: $excluded SNPs have been excluded from analysis! Check file $excluded_snp for additional information.\n";
	print "[revise-gv-files.pl]   **************\n";
}
print "[revise-gv-files.pl]   $kept kept SNPs written to $output_snp\n";

# revising insertions
print "[revise-gv-files.pl] Revising insertions... on ";
system("date");

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
		elsif (!exists $contigs{$line[0]})
		{
			print EXC_INS $_,"\tinvalid contig name\n";
			$excluded ++;
		}
		elsif ($line[1] > $contigs{$line[0]})
		{
			print EXC_INS $_,"\tposition exceeds contig length\n";
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

if ($excluded > 0)
{
	print "[revise-gv-files.pl]   **************\n";
	print "[revise-gv-files.pl]   *** WARNING: $excluded insertions have been excluded from analysis! Check file $excluded_ins for additional information.\n";
	print "[revise-gv-files.pl]   **************\n";
}
print "[revise-gv-files.pl]   $kept kept insertions written to $output_ins\n";

# revising deletions
print "[revise-gv-files.pl] Revising deletions... on ";
system("date");

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
		elsif (!exists $contigs{$line[0]})
		{
			print EXC_DEL $_,"\tinvalid contig name\n";
			$excluded ++;
		}
		elsif ($line[1] > $contigs{$line[0]})
		{
			print EXC_DEL $_,"\tposition exceeds contig length\n";
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

if ($excluded > 0)
{
	print "[revise-gv-files.pl]   **************\n";
	print "[revise-gv-files.pl]   *** WARNING: $excluded deletions have been excluded from analysis! Check file $excluded_del for additional information.\n";
	print "[revise-gv-files.pl]   **************\n";
}
print "[revise-gv-files.pl]   $kept kept deletions written to $output_del\n";

print "[revise-gv-files.pl] Done at ";
system("date");
