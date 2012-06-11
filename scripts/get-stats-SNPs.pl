#!/usr/bin/perl

use strict;

print "[get-stats-SNPs.pl] Start executing script on ";
system(date);

open(D,$ARGV[0]) || die "[get-stats-SNPs.pl] $!: $ARGV[0]\n";
my $out_dir = $ARGV[1] or die "[get-stats-SNPs.pl] output directory not specified\n";

open(CB,">$out_dir/VA_SNPs/codon_bias.stat") || die "[get-stats-SNPs.pl] $!: $out_dir/VA_SNPs/codon_bias.stat\n";
open(TYPE,">$out_dir/VA_SNPs/categorized_SNPs.stat") || die "[get-stats-SNPs.pl] $!: $out_dir/VA_SNPs/categorized_SNPs.stat\n";
open(AMB,">$out_dir/VA_SNPs/ambiguous_SNPs.list") || die "[get-stats-SNPs.pl] $!: $out_dir/VA_SNPs/ambiguous_SNPs.list\n";


my %count_syn=();
my %count_mis_sense=();
my %count_non_sense=();
my %count_extension=();
my %count_codon_bias=();
my %count_splice_junction=();
my %count_intronic_intergenic=();

my %snp2type=();
my %to_exclude_ambiguous=();

my @file;

while(<D>)
{
	chomp($_);
	push @file,$_;
}
close(D);

for(my $i=0;$i<@file;$i++)
{
        my @line=split(/\t/,$file[$i]);
        my $coord=$line[0];
        my $type=$line[2];
	if(!(defined $snp2type{$coord}))
	{
		$snp2type{$coord}=$type;
	}
	else
	{
		if($snp2type{$coord} ne $type)
		{
			$to_exclude_ambiguous{$coord}++;
		}
	}
}

my $ambs = 0;
for(my $i=0;$i<@file;$i++)
{
	my @line=split(/\t/,$file[$i]);
	my $coord=$line[0];
	my $type=$line[2];

	if(defined $to_exclude_ambiguous{$coord})
	{
		print AMB $file[$i],"\n";
		$ambs ++;
		next;
	}
	my $codon_pos= $line[7];
	
	if($type eq 'extension')
	{
		$count_extension{$coord}++;	
		$count_codon_bias{"$codon_pos\tnon_synonymous"}++;
	}
	elsif($type eq 'non_sense')
	{
		$count_non_sense{$coord}++;
		$count_codon_bias{"$codon_pos\tnon_synonymous"}++;
	}
	elsif($type eq 'mis_sense')
	{
		$count_mis_sense{$coord}++;
		$count_codon_bias{"$codon_pos\tnon_synonymous"}++;
	}
	elsif($type eq 'synonymous')
	{
		$count_syn{$coord}++;
		$count_codon_bias{"$codon_pos\tsynonymous"}++;
	}
	elsif($type eq 'splice_junction')
	{
		$count_splice_junction{$coord}++;
	}
	elsif($type eq 'intronic/intergenic')
	{
		$count_intronic_intergenic{$coord}++;
	}
}
close(AMB);
print "[get-stats-SNPs.pl] $ambs lines written to $out_dir/VA_SNPs/ambiguous_SNPs.list\n";


my $cbs = 0;
for my $key (sort keys %count_codon_bias)
{
	print CB $key,"\t$count_codon_bias{$key}\n";
	$cbs ++;
}
close(CB);
print "[get-stats-SNPs.pl] $cbs lines written to $out_dir/VA_SNPs/codon_bias.stat\n";

print TYPE keys(%count_extension) . " extension SNPs\n";
print TYPE keys(%count_non_sense) . " non_sense SNPs\n";
print TYPE keys(%count_mis_sense) . " mis_sense SNPs\n";
print TYPE keys(%count_syn) . " synonymous SNPs\n";
print TYPE keys(%count_splice_junction) . " splice junction SNPs\n";
print TYPE keys(%count_intronic_intergenic) . " intronic\/intergenic SNPs\n";

my %chrom=();

print TYPE "\n";
print TYPE "Chromosomal Level Summary\n";
print TYPE "-------------------------\n";

for my $key (sort keys %count_extension)
{
	$key=~/(\S+)\:\d+/;
	$chrom{$1}++;
}

print TYPE "\n";
for my $key (sort keys %chrom)
{
	print TYPE "$key\textension\t$chrom{$key}\n";
}
print TYPE "\n";
my %chrom=();

for my $key (sort keys %count_non_sense)
{ 
        $key=~/(\S+)\:\d+/;
        $chrom{$1}++;
}

for my $key (sort keys %chrom)
{ 
        print TYPE "$key\tnon_sense\t$chrom{$key}\n";
}
print TYPE "\n";
my %chrom=();

for my $key (sort keys %count_mis_sense)
{ 
        $key=~/(\S+)\:\d+/;
        $chrom{$1}++;
}

for my $key (sort keys %chrom)
{ 
        print TYPE "$key\tmis_sense\t$chrom{$key}\n";
}
print TYPE "\n";
my %chrom=();

for my $key (sort keys %count_syn)  
{ 
        $key=~/(\S+)\:\d+/;
        $chrom{$1}++;
}

for my $key (sort keys %chrom)
{ 
        print TYPE "$key\tsynonymous\t$chrom{$key}\n";
}

print TYPE "\n";
my %chrom=();

for my $key (sort keys %count_splice_junction)
{
        $key=~/(\S+)\:\d+/;
        $chrom{$1}++;
}

for my $key (sort keys %chrom)
{
        print TYPE "$key\tsplice_junction\t$chrom{$key}\n";
}
print TYPE "\n";
my %chrom=();

for my $key (sort keys %count_intronic_intergenic)
{
        $key=~/(\S+)\:\d+/;
        $chrom{$1}++;
}

for my $key (sort keys %chrom)
{
        print TYPE "$key\tintronic\/intergenic\t$chrom{$key}\n";
}

close(TYPE);
print "[get-stats-SNPs.pl] Categorized SNP statistics written to $out_dir/VA_SNPs/categorized_SNPs.stat\n";

print "[get-stats-SNPs.pl] Done at ";
system(date);
