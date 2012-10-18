#!/usr/bin/perl

use strict;

use Bio::SeqIO;
use POSIX qw(ceil floor);
use File::Basename;
use List::Util qw[min max];

print "[parse2circos.pl] Start executing script on ";
system("date");

sub print_circos
{
	my $coord_gv_ref = shift;
	my $chrom_ref = shift;
	my $bin_size = shift or die "[parse2circos.pl] print_circos(): bin_size not specified\n";
	my $file = shift or die "[parse2circos.pl] print_circos(): filename not specified\n";

	my %gvs=%$coord_gv_ref;
	my %chroms = %$chrom_ref;

	open(CIRCOS,">$file") || die "[parse2circos.pl] $!: $file\n";

	my $line = 0;
	for my $key (sort keys %chroms)
	{
       	for (my $bin = 0; $bin <= int($chroms{$key}/$bin_size); $bin++)
        {   
        	my $cand = "$key\:$bin";
            if(defined $gvs{$cand})
        	{          
	        	my $start = $bin*$bin_size+1;
    	    	my $end = ($bin+1)*$bin_size;
    	    	$end = $chroms{$key} if ($end > $chroms{$key});
				my $count = $gvs{$cand};
				print CIRCOS "$key\t$start\t$end\t$count\n";
				$line ++;
        	}
        }
	}
	close(CIRCOS);
	
	print "[parse2circos.pl] $line lines written to $file\n";
}

open(SNP,$ARGV[0]) || die "[parse2circos.pl] $!: $ARGV[0]\n";	#list of filtered SNPs
open(INS,$ARGV[1]) || die "[parse2circos.pl] $!: $ARGV[1]\n";	#list of filtered insertions
open(DEL,$ARGV[2]) || die "[parse2circos.pl] $!: $ARGV[2]\n";	#list of filtered deletions
open(EXON,$ARGV[3]) || die "[parse2circos.pl] $!: $ARGV[3]\n";	#list of exons
my $out_dir = $ARGV[5] or die "[parse2circos.pl] output directory not specified\n";

my $ref=Bio::SeqIO->new('-format'=>'Fasta','-file'=>$ARGV[4]);	#fasta sequence of genome

my $genome_length=0;
my %chrom_size=();

while (my $seq_obj=$ref->next_seq)
{
	$genome_length+=$seq_obj->length();
	$chrom_size{$seq_obj->id()}=$seq_obj->length();
}
my $bin_size = ceil(0.001*$genome_length);

my $snp_circos = $ARGV[0]  . '.circos';
my $ins_circos = $ARGV[1]  . '.circos';
my $del_circos = $ARGV[2]  . '.circos';
my $exon_circos= $ARGV[3]  . '.circos';
$exon_circos = "$out_dir/transcripts/".basename($exon_circos);

my %snp=();
while(<SNP>)
{
	chomp($_);
	my @line=split(/\t/,$_);
	my $chrom=$line[0];
	my $start=$line[1];
	my $bin=int($start/$bin_size);
	my $cand="$chrom\:$bin";
	$snp{$cand}++;
}  
close(SNP);

print_circos(\%snp,\%chrom_size,$bin_size,$snp_circos);

my %ins=();
my %dels=();

while(<INS>)
{
	chomp($_);
	my @line=split(/\t/,$_);
	my $chrom=$line[0];
	my $start=$line[1];
	my $end=$line[2];
	my $start_bin=int($start/$bin_size);
	my $end_bin=int($end/$bin_size);
	my $len_ins=length($line[3]);
		
	$ins{"$chrom\:$start_bin"}+=$len_ins;

	#insertion assoc with deletion
	if($end - $start > 1)
	{
		for(my $bin = $start_bin; $bin <= $end_bin; $bin++)
		{
			my $bin_start = $bin*$bin_size+1;
			my $bin_end = ($bin+1)*$bin_size;
			my $overlap = min($bin_end,$end-1)-max($bin_start,$start+1)+1;
			$dels{"$chrom\:$bin"} += $overlap;
		}
	}
}
close(INS);

print_circos(\%ins,\%chrom_size,$bin_size,$ins_circos);

while(<DEL>)
{
	chomp($_);
	my @line=split(/\t/,$_);
    my $chrom=$line[0];
    my $start=$line[1];
    my $end=$line[2];
	my $start_bin=int($start/$bin_size);
	my $end_bin=int($end/$bin_size);

	for(my $bin = $start_bin; $bin <= $end_bin; $bin++)
	{
		my $bin_start = $bin*$bin_size+1;
		my $bin_end = ($bin+1)*$bin_size;
		my $overlap = min($bin_end,$end)-max($bin_start,$start)+1;
		$dels{"$chrom\:$bin"} += $overlap;
	}
}
close(DEL);

print_circos(\%dels,\%chrom_size,$bin_size,$del_circos);

my %exons=();
while(<EXON>)      
{        
    chomp($_);   
	next if ($_=~/^\#/);
	next if ($_!~/\scds\s+\d+\s+\d+/); # 2012-10-04 | CF | ignore non-CDS entries   
    my @line=split(/\t/,$_);        
    my $chrom=$line[0];        
    my $start=$line[3];        
    my $end=$line[4];        
	my $start_bin=int($start/$bin_size);
	my $end_bin=int($end/$bin_size);

	for(my $bin = $start_bin; $bin <= $end_bin; $bin++)
	{
		my $bin_start = $bin*$bin_size+1;
		my $bin_end = ($bin+1)*$bin_size;
		my $overlap = min($bin_end,$end)-max($bin_start,$start)+1;
		$exons{"$chrom\:$bin"} += $overlap;
#		if ($chrom eq "I" and ($bin == 0 or $bin == 1))
#		{
#			print "$chrom\t$start\t$end\t$bin_start\t$bin_end\t$bin\t$overlap\n";
#		}
	}
}
close(EXON);

print_circos(\%exons,\%chrom_size,$bin_size,$exon_circos); 

print "[parse2circos.pl] Done at ";
system("date");
