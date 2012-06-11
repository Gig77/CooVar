#!/usr/bin/perl

use strict vars;
use strict subs;
use Bio::SeqIO;
use POSIX qw(ceil floor);

sub print_circos
{
	my $coord_gv_ref = shift;
	my $chrom_ref = shift;
	my $window_pace = shift;
	my $label = shift;

	my %gvs=%$coord_gv_ref;
	my %chroms = %$chrom_ref;

	for my $key  (sort keys %chroms)
	{

		my $window_count=0;
        	my $start = 1;
        	for (my $i=1;$i<=$chroms{$key};$i++)
        	{   
                	my $cand="$key\:$i";

	                if(defined $gvs{$cand})
        	        {          
				$window_count+=$gvs{$cand};
                	}        

	                if($i%$window_pace == 0)
        	        {
                	        my $end=$start+$window_pace-1;
                	        if($end>$chroms{$key})
                        	{
                                	$end=$chroms{$key};
                        	}

	                        if($window_count > 0)
        	                {
                	                print {$label} "$key\t$start\t$end\t$window_count\n";
                        	}
                        	$window_count=0;  
                        	$start+=$window_pace;   
                	}
        	}
        	my $end=$start+$window_pace-1;
        	if($end>$chroms{$key})
        	{
                	$end=$chroms{$key};
        	}

	        if($window_count > 0)
        	{
                	print {$label} "$key\t$start\t$end\t$window_count\n";
        	}
        	$start+=$window_pace;
	}
	return;
}


open(SNP,$ARGV[0]) || die "$!\n";	#list of filtered SNPs
open(INS,$ARGV[1]) || die "$!\n";	#list of filtered insertions
open(DEL,$ARGV[2]) || die "$!\n";	#list of filtered deletions
open(EXON,$ARGV[3]) || die "$!\n";	#list of exons

my $ref=Bio::SeqIO->new('-format'=>'Fasta','-file'=>$ARGV[4]);	#fasta sequence of genome

my $genome_length=0;
my %chrom_size=();

while (my $seq_obj=$ref->next_seq)
{
        $genome_length+=$seq_obj->length();
	$chrom_size{$seq_obj->id()}=$seq_obj->length();
}
my $window_size = ceil(0.001*$genome_length);

my $snp_circos = $ARGV[0]  . '.circos';
my $ins_circos = $ARGV[1]  . '.circos';
my $del_circos = $ARGV[2]  . '.circos';
my $exon_circos= $ARGV[3]  . '.circos';

open(CIRCOS_SNP,">$snp_circos");
open(CIRCOS_INS,">$ins_circos");       
open(CIRCOS_DEL,">$del_circos");         
open(CIRCOS_EX,">$exon_circos");

my %snp=();

while(<SNP>)
{
	chomp($_);
	my @line=split(/\t/,$_);
	my $chrom=$line[0];
	my $start=$line[1];
	my $cand="$chrom\:$start";
	$snp{$cand}=1;
}  
close(SNP);

print_circos(\%snp,\%chrom_size,$window_size,"CIRCOS_SNP");

for my $key (keys %snp)
{
	delete($snp{$key});
}

my %ins=();
my %dels=();

while(<INS>)
{
	chomp($_);
	my @line=split(/\t/,$_);
	my $chrom=$line[0];
	my $start=$line[1];
	my $end=$line[2];
	my $len_ins = length($line[3]);
	$ins{"$chrom\:$start"}+=$len_ins;
	#insertion assoc with deletion
	if($end - $start > 1)
	{
	        for (my $i=$start+1;$i<=$end-1;$i++)        
        	{        
                	$dels{"$chrom\:$i"}=1;     
        	}
	}
}
close(INS);

print_circos(\%ins,\%chrom_size,$window_size,"CIRCOS_INS");

for my $key (keys %ins) 
{ 
        delete($ins{$key});
}

while(<DEL>)
{
	chomp($_);
	my @line=split(/\t/,$_);
        my $chrom=$line[0];
        my $start=$line[1];
        my $end=$line[2];
        for (my $i=$start;$i<=$end;$i++)
        {
		$dels{"$chrom\:$i"}=1;
        }
}
close(DEL);

print_circos(\%dels,\%chrom_size,$window_size,"CIRCOS_DEL");

for my $key (keys %dels)
{ 
        delete($dels{$key});
}

my %exons=();
while(<EXON>)      
{        
        chomp($_);   
	next if ($_=~/^\#/);     
        my @line=split(/\t/,$_);        
        my $chrom=$line[0];        
        my $start=$line[3];        
        my $end=$line[4];        
        for (my $i=$start;$i<=$end;$i++)        
        {        
                $exons{"$chrom\:$i"}=1;    
        }
}
close(EXON);

print_circos(\%exons,\%chrom_size,$window_size,"CIRCOS_EX"); 

for my $key (keys %exons)
{            
        delete($exons{$key});
}    

close(CIRCOS_SNP);
close(CIRCOS_INS);
close(CIRCOS_DEL);
close(CIRCOS_EX);
