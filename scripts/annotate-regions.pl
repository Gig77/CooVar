#!/usr/bin/perl

use strict;
use Set::IntervalTree;

print "[annotate-regions.pl] Start executing script on ";
system("date");

my $categorized_GVs = $ARGV[0] or die "[annotate-regions.pl] input file with categorized GVs not specified\n";
my $region_file = $ARGV[1] or die "[annotate-regions.pl] input file with genomic regions not specified\n";
my $out_file = $ARGV[2] or die "[annotate-regions.pl] output file not specified\n";

open(GV, $categorized_GVs) || die "[annotate-regions.pl] $!: $categorized_GVs\n";
open(REGION, $region_file) || die "[annotate-regions.pl] $!: $region_file\n";
open(OUT, ">$out_file") || die "[annotate-regions.pl] $!: $out_file\n";

my $tree = Set::IntervalTree->new;
  
my $num_regions = 0;
while(<REGION>)
{
	chomp;
	next if (/^[#\s]/);
	
	my $line = $_;
	next if (!$line);
	
	my ($chr, $src, $type, $start, $end, $dummy1, $strand, $dummy2, $descr) = split("\t");
	my ($r_id) = $descr =~ /[\s;]?ID=([^;\n]+)/i;

	die("[annotate-regions.pl] Error parsing following region entry:\n$line\n")
		if (!$start or !$end or !$r_id);

	$end += 1 if ($start == $end);  # interval minimum size is 1
	
	# store argument in variable first, then pass on to C function call; 
	# otherwise arguments of function call get mixed up somehow and it will not work
	my $id = "$chr\t$r_id";  
	$tree->insert($id, $start, $end);
		
	$num_regions++;
}
print "[annotate-regions.pl] $num_regions regions read from file $region_file\n";

# annotate GVs with region IDs
my $num_gvs = 0;
my $num_annotated_gvs = 0;
while(<GV>) 
{
	chomp;
	my $line = $_;

	if ($line =~ /^[#\s]/)
	{
		print OUT "$line\n";
		next;
	}
			
	my ($chr, $src, $type, $start, $end) = split("\t", $line);

	die("[annotate-regions.pl] Error parsing following GV entry:\n$line\n")
		if (!$start or !$end);

	# get overlapping regions from interval tree
	my $regions = "";
	$end += 1 if ($start == $end);  # interval minimum size is 1
	my $overlaps = $tree->fetch($start, $end); 
	
	my %annotated;
	foreach my $region (@$overlaps)
	{
		my ($r_chr, $r_id) = split("\t", $region);
		next if ($chr ne $r_chr);  # ignore if overlap is on different chromosome
		next if ($annotated{$r_id});   # skip if this region has already been annotated for this GV
		
		$regions .= "," if ($regions);
		$regions .= $r_id;
		
		$annotated{$r_id} = 1;
	}

	print OUT $line;
	print OUT ";region=$regions" if ($regions);
	print OUT "\n";
	
	$num_gvs ++;
	$num_annotated_gvs ++ if ($regions);
}

print "[annotate-regions.pl] $num_gvs GVs written to $out_file\n";
print "[annotate-regions.pl] $num_annotated_gvs GVs overlapped with input regions.\n";
	
close(GV);
close(REGION);
close(OUT);

#print STDERR "[annotate-regions.pl] $overlapping_regions regions overlap with GVs.\n";

print "[annotate-regions.pl] Done at ";
system("date");
