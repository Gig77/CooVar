#!/usr/bin/perl

#------------------------------------------------------------------------------------------------
# DESCR: annotate GVs with genomic regions they overlap with
# INPUT: 1) list of GVs in GVF format (categorized-gvs.gvf output file, read from STDIN) 
#        2) genomic regions in GFF3 format (protein2genome.pl output, program parameter)  
# OUTPUT: GVs with annotated regions in GVF format (written to STDOUT)
#
# SYNOPSIS:
#
#   cat categorized-gvs.gvf | ./annotate-regions.pl mapped_protein_domains.gff3 > categorized-gvs.w_regions.gvf
#
# Example region input lines (produced by protein2genome.pl):
#
#   V	Motif_homol	interpro	7642716	7642826	.	-	.	ID=INTERPRO:IPR003582(CDS:F07C4.11);p_start=3;p_end=39
#   V	Motif_homol	interpro	7642542	7642670	.	-	.	ID=INTERPRO:IPR003582(CDS:F07C4.11);p_start=55;p_end=97
#
# These two entries specify that INTERPRO domain IPR003582 of protein-coding transcript CDS:4R79.1a maps
# to two genomic regions, one at chromosome V from 7642716 to 7642826, and one at a neighboring location on
# chromosome V from 7642542 to 7642670. The mapping of a protein domain to more than one genomic location
# is possible because a protein contains more than one copies of the domain (as shown in this case) or because
# a single domain spans a splice junction and is therefore split-mapped to multiple genomic locations. 
# Note that the format of the ID entry does not follow any special syntax and is just interpreted as 
# the name of the genomic region that is copied over into the annotated GVF output file. Also 'p_start' 
# and 'p_end' specifing start and end of this region in protein coordinates are optional.
#
# Example output line:
#
#   V	CooVar	SNV	7642800	7642800	.	+	.	ID=snp_245887;Variant_seq=G;Reference_seq=A;Variant_type=synonymous_codon;Variant_effect=synonymous_codon 0 mRNA CDS:F07C4.11;Note=tgt>tgC_C>C_aa11_codon_loc3;region=INTERPRO:IPR003582(CDS:F07C4.11)
#
# Note the "region=INTERPRO:IPR003582(CDS:F07C4.11)" annotation at the end of the line,
# which was not present before in the GVF input and specifies with which genomic regions this paritcular GV 
# (a SNP in this case) overlaps. If multiple regions are found to overlap with a GV
# they will be separated by a comma. 
#
# AUTHOR: Christian Frech
#------------------------------------------------------------------------------------------------

use strict;
use Set::IntervalTree;

print STDERR "[annotate-regions.pl] Start executing script on ".localtime()."\n";

my $region_file = $ARGV[0] or die "[annotate-regions.pl] input file with genomic regions not specified\n";

my $tree = Set::IntervalTree->new;
  
my $num_regions = 0;
open(REGION, $region_file) || die "[annotate-regions.pl] $!: $region_file\n";
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
close(REGION);

print STDERR "[annotate-regions.pl] $num_regions regions read from file $region_file\n";

# annotate GVs with region IDs
my $num_gvs = 0;
my $num_annotated_gvs = 0;
while(<STDIN>) 
{
	chomp;
	my $line = $_;

	if ($line =~ /^[#\s]/)
	{
		print "$line\n";
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

	print $line;
	print ";region=$regions" if ($regions);
	print "\n";
	
	$num_gvs ++;
	$num_annotated_gvs ++ if ($regions);
}

print STDERR "[annotate-regions.pl] $num_gvs GVs written.\n";
print STDERR "[annotate-regions.pl] $num_annotated_gvs GVs overlapped with input regions.\n";
	
print STDERR "[annotate-regions.pl] Done at ".localtime()."\n";
