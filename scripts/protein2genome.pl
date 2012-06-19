#!/usr/bin/perl

#------------------------------------------------------------------------------------------------
# DESCR: map protein coordinates back to genome
# INPUT: 1) list of transcript IDs and protein coordinates (read from STDIN) 
#        2) GFF file with transcript CDS definitions  
# OUTPUT: input lines annotated with genomic coordinates (written to STDOUT)
#
# SYNOPSIS:
#
#     $ echo "4R79.1a 10 30 protein domain 1" | perl protein2genome.pl cds.gff
#     4R79.1a     10      30      IV      17490363        17490425        protein domain 1
#
#     Each line in the output consists of the first three columns from the input 
#     line, followed by the mapped chromosome name, genomic start coordinate, 
#     and genomic end coordinate. All remaining columns are copied from the input line. 
#------------------------------------------------------------------------------------------------

use strict;

my %transcript_exons;
my %transcript_strand;
my %transcript_contig;
#my $exons = [ [100,200] , [300,400] , [500,600] ];
#my $tr = "CDS:AC3.1";
#my $tr = "CDS:4R79.1a";

#my $exon_file = "/home/cfa24/scripts_data/hawaiian/variants/ALL_VARIANTS/VARIANT_ANALYZER_INPUT/exons_WS210_nuclear_DNA.gff3";
my $gff_file = $ARGV[0] or die "GFF file with CDS definitions not specified\n";
load_cds($gff_file);

# read protein coordinates from STDIN
while (<STDIN>)
{
	chomp;
	my $line = $_;
	my ($tr, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = split("\t", $line);
	
	if (!$tr or !$start or !$end or !$attributes)
	{
		warn "[protein2genome.pl] WARNING: Could not parse input line:\n$line\n";
		next;
	}
	
	my ($id) = $attributes =~ /^ID=([^;\s]+)/i;
	if (!$id)
	{
		warn "[protein2genome.pl] WARNING: Could not parse ID from input line:\n$line\n";
		next;	
	}
	$attributes =~ s/^ID=[^;]+;?//i;

	if (!exists $transcript_contig{$tr})
	{
		warn "[protein2genome.pl] WARNING: Transcript ID '$tr' not found in GFF file $gff_file\n";
		next;
	}
		
	my ($nt_start, $nt_end) = aa2nt([$start, $end], $transcript_exons{$tr}, $transcript_strand{$tr});
	if ($nt_start < 0 or $nt_end < 0)
	{
		warn "[protein2genome.pl] WARNING: Could not map aa coordinate $start for following input line:\n" if ($nt_start < 0);
		warn "[protein2genome.pl] WARNING: Could not map aa coordinate $end for following input line:\n" if ($nt_end < 0);
		warn "$line\n";
		next;
	}
	
	($nt_start, $nt_end) = ($nt_end, $nt_start) if ($nt_start > $nt_end);

	print "$transcript_contig{$tr}\t$source\t$feature\t$nt_start\t$nt_end\t$score\t$transcript_strand{$tr}\t$frame\tID=$id($tr);";
	print "$attributes;" if ($attributes);
	print "p_start=$start;p_end=$end\n";
}

sub load_cds
{
	my $gff = shift;
	open(D,$gff) or die "could not read file $gff\n";
	while(<D>)
	{
		chomp($_);
		next if ($_=~/^\#/);
		my @line=split(/\t/,$_);
		my @data=split(/\;/,$line[8]);
		
		if (scalar(@data) == 0)
		{
			warn "[protein2genome.pl] WARNING: Could not parse GFF input line: $_\n";
			next;		
		}
		
		my $id;
		for(my $i=0;$i<@data;$i++)
		{
			next if ($data[$i]!~/(ID\=|Parent\=|transcript_id\ )([^;\s]+)/);
			$id = $2;
			$id=~s/\"//g;
			last;
		}
		my ($start,$end);
		
		if($line[3] <= $line[4])
		{
			$start = $line[3];
			$end = $line[4];
		}
		else
		{
			$start = $line[4];
			$end = $line[3];
		}

		$transcript_contig{$id} = $line[0]; 
		$transcript_strand{$id} = ($line[6] eq '+' or $line[6] == 1) ? '+' : '-'; 
		if (!exists $transcript_exons{$id})
		{
			$transcript_exons{$id} = [[$start, $end]];
		}
		else
		{
			push(@{$transcript_exons{$id}}, [$start, $end]);
		}		
	}
	close(D);
	
}

sub aa2nt 
{
	my $aa_coordinates = shift;  # scalar or array-ref with amino acid coordinate(s) (1-based)
	my $exons = shift;           # array-ref of exons, each representing an array-ref of a start/end coordinate pair (e.g. [ [100,200] , [300,400] , [500,600] ])
	my $strand = shift;          # strand; 1 or '+' --> plus strand; else negative strand
	my $codon_pos = shift;       # codon position (optional); map to 'start' or 'end' of codon; default is 'start'; 
								 #   if two input coordinates are specified, by default the first coordinate is mapped 
								 #   to the codon start position and the second coordinate is mapped to the codon end position

	die("[protein2genome.pl] FATAL: aa coordinate not specified") if (!defined $aa_coordinates or ref($aa_coordinates) == 'ARRAY' and @$aa_coordinates == 0);
	die("[protein2genome.pl] FATAL: exons not array reference") if (!defined $exons or ref($exons) ne "ARRAY");
	die("[protein2genome.pl] FATAL: empty exon array") if (@$exons == 0);
	die("[protein2genome.pl] FATAL: strand not specified") if (!defined $strand);
	die("[protein2genome.pl] FATAL: invalid codon position") if (defined $codon_pos and $codon_pos ne 'start' and $codon_pos ne 'end');
	
	# create hash map with all possible protein coordinates mapped to nucleotide coordinates
	my %map;
	if ( $strand > 0 or $strand eq '+' ) 
	{
		my ( $ppos, $residue ) = ( 1, 0 );
		my @sorted_exons = sort { $a->[0] <=> $b->[0] } (@$exons);
		foreach my $e (@sorted_exons) {
			$map{ $ppos - 1 }{'end'} = $e->[0] + $residue - 1 if ($residue);
			for (my $ntpos = $e->[0] + $residue; $ntpos <= $e->[1]; $ntpos += 3) {
				$map{$ppos}{'start'} = $ntpos;
				$map{$ppos}{'end'}   = $ntpos + 2;
				$ppos++;
				$residue = $ntpos + 2 - $e->[1];
			}
		}
	}
	else 
	{
		my ( $ppos, $residue ) = ( 1, 0 );
		my @sorted_exons = reverse sort { $a->[0] <=> $b->[0] } (@$exons);
		foreach my $e (@sorted_exons) {
			$map{ $ppos - 1 }{'end'} = $e->[1] - $residue + 1 if ($residue);
			for (my $ntpos = $e->[1] - $residue; $ntpos >= $e->[0]; $ntpos -= 3)
			{
				$map{$ppos}{'start'} = $ntpos;
				$map{$ppos}{'end'}   = $ntpos - 2;
				$ppos++;
				$residue = $e->[0] - ($ntpos - 2);
			}
		}
	}
	
	# lookup mapping, return mapped coordinates
	my (@aa_coords, @nt_coords);
	push(@aa_coords, ref($aa_coordinates) eq 'ARRAY' ? @$aa_coordinates : $aa_coordinates);
	for(my $i = 0; $i < @aa_coords; $i ++)
	{
		die("[protein2genome.pl] FATAL: input aa coordinate undefined\n")
			if (!$aa_coords[$i]);
			
		my $nt;
		if ($codon_pos) {
			$nt = $map{$aa_coords[$i]}{$codon_pos};
		}
		elsif (@aa_coords == 2 and $i == 0)	{
			$nt = $map{$aa_coords[$i]}{'start'};
		}
		elsif(@aa_coords == 2 and $i == 1) {
			$nt = $map{$aa_coords[$i]}{'end'};
		}
		else {
			$nt = $map{$aa_coords[$i]}{'start'};			
		}

		$nt = -1 if (!$nt);
		push(@nt_coords, $nt);
	}

	return (wantarray) ? @nt_coords : $nt_coords[0];
}
