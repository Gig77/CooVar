#!/usr/bin/perl

use strict;
use Bio::Index::Fasta;
use Bio::Seq;
use Bio::SeqUtils;
use Bio::SeqIO;

sub get_pep
{
	my $cdna_seq = shift;
	my $seqOb = Bio::Seq->new(-seq => $cdna_seq);
        my @tr3frame = Bio::SeqUtils->translate_3frames($seqOb);
        my $pepSeq = $tr3frame[0];
        my $peptide_seq = $pepSeq->seq;
	return $peptide_seq;
}

sub get_cDNA
{
	my $exon = shift;
	my $seq = shift;
	my $strand = shift;
        my @line=split(/\t/,$exon);
        my $cdna_seq="";
	my $pseudo_fasta = "";
        my $start_CDS;
        my $end_CDS;

	my $sj="";

        for(my $i=0;$i<@line;$i++)
        {
		my @exons=split(/\-/,$line[$i]);
		
                my $dna = substr $seq, $exons[0]-1, abs ($exons[1] - $exons[0]+1);

                if($i==0)
                {
                	$start_CDS = $exons[0];
                }
                if($i == $#line)
                {
                	$end_CDS = $exons[1];
                }

                if($strand eq '+')
                {
                	$cdna_seq.= $dna;
			$pseudo_fasta.="\(" . $line[$i] . "\)\ " .  $dna . "\n";
                }
                else
                {
                	$dna = reverse($dna);
                        $dna=~tr/[ACTGactg]/[TGACtgac]/;
                        $cdna_seq= $dna . $cdna_seq;
			$pseudo_fasta="\($exons[1]\-$exons[0]\)\ " .  $dna . "\n" . $pseudo_fasta;
                } 

		#now lets get SJ info
		if($i<scalar(@line))
		{
                        if($i>0)
                        {
                                my $sj_seq = substr $seq, $exons[0]-3 , 2;
				my $label="";
                                if($strand eq '-')
                                {
                                        $sj_seq=~tr/[ACTGactg]/[TGACtgac]/;
                                        $sj_seq=reverse($sj_seq);
					$label = "donor";
                                }
				else
				{
					$label = "acceptor";
				}
                                my $start_sj = $exons[0]-2;
                                my $end_sj = $exons[0] -1;
                                $sj.="$start_sj\.\.$end_sj\;$sj_seq\;$strand\;$label ";
                        }


			if($i<scalar(@line) - 1)
			{
				my $sj_seq = substr $seq, $exons[1], 2;
				my $label = "";
				if($strand eq '-')
				{
					$sj_seq=~tr/[ACTGactg]/[TGACtgac]/;
					$sj_seq=reverse($sj_seq);
					$label = "acceptor";
				}
				else
				{
					$label = "donor";
				}
				my $start_sj=$exons[1] +1;
				my $end_sj = $exons[1] + 2;
				$sj.="$start_sj\.\.$end_sj\;$sj_seq\;$strand\;$label ";
			}
		}            
	}
	
		
	my $all = "$start_CDS\t$end_CDS\t$cdna_seq\t$pseudo_fasta\t$sj";

	return $all;
}

#open(D,$ARGV[0]) || die "$!\n"; #gff_file
my $fasta=$ARGV[1];
my $Index_File_Name = $fasta . 'tmp_index';
my $inx = Bio::Index::Fasta->new('-filename' => $Index_File_Name,'-write_flag' => 1);

$inx->make_index($fasta);

my %exon=();
my %exon_contig=();
my %exon_strand=();

system("sort -k 4 -n $ARGV[0] > sorted_gff3.tmp");

#gff3 file is sorted according to start 
open(D,"sorted_gff3.tmp") || die "$!\n";

open(SJ,">splice_junctions.tmp");

my %chrom=();
my %seen=();


while(<D>)
{
	chomp($_);
	next if ($_=~/^\#/);
	my @line=split(/\t/,$_);
	my @data=split(/\;/,$line[8]);
	my $id;
	for(my $i=0;$i<@data;$i++)
	{
		next if ($data[$i]!~/(ID\=|Parent\=|transcript_id\ )(\S+)/);
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

	$exon{$id}=$exon{$id} . $start . '-' . $end . "\t";
	$exon_contig{$id}=$line[0];
	$exon_strand{$id}=$line[6];
	
	next if (defined $seen{$id});
	$chrom{$line[0]}.=$id . " ";
	$seen{$id}++;
}
close(D);

open(ORI_CDNA,">reference_cDNA.fasta");
open(ORI_PEP,">reference_peptides.fasta");
open(ORI_PFA,">reference_cDNA.exons");
open(GFF3_TRANS,">transcripts.gff3.tmp");

print GFF3_TRANS "\#\#gff-version 3\n";

for my $chr (sort keys %chrom)
{
	my $seq = $inx->get_Seq_by_id($chr);

	print GFF3_TRANS "\#\#sequence\-region\ $chr\ 1\ " , length($seq->seq) , "\n";
	my @transcripts = split(/\s+/,$chrom{$chr});

	print STDERR scalar(@transcripts) , " transcripts in chromosome $chr\n";

	print STDERR "Starting $chr at ";
	system(date);

	for (my $i=0;$i<@transcripts;$i++)
	{
		next if ($exon_contig{$transcripts[$i]} ne $chr);
		my @line=split(/\t/,$exon{$transcripts[$i]});
		my $cdna=get_cDNA($exon{$transcripts[$i]},$seq->seq,$exon_strand{$transcripts[$i]});
		my @info = split(/\t/,$cdna);

		print ORI_CDNA '>',"$transcripts[$i]\t$exon_contig{$transcripts[$i]}\t$info[0]\t$info[1]\t$exon_strand{$transcripts[$i]}\n";
		print ORI_CDNA lc($info[2]) , "\n";
		print ORI_PEP '>',"$transcripts[$i]\t$exon_contig{$transcripts[$i]}\t$info[0]\t$info[1]\t$exon_strand{$transcripts[$i]}\n"; 
        	print ORI_PEP get_pep($info[2]), "\n";

		print ORI_PFA '>',"$transcripts[$i]\t$exon_contig{$transcripts[$i]}\t$info[0]\t$info[1]\t$exon_strand{$transcripts[$i]}\n";
		print ORI_PFA lc($info[3]);

		print GFF3_TRANS "$exon_contig{$transcripts[$i]}\tvariant_analyzer\tmRNA\t";
		print GFF3_TRANS "$info[0]\t$info[1]\t\.\t$exon_strand{$transcripts[$i]}\t\.\tID\=$transcripts[$i]\n";
		
		for(my $j=0;$j<@line;$j++)
		{
			print GFF3_TRANS "$exon_contig{$transcripts[$i]}\tvariant_analyzer\tcoding_exon\t";
			$line[$j]=~/(\d+)\-(\d+)/;
			my $start_ex = $1;
			my $end_ex = $2;
			print GFF3_TRANS "$start_ex\t$end_ex\t\.\t$exon_strand{$transcripts[$i]}\t\.\tParent\=$transcripts[$i]\n";
		}


		my @splice_junctions = split(/\s+/,$info[4]);

		for(my $j=0;$j<@splice_junctions;$j++)
		{
			$splice_junctions[$j]=~s/\;/\t/g;
			print SJ $transcripts[$i],"\t", "$chr\:$splice_junctions[$j]","\n";
		}
	}
	print STDERR "Done with $chr at ";
	system(date);
}
close(ORI_CDNA);
close(ORI_PEP);
close(ORI_PFA);
close(SJ);
close(GFF3_TRANS);
system("rm $Index_File_Name sorted_gff3\.tmp");
