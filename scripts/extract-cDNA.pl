use strict;

use Bio::Seq;
use Bio::SeqUtils;
use Bio::SeqIO;
use Bio::DB::Fasta;

print "[extract-cDNA.pl] Start executing script on ";
system("date");

my $out_dir = $ARGV[2] or die "[extract-cDNA.pl] output directory not specified\n";

#my %invalidate_first_codon;
my $db = Bio::DB::Fasta->new($ARGV[1])
	or die ("[extract-cDNA.pl] ERROR: Could not index/access FASTA file $ARGV[1].");

sub get_pep
{
	my $cdna_seq = shift;
	my $seqOb = Bio::Seq->new(-seq => $cdna_seq);
#	my @tr3frame = Bio::SeqUtils->translate_3frames($seqOb);
#	my $pepSeq = $tr3frame[0];
	my $pepSeq = $seqOb->translate;
	my $peptide_seq = $pepSeq->seq;
	return $peptide_seq;
}

sub get_cDNA
{
	my $exon = shift;
	my $chr = shift;
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

#               my $dna = substr $seq, $exons[0]-1, abs ($exons[1] - $exons[0]+1);
				my $dna = $db->seq($chr, $exons[0] => $exons[1])
					or die ("[extract-cDNA.pl] ERROR: Could not determine sequence for $chr:$exons[0]-$exons[1]\n");

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
 #               	$dna =~ s/^.../NNN/ if ($invalidate_first_codon{$line[$i]});
                	$cdna_seq.= $dna;
					$pseudo_fasta.="\(" . $line[$i] . "\)\ " .  $dna . "\n";
                }
                else
                {
#                	$dna =~ s/...$/NNN/ if ($invalidate_first_codon{$line[$i]});
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
#                               my $sj_seq = substr $seq, $exons[0]-3 , 2;
								my $sj_seq = $db->seq($chr, $exons[0]-2, $exons[0]-1)
									or die ("[extract-cDNA.pl] ERROR: Could not determine sequence for $chr:".($exons[0]-2)."-".($exons[0]-1)."\n");
								
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
#				my $sj_seq = substr $seq, $exons[1], 2;
				my $sj_seq = $db->seq($chr, $exons[1]+1, $exons[1]+2)
					or die ("[extract-cDNA.pl] ERROR: Could not determine sequence for $chr:".($exons[1]+1)."-".($exons[1]+2)."\n");
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

my %exon=();
my %exon_contig=();
my %exon_strand=();
my %exon_gene=();
my %exon_phase=();

# filter input GFF/GTF file for CDS features
print "[extract-cDNA.pl] Sorting input GFF file $ARGV[0] ...\n";
my $cmd = 'grep -iP "\scds\s+\d+\s+\d+"'." $ARGV[0] | sort -k 4 -n > $out_dir/sorted_gff3.tmp";
print "[extract-cDNA.pl]   Executing command: $cmd\n";
system($cmd) == 0 
  or die ("[extract-cDNA.pl] ERROR executing command: $cmd\n");

# make sure CDS features were found
my $line_count =  `wc -l < $out_dir/sorted_gff3.tmp`;
die ("[extract-cDNA.pl] ERROR: Could not find coding sequence (CDS) entries in file $ARGV[0]. Please make sure that this files contains features of type 'CDS' (see README).\n")
	if ($line_count == 0);

#gff3 file is sorted according to start 
open(D,"$out_dir/sorted_gff3.tmp") || die "[extract-cDNA.pl] $!: $out_dir/sorted_gff3.tmp\n";

open(SJ,">$out_dir/CV_Intermediate_Files/splice_junctions.tmp") || die "[extract-cDNA.pl] $!: $out_dir/CV_Intermediate_Files/splice_junctions.tmp\n";

my %chrom=();
my %seen=();


print "[extract-cDNA.pl] Parsing GFF file $out_dir/sorted_gff3.tmp ...\n";
while(<D>)
{
	chomp($_);
	next if ($_=~/^\#/);
	my @line=split(/\t/,$_);
	my @data=split(/\;/,$line[8]);
	
	if (scalar(@data) == 0)
	{
		print "[extract-cDNA.pl] WARNING: could not parse input line from file $ARGV[0]: $_\n";
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
	$exon{$id}=$exon{$id} . $start . '-' . $end . "\t";
	$exon_contig{$id}=$line[0];
	$exon_strand{$id}=$line[6];
	$exon_phase{$id}=$exon_phase{$id} . $line[7] . "\t";
	
	# remember exon/gene association (if specified)
	my ($geneid) = $line[8] =~ /gene=([^;]+)/i;
	($geneid) = $line[8] =~ /gene_id \"([^\"]+)\"/i if (!$geneid); 
	$exon_gene{$id}=$geneid;
	
	next if (defined $seen{$id});
	$chrom{$line[0]}.=$id . " ";
	$seen{$id}++;
}
close(D);
system("rm $out_dir/sorted_gff3.tmp");

# 2012-10-09 | CF | extend transcript coordinate by three if last CDS is not a stop codon
# 2012-10-09 | CF | adjust for phase of first exon to ensure correct translation
print "[extract-cDNA.pl] Ensuring inclusion of terminal stop codons and adjusting for phase of initial exons\n";
my ($adjusted_phase, $adjusted_stop) = (0, 0);
foreach my $k (keys(%exon))
{
	my @exons = split("\t", $exon{$k});
	my @phases = split("\t", $exon_phase{$k});
	
	die "[extract-cDNA.pl] ERROR: Could not parse exons for transcript $k\n"
		if (@exons == 0);
	
	my $strand = $exon_strand{$k}; 
	die "[extract-cDNA.pl] ERROR: Could not parse strand for transcript $k\n"
		if ($strand ne '+' and $strand ne '-');
		
	if ($exon_strand{$k} eq '+')
	{
		my $last = $exons[@exons-1];
		my ($start, $end) = split("-", $last);
		die "[extract-cDNA.pl] ERROR: Could not parse exon coordinates for transcript $k: $last\n"
			if (!$start or !$end);
		my $stop_cds = $db->seq($exon_contig{$k}, $end-2 => $end+3);
		if ($stop_cds !~ /(tag|taa|tga).../i and $stop_cds =~ /...(tag|taa|tga)/i)
		{
			$exons[@exons-1] = "$start-".($end+3); # extend by 3 to include stop
			$exon{$k} = join("\t", @exons);
			$adjusted_stop ++;
		}
		if ($phases[0] == 1 || $phases[0] == 2)
		{
#			print "Adjusting for phase for $k (+)...\n";
			($start, $end) = split("-", $exons[0]);
			$exons[0] = ($start-(3-$phases[0]))."-$end";
			$exon{$k} = join("\t", @exons);
#			$invalidate_first_codon{$exons[0]} = 1;
			$adjusted_phase ++;
		}
	}
	else
	{
		my $last = $exons[0];
		my ($start, $end) = split("-", $last);
		die "[extract-cDNA.pl] ERROR: Could not parse exon coordinates for transcript $k: $last\n"
			if (!$start or !$end);
		my $stop_cds = reverse($db->seq($exon_contig{$k}, $start-3, $start+2));
        $stop_cds=~tr/[ACTGactg]/[TGACtgac]/;
		if ($stop_cds !~ /(tag|taa|tga).../i and $stop_cds =~ /...(tag|taa|tga)/i)
		{
			$exons[0] = ($start-3)."-$end"; # extend by 3 to include stop
			$exon{$k} = join("\t", @exons);
			$adjusted_stop ++;
		}
		if ($phases[@phases-1] == 1 || $phases[@phases-1] == 2)
		{
#			print "Adjusting for phase for $k (-)...\n";
			($start, $end) = split("-", $exons[@exons-1]);
			$exons[@exons-1] = "$start-".($end+(3-$phases[@phases-1]));
			$exon{$k} = join("\t", @exons);
#			$invalidate_first_codon{$exons[@exons-1]} = 1;
			$adjusted_phase ++;
		}
	}
}

print "[extract-cDNA.pl]   Coordinates of $adjusted_phase transcripts were phase-adjusted to ensure their proper translation\n"
	if ($adjusted_phase > 0);
print "[extract-cDNA.pl]   Coordinates of $adjusted_stop transcripts were adjusted to include terminal stop codons\n"
	if ($adjusted_stop > 0);

open(ORI_CDNA,">$out_dir/CV_Transcripts/reference_cDNA.fasta") || die "[extract-cDNA.pl] $!: $out_dir/CV_Transcripts/reference_cDNA.fasta\n";;
open(ORI_PEP,">$out_dir/CV_Transcripts/reference_peptides.fasta") || die "[extract-cDNA.pl] $!: $out_dir/CV_Transcripts/reference_peptides.fasta\n";;
open(ORI_PFA,">$out_dir/CV_Transcripts/reference_cDNA.exons") || die "[extract-cDNA.pl] $!: $out_dir/CV_Transcripts/reference_cDNA.exons\n";;
open(GFF3_TRANS,">$out_dir/CV_Intermediate_Files/transcripts.gff3.tmp") || die "[extract-cDNA.pl] $!: $out_dir/CV_Intermediate_Files/transcripts.gff3.tmp\n";;

print GFF3_TRANS "\#\#gff-version 3\n";

my ($cdnas, $peptides, $exons, $junctions) = (0, 0, 0, 0);
for my $chr (sort keys %chrom)
{
	print GFF3_TRANS "\#\#sequence\-region\ $chr\ 1\ ", $db->length($chr) , "\n";
	my @transcripts = split(/\s+/,$chrom{$chr});

	print "[extract-cDNA.pl] Extracting ", scalar(@transcripts) , " transcripts on chromosome $chr on ";
	system("date");

	for (my $i=0;$i<@transcripts;$i++)
	{
		next if ($exon_contig{$transcripts[$i]} ne $chr);
		my @line=split(/\t/,$exon{$transcripts[$i]});
		my $cdna=get_cDNA($exon{$transcripts[$i]},$chr,$exon_strand{$transcripts[$i]});
		my @info = split(/\t/,$cdna);

		print ORI_CDNA '>',"$transcripts[$i]\t$exon_contig{$transcripts[$i]}\t$info[0]\t$info[1]\t$exon_strand{$transcripts[$i]}\n";
		print ORI_CDNA lc($info[2]) , "\n";
		$cdnas ++;
		print ORI_PEP '>',"$transcripts[$i]\t$exon_contig{$transcripts[$i]}\t$info[0]\t$info[1]\t$exon_strand{$transcripts[$i]}\n"; 
        print ORI_PEP get_pep($info[2]), "\n";
        $peptides ++;

		print ORI_PFA '>',"$transcripts[$i]\t$exon_contig{$transcripts[$i]}\t$info[0]\t$info[1]\t$exon_strand{$transcripts[$i]}\n";
		print ORI_PFA lc($info[3]);
		$exons ++;

		print GFF3_TRANS "$exon_contig{$transcripts[$i]}\tCooVar\tmRNA\t";
		print GFF3_TRANS "$info[0]\t$info[1]\t\.\t$exon_strand{$transcripts[$i]}\t\.\tID\=$transcripts[$i]";
		print GFF3_TRANS ";gene=$exon_gene{$transcripts[$i]}" if (defined $exon_gene{$transcripts[$i]});
		print GFF3_TRANS "\n";
		
		for(my $j=0;$j<@line;$j++)
		{
			print GFF3_TRANS "$exon_contig{$transcripts[$i]}\tCooVar\tCDS\t";
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
			$junctions ++;
		}
	}
	print "[extract-cDNA.pl]   Done with chromosome $chr on ";
	system("date");
}
close(ORI_CDNA);
close(ORI_PEP);
close(ORI_PFA);
close(SJ);
close(GFF3_TRANS);

print "[extract-cDNA.pl] $cdnas cDNAs written to $out_dir/CV_Transcripts/reference_cDNA.fasta\n";
print "[extract-cDNA.pl] $peptides peptides written to $out_dir/CV_Transcripts/reference_peptides.fasta\n";
print "[extract-cDNA.pl] $exons exons written to $out_dir/CV_Transcripts/reference_cDNA.exons\n";
print "[extract-cDNA.pl] $junctions splice junctions written to $out_dir/CV_Intermediate_Files/splice_junctions.tmp\n";

print "[extract-cDNA.pl] Done at ";
system("date");
