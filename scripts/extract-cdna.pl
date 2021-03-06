use strict;

use Bio::Seq;
use Bio::SeqUtils;
use Bio::SeqIO;
use Bio::DB::Fasta;

print "[extract-cdna.pl] Start executing script on ";
print localtime()."\n";

my $out_dir = $ARGV[2] or die "[extract-cdna.pl] output directory not specified\n";
my $feature_type = $ARGV[3] or die "[extract-cdna.pl] feature type not specified\n"; 
my $feature_source = $ARGV[4]; # optional 

#my %invalidate_first_codon;
my $db = Bio::DB::Fasta->new($ARGV[1])
	or die ("[extract-cdna.pl] ERROR: Could not index/access FASTA file $ARGV[1].");

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
					or die ("[extract-cdna.pl] ERROR: Could not determine sequence for $chr:$exons[0]-$exons[1]\n");

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
									or die ("[extract-cdna.pl] ERROR: Could not determine sequence for $chr:".($exons[0]-2)."-".($exons[0]-1)."\n");
								
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
					or die ("[extract-cdna.pl] ERROR: Could not determine sequence for $chr:".($exons[1]+1)."-".($exons[1]+2)."\n");
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

sub filter_cds
{
	my $in = shift;
	my $out = shift;
	
	# filter input GFF/GTF file for CDS features
	my @cds;
	open(IN, "$in") or die ("[extract-cdna.pl] ERROR Could not open file $in\n");
	while(<IN>)
	{
		next if (/^\#/);
		if (/^([^\t]+)\t([^\t]+)\t($feature_type)\t(\d+)/i)
		{
			next if ($feature_source and $feature_source ne $2);
			push(@cds, [$1, $4, $_]);
		}
	}
	close(IN);

	die ("[extract-cdna.pl] ERROR: Could not find coding sequence (CDS) entries in file $in. Please make sure that this files contains features of type '$feature_type' ".($feature_source ? "and of source '$feature_source' " : "")."(see README).\n")
		if (@cds == 0);
	
	# sort entries
	my @sorted_cds = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @cds;
	
	# write filtered and sorted CDS to output file
	open(OUT, ">$out") or die ("[extract-cdna.pl] ERROR Could not write file $out\n");
	foreach my $e (@sorted_cds)
	{
		print OUT $e->[2];
	}
	close(OUT);
}

#open(D,$ARGV[0]) || die "$!\n"; #gff_file
my $fasta=$ARGV[1];

my %exon=();
my %exon_contig=();
my %exon_strand=();
my %exon_gene=();
my %exon_phase=();

# filter input GFF/GTF file for CDS features
print "[extract-cdna.pl] Filtering and sorting input GFF/GTF file $ARGV[0] ...\n";
filter_cds($ARGV[0], "$out_dir/intermediate-files/sorted_gff3.tmp");

#gff3 file is sorted according to start 
open(D,"$out_dir/intermediate-files/sorted_gff3.tmp") || die "[extract-cdna.pl] $!: $out_dir/intermediate-files/sorted_gff3.tmp\n";

open(SJ,">$out_dir/intermediate-files/splice_junctions.tmp") || die "[extract-cdna.pl] $!: $out_dir/intermediate-files/splice_junctions.tmp\n";

my %chrom=();
my %seen=();


print "[extract-cdna.pl] Parsing GFF file $out_dir/intermediate-files/sorted_gff3.tmp ...\n";
while(<D>)
{
	chomp($_);
	my @line=split(/\t/,$_);

	my ($id) = $line[8] =~ /[;\s]parent=([^;\s]+)/i; # check for parent ID first
	(my $tmp, $id) = $line[8] =~ /[;\s](ID=|transcript_id )([^;\s]+)/i if (!$id);
	$id=~s/\"//g;
	if (!$id)
	{
		print "[extract-cdna.pl] WARNING: could not parse transcript id from following line:\n$line[8]\n";
		next;		
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
#unlink("$out_dir/intermediate-files/sorted_gff3.tmp");

# 2012-10-09 | CF | extend transcript coordinate by three if last CDS is not a stop codon
# 2012-10-09 | CF | adjust for phase of first exon to ensure correct translation
print "[extract-cdna.pl] Ensuring inclusion of terminal stop codons and adjusting for phase of initial exons\n";
my ($adjusted_phase, $adjusted_stop) = (0, 0);
foreach my $k (keys(%exon))
{
	my @exons = split("\t", $exon{$k});
	my @phases = split("\t", $exon_phase{$k});
	
	die "[extract-cdna.pl] ERROR: Could not parse exons for transcript $k\n"
		if (@exons == 0);
	
	my $strand = $exon_strand{$k}; 
	die "[extract-cdna.pl] ERROR: Could not parse strand for transcript $k\n"
		if ($strand ne '+' and $strand ne '-');
		
	if ($exon_strand{$k} eq '+')
	{
		my $last = $exons[@exons-1];
		my ($start, $end) = split("-", $last);
		die "[extract-cdna.pl] ERROR: Could not parse exon coordinates for transcript $k: $last\n"
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
		die "[extract-cdna.pl] ERROR: Could not parse exon coordinates for transcript $k: $last\n"
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

print "[extract-cdna.pl]   Coordinates of $adjusted_phase transcripts were phase-adjusted to ensure their proper translation\n"
	if ($adjusted_phase > 0);
print "[extract-cdna.pl]   Coordinates of $adjusted_stop transcripts were adjusted to include terminal stop codons\n"
	if ($adjusted_stop > 0);

open(ORI_CDNA,">$out_dir/transcripts/reference_cdna.fasta") || die "[extract-cdna.pl] $!: $out_dir/transcripts/reference_cdna.fasta\n";;
open(ORI_PEP,">$out_dir/transcripts/reference_peptides.fasta") || die "[extract-cdna.pl] $!: $out_dir/transcripts/reference_peptides.fasta\n";;
open(ORI_PFA,">$out_dir/transcripts/reference_cdna.exons") || die "[extract-cdna.pl] $!: $out_dir/transcripts/reference_cdna.exons\n";;
open(GFF3_TRANS,">$out_dir/intermediate-files/transcripts.gff3.tmp") || die "[extract-cdna.pl] $!: $out_dir/intermediate-files/transcripts.gff3.tmp\n";;

print GFF3_TRANS "\#\#gff-version 3\n";

my ($cdnas, $peptides, $exons, $junctions) = (0, 0, 0, 0);
for my $chr (sort keys %chrom)
{
	print GFF3_TRANS "\#\#sequence\-region\ $chr\ 1\ ", $db->length($chr) , "\n";
	my @transcripts = split(/\s+/,$chrom{$chr});

	print "[extract-cdna.pl] Extracting ", scalar(@transcripts) , " transcripts on chromosome $chr on ";
	print localtime()."\n";

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
			$exons ++;
		}


		my @splice_junctions = split(/\s+/,$info[4]);

		for(my $j=0;$j<@splice_junctions;$j++)
		{
			$splice_junctions[$j]=~s/\;/\t/g;
			print SJ $transcripts[$i],"\t", "$chr\:$splice_junctions[$j]","\n";
			$junctions ++;
		}
	}
	print "[extract-cdna.pl]   Done with chromosome $chr on ";
	print localtime()."\n";
}
close(ORI_CDNA);
close(ORI_PEP);
close(ORI_PFA);
close(SJ);
close(GFF3_TRANS);

print "[extract-cdna.pl] $cdnas cDNAs written to $out_dir/transcripts/reference_cdna.fasta\n";
print "[extract-cdna.pl] $peptides peptides written to $out_dir/transcripts/reference_peptides.fasta\n";
print "[extract-cdna.pl] $exons exons written to $out_dir/transcripts/reference_cdna.exons\n";
print "[extract-cdna.pl] $junctions splice junctions written to $out_dir/intermediate-files/splice_junctions.tmp\n";

print "[extract-cdna.pl] Done at ";
print localtime()."\n";
