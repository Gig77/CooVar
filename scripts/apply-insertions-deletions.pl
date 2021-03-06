#!/usr/bin/perl

use strict;
use Bio::Seq;
use Bio::SeqUtils;
use Bio::SeqIO;
use POSIX qw(ceil floor);
use Set::IntervalTree;

print "[apply-insertions-deletions.pl] Start executing script on ";
print localtime()."\n";

my %ref_cdnas=();

sub evaluate_cDNA
{
	my $seq = shift;
	my $name = shift;
	my @result;

	#print STDERR "Checking $name with sequence $seq\n";

	if(length($seq) < 3)
    {
		$result[0]="";
		$result[1]="FULLY_DELETED";
		$result[2]="N\.A\.";
		return @result;
    }

    my $seqOb = Bio::Seq->new(-seq => $seq);
    my @tr3frame = Bio::SeqUtils->translate_3frames($seqOb);
    my $pepSeq = $tr3frame[0];
    my $peptide_seq = $pepSeq->seq;

	$result[0]=$peptide_seq;
	my $peptide_length = length($peptide_seq);

    if($peptide_seq=~/^([^\*]*)\*\S+$/) # internal stop codon in variant peptide sequence?
    {
		my $pre = $1;
		my $pre_length = length($1);
	
#		# 2012-12-12 | CF | check also reference peptide sequence for internal stops
		my $ref_cDNA = $ref_cdnas{$name}
			or die "[apply-insertions-deletions.pl] ERROR: Could not determine reference cDNA for transcript $name\n";
		my $refSeqOb = Bio::Seq->new(-seq => $ref_cDNA);
		my $refPepSeq = $refSeqOb->translate->seq;
#		print "$ref_cDNA\n";
#		print "REF: $refPepSeq\n";
#		print "ALT: $peptide_seq\n";
		if ($refPepSeq =~ /^([^\*]*)\*\S+$/)
		{
			if (length($1) <= $pre_length)
			{
				$result[1]="ORF_PRESERVED";
				$result[2]="N\.A\.";
				return @result;				
			}
		}

		my $perc = 100*($pre_length+1)/$peptide_length;
		if($perc=~/(\d+\.\d\d)/)
		{
			$perc = $1;
		}
        $result[1]= "ORF_DISRUPTED";
		$result[2]=$perc;
	}
	else
	{
		$result[1]="ORF_PRESERVED";
		$result[2]="N\.A\.";
	}
	
	return @result;
}

# hash of Set::InvervalTree, one per chromosome, to efficiently query for overlaps between exons and GVs
my (%del_tree, %ins_tree); 

sub get_overlap_del
{
	my $chr = shift;
	my $transcript_exons = shift;

    my $overlapping_GVs = "";

	my @exons = split(/\n/,$transcript_exons);
	for (my $i=0;$i<@exons;$i++)
	{
		$exons[$i]=~/\((\d+)\-(\d+)\)\s+/;
		my $exon_start = $1;
		my $exon_end = $2;

		($exon_start, $exon_end) = ($exon_end, $exon_start) 
			if ($exon_start > $exon_end); # happens for transcripts in the - strand
		$exon_end ++;  # to find GVs overlapping EXACTLY with exon end coordinate

		next if (!$del_tree{$chr});
		my $overlaps = $del_tree{$chr}->fetch($exon_start, $exon_end); # WARNING: do not include ANY expression in parameter list; confuses C function call 
		foreach my $overlap (@$overlaps)
		{
#			print "OVERLAP: $del_tree{$chr};chr=$chr exon_start=$exon_start exon_end=$exon_end overlap=$overlap\n";
			my ($gv_start, $gv_end) = split("\t", $overlap);
	        $overlapping_GVs .= "$gv_start\.\.$gv_end ";
		}
	}
    
    return $overlapping_GVs;
}

sub get_overlap_ins
{
	my $chr = shift;
	my $transcript_exons = shift;

	my $overlapping_GVs = "";

	my @exons = split(/\n/,$transcript_exons);
	for (my $i=0;$i<@exons;$i++)
	{
		$exons[$i]=~/\((\d+)\-(\d+)\)\s+/;
		my $exon_start = $1;
		my $exon_end = $2;

		($exon_start, $exon_end) = ($exon_end, $exon_start) 
			if ($exon_start > $exon_end); # happens for transcripts in the - strand
#		$exon_end ++;  # to find GVs overlapping EXACTLY with exon end coordinate
		
		next if (!$ins_tree{$chr});
		my $overlaps = $ins_tree{$chr}->fetch($exon_start, $exon_end); # WARNING: do not include ANY expression in parameter list; confuses C function call 
		foreach my $overlap (@$overlaps)
		{
			my ($gv_start, $gv_end) = split("\t", $overlap);
	        $overlapping_GVs .= "$gv_start ";
		}
	}

	return $overlapping_GVs;
}

my $ref_pseudo_fasta = $ARGV[0];
my $pseudo_fasta = $ARGV[1];
my $categorized_snps = $ARGV[2];
my $insertion_list = $ARGV[3];
my $deletion_list = $ARGV[4];
my $gvf_gvs = $ARGV[5];
my $trans_gff3 = $ARGV[6];
my $sj_file = $ARGV[7];
my $out_dir = $ARGV[8] or die "[apply-insertions-deletions.pl] output directory not specified\n";

open(cDNA_ref,$ref_pseudo_fasta) || die "[apply-insertions-deletions.pl] $!: $ref_pseudo_fasta\n";
open(cDNA_SNP,$pseudo_fasta) || die "[apply-insertions-deletions.pl] $!: $pseudo_fasta\n";
open(SNP,$categorized_snps) || die "[apply-insertions-deletions.pl] $!: $categorized_snps\n";
open(INS,$insertion_list) || die "[apply-insertions-deletions.pl] $!: $insertion_list\n";
open(DEL,$deletion_list) || die "[apply-insertions-deletions.pl] $!: $deletion_list\n";
open(SJ,$sj_file) || die "[apply-insertions-deletions.pl] $!: $sj_file\n";
open(GVF,">>$gvf_gvs") || die "[apply-insertions-deletions.pl] $!: $gvf_gvs\n";
open(TRANS_GFF3,$trans_gff3) || die "[apply-insertions-deletions.pl] $!: $trans_gff3\n";

my %ref_exons=();
my $ref_id;

print "[apply-insertions-deletions.pl] Reading file $ref_pseudo_fasta ...\n";
while(<cDNA_ref>)
{
        chomp($_);
        if($_=~/^\>(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
        {
                $ref_id=$1;
                my $chrom=$2;
                my $start = $3;
                my $end = $4;
                my $strand = $5;
        }
	else
	{
                $ref_exons{$ref_id}.=$_ . "\n";
                /\s(.*)/;
                $ref_cdnas{$ref_id} .= $1;
    }
}
close(cDNA_ref);

print "[apply-insertions-deletions.pl] Reading splice junctions from file $sj_file ...\n";
my %sj=();
while(<SJ>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	$sj{$line[0]} .= $line[1] . "\;$line[2]\;$line[4] ";
}
close(SJ);


my %trans_exons= ();
my %trans_start= ();
my %trans_end  = ();
my %trans_chrom= ();
my %trans_strand=();

my $id;

my %chrom_trans=();

print "[apply-insertions-deletions.pl] Reading file $pseudo_fasta ...\n";
while(<cDNA_SNP>)
{
        chomp($_);
        if($_=~/^\>(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
        {
                $id=$1;
                my $chrom=$2;
                my $start = $3;
                my $end = $4;
                my $strand = $5;
                $trans_chrom{$id} = $chrom;
                $trans_start{$id} = $start;
                $trans_end{$id} = $end;
                $trans_strand{$id} = $strand;
		$chrom_trans{$chrom} = [] if (!exists $chrom_trans{$chrom});
		push(@{$chrom_trans{$chrom}}, $id);
#		$chrom_trans{$chrom}.=$id . " ";
        }
        else
        {
                $trans_exons{$id}.=$_ . "\n";
        }
}
close(cDNA_SNP);

my %trans2snp=();
my %trans2splice=();

print "[apply-insertions-deletions.pl] Reading file $categorized_snps ...\n";
while(<SNP>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	my $snp="";
	$line[3]=~/(\S)(\S)(\S)\>(\S)(\S)(\S)/;
	my $ref1=$1;
	my $ref2=$2;
	my $ref3=$3;

	my $tar1=$4;
	my $tar2=$5;
	my $tar3=$6;

	if($line[7] == 1)
	{
		$snp="$ref1\>$tar1";
	}
	elsif($line[7] == 2)
	{
		$snp="$ref2\>$tar2";
	}
	elsif($line[7] == 3)
	{
		$snp="$ref3\>$tar3";
	}
	$trans2snp{$line[1]}.=$line[0] . "\;" .  $snp . " ";
	$trans2splice{$line[1]}.=$line[0] . "\;" .  $snp . " " if ($line[2] =~ /splice_junction/); # 2012-06-19 | CF |
}
close(SNP);

my %to_ins=();
#my %to_del=();
#my %range_del=();
my %chrom_ins=();
my %chrom_del=();


my (@dist_ins,@dist_del,@dist_ins_exons,@dist_del_exons);

print "[apply-insertions-deletions.pl] Reading file $insertion_list ...\n";
while(<INS>)
{
	chomp($_);
	my @line = split(/\t/,$_);
	$to_ins{"$line[0]\:$line[1]"}="$line[2]\-$line[3]";

	$dist_ins[length($line[3])]++;

	if(abs($line[2] - $line[1]) > 1)
	{
		my $start_del = $line[1] +1;
		my $end_del = $line[2] - 1;
#		$range_del{"$line[0]\:$start_del\.\.$end_del"}++;


		my $len_del = $end_del - $start_del + 1;
		$dist_del[$len_del]++;

#		for (my $i=$start_del;$i<=$end_del;$i++)
#		{
#			$to_del{"$line[0]\:$i"}++;
#		}
		
		# interval must be at least 1bp, otherwise overlaps with exon start coordinate will be missed
		$del_tree{$line[0]} = Set::IntervalTree->new if (!$del_tree{$line[0]});
		my ($id, $start, $end) = ("$start_del\t$end_del", $start_del, $end_del + 1);
		$del_tree{$line[0]}->insert($id, $start, $end); # WARNING: do not include ANY expression in parameter list; confuses C function call and strange stuff is inserted into tree!  

		$chrom_del{$line[0]} = [] if (!exists $chrom_del{$line[0]});
		push(@{$chrom_del{$line[0]}}, "$start_del\.\.$end_del");
	}

	# interval must be at least 1bp, otherwise overlaps with exon start coordinate will be missed
 	$ins_tree{$line[0]} = Set::IntervalTree->new if (!$ins_tree{$line[0]});
	my ($id, $start, $end) = ("$line[1]\t$line[2]", $line[1], $line[1] + 1);
	$ins_tree{$line[0]}->insert($id, $start, $end); # WARNING: do not include ANY expression in parameter list; confuses C function call and strange stuff is inserted into tree!
	
	$chrom_ins{$line[0]} = [] if (!exists $chrom_ins{$line[0]});
	push(@{$chrom_ins{$line[0]}}, "$line[1]\.\.$line[2]\-$line[3]");
}
close(INS);

print "[apply-insertions-deletions.pl] Reading file $deletion_list ...\n";
while(<DEL>)
{
	chomp($_);
	my @line = split(/\t/,$_);
        
	my $len_del = $line[2] - $line[1] + 1;
        $dist_del[$len_del]++;

	# interval must be at least 1bp, otherwise overlaps with exon start coordinate will be missed
	$del_tree{$line[0]} = Set::IntervalTree->new if (!$del_tree{$line[0]});
	my ($id, $start, $end) = ("$line[1]\t$line[2]", $line[1], $line[2] + 1);
	$del_tree{$line[0]}->insert($id, $start, $end); # WARNING: do not include ANY expression in parameter list; confuses C function call and strange stuff is inserted into tree! 

#	$range_del{"$line[0]\:$line[1]\.\.$line[2]"}++;
#	for(my $i=$line[1];$i<=$line[2];$i++)
#	{
#		$to_del{"$line[0]\:$i"}++;
#	}

	$chrom_del{$line[0]} = [] if (!exists $chrom_del{$line[0]});
	push(@{$chrom_del{$line[0]}}, "$line[1]\.\.$line[2]");
}
close(DEL);

open(MOD_PFA,">$out_dir/transcripts/variant_cdna.exons") || die "[apply-insertions-deletions.pl] $!: $out_dir/transcripts/variant_cdna.exons\n";;
open(MOD_cDNA,">$out_dir/transcripts/variant_cdna.fasta") || die "[apply-insertions-deletions.pl] $!: $out_dir/transcripts/variant_cdna.fasta\n";
open(MOD_PEP,">$out_dir/transcripts/variant_peptides.fasta") || die "[apply-insertions-deletions.pl] $!: $out_dir/transcripts/variant_peptides.fasta\n";
open(GENE_SUMMARY,">$out_dir/intermediate-files/variant.summary") || die "[apply-insertions-deletions.pl] $!: $out_dir/intermediate-files/variant.summary\n";
open(STAT,">$out_dir/transcripts/variant.stat") || die "[apply-insertions-deletions.pl] $!: $out_dir/transcripts/variant.stat\n";
open(ALIGNMENT,">$out_dir/transcripts/reference_variant.alignment") || die "[apply-insertions-deletions.pl] $!: $out_dir/transcripts/reference_variant.alignment\n";
open(DIST_INS,">$out_dir/insertions/distribution_insertions.out") || die "[apply-insertions-deletions.pl] $!: $out_dir/insertions/distribution_insertions.out\n";
open(DIST_DEL,">$out_dir/deletions/distribution_deletions.out") || die "[apply-insertions-deletions.pl] $!: $out_dir/deletions/distribution_deletions.out\n";

my %count_orf_status=();
my @dist_disrupted;
my %trans2ins=();
my %trans2del=();
my %all_transcripts=();


my $ins_id = 1;
my $del_id = 1;
my %trans2ins_id=();
my %trans2del_id=();

my %seen_ins=();
my %seen_del=();
my %transcript2fate=();

my ($gvfs, $mod_pfas, $mod_cdnas, $summaries, $mod_peps, $alignments) = (0, 0, 0, 0, 0, 0);
for my $key (keys %chrom_trans)
{
	print "[apply-insertions-deletions.pl] Processing chromosome $key on ";
    print localtime()."\n";
    
#	my @transcripts = split(/\s+/,$chrom_trans{$key});
#	my @insertions = split(/\s+/,$chrom_ins{$key});
#	my @deletions = split(/\s+/,$chrom_del{$key});

	my %seen_dist_ins=();
	my %seen_dist_del=();

	my %ins2gvf_variant=();
	my %del2gvf_variant=();

	for (my $i=0; $chrom_trans{$key} and $i < @{$chrom_trans{$key}};$i++)
	{
		my $transcript = $chrom_trans{$key}->[$i];
		$all_transcripts{$transcript}++;
		print ALIGNMENT "\>$transcript\t$key\t$trans_start{$transcript}\t$trans_end{$transcript}\t$trans_strand{$transcript}\n";
		$alignments ++;
		#print STDERR "Checking transcript $transcript\t";
		#returns insertion coordinates that overlap with this transcript's exon(s)
		my $overlap_ins = get_overlap_ins($key, $trans_exons{$transcript});

		#returns deletion coordinates that overlap with this transcript's exon(s)
		my $overlap_del	= get_overlap_del($key, $trans_exons{$transcript});

		print MOD_PFA "\>$transcript\t$trans_chrom{$transcript}\t$trans_start{$transcript}\t";
		print MOD_PFA "$trans_end{$transcript}\t$trans_strand{$transcript}\n";
		$mod_pfas ++;

        print MOD_cDNA "\>$transcript\t$trans_chrom{$transcript}\t$trans_start{$transcript}\t";
        print MOD_cDNA "$trans_end{$transcript}\t$trans_strand{$transcript}\n";
        $mod_cdnas ++;

		print GENE_SUMMARY "$transcript\t$trans2snp{$transcript}\t";

		#print STDERR "$transcript\t$overlap_ins\t$overlap_del\n";

        my @aux_ref_exons = split(/\n/,$ref_exons{$transcript});
		@aux_ref_exons = reverse(@aux_ref_exons);

      	my @exons = split(/\n/,$trans_exons{$transcript});
        @exons = reverse(@exons);       #here I have the exons going from 3' to 5'

		my %to_delete=();
		my %to_insert=();

		my @dels = split(/\s+/,$overlap_del);
		my @ins = split(/\s+/,$overlap_ins);

		for (my $j=0;$j<@ins;$j++)
		{
#			print "BLA: $key\:$ins[$j]\n";
			my $aux_coord = "$key\:$ins[$j]";
			$seen_ins{$aux_coord}++;
			$trans2ins{$transcript}=1;
			my $end_ins;
			if($trans_strand{$transcript} eq '-')
			{
				$to_ins{"$key\:$ins[$j]"}=~/^(\d+)\-(\S+)$/;
				$end_ins = $1;
				my $aux_ins = $2;
				$aux_ins = reverse($aux_ins);
				$aux_ins=~tr/[ACTGactg]/[TGACtgac]/;
				$to_insert{$ins[$j]} = $aux_ins;
			}
			else
			{
				$to_ins{"$key\:$ins[$j]"}=~/^(\d+)\-(\S+)$/;
				$end_ins = $1;
				$to_insert{$ins[$j]} = $2;
			}
			print GENE_SUMMARY "$key\:$ins[$j]\.\.",$to_ins{"$key\:$ins[$j]"}," ";

			next if (defined $seen_dist_ins{"$key\:$ins[$j]"});
			$seen_dist_ins{"$key\:$ins[$j]"}++;
			my $len_ins = length($to_insert{$ins[$j]});
			$dist_ins_exons[$len_ins]++;
		}
		print GENE_SUMMARY "\t";
		
        for (my $j=0;$j<@dels;$j++)
        {
			$trans2del{$transcript}=1;
			$dels[$j]=~/(\d+)\.\.(\d+)/;
			my $start_del = $1;
			my $end_del = $2;

			$seen_del{"$key:$start_del\.\.$end_del"}++;

			for(my $k=$start_del;$k<=$end_del;$k++)
			{
				$to_delete{$k}++;
			}
			print GENE_SUMMARY "$key\:$dels[$j] ";
        
            next if (defined $seen_dist_del{"$key\:$dels[$j]"}); 
            $seen_dist_del{"$key\:$dels[$j]"}++;
            my $len_del = $end_del - $start_del +1;
            $dist_del_exons[$len_del]++;
        }

		print GENE_SUMMARY "\t";
		
		my $len_before=0;
		my $len_after =0;

		my @ref_ali_exons;
		my @tar_ali_exons;

      	for(my $j=0;$j<@exons;$j++)
      	{
            my @ref_alignment;
        	my @tar_alignment;

            my ($exon_prefix,$exon_start,$exon_end, $exon_seq);
            $exons[$j]=~/^(\((\d+)\-(\d+)\))\s+(\S+)/;
            $exon_prefix=$1;
            $exon_start = $2;
            $exon_end = $3;
            $exon_seq = $4;

			$aux_ref_exons[$j]=~/^\(\d+\-\d+\)\s+(\S+)/;
			my $ref_exon_seq = $1;
			my @ref_nts = split(//,$ref_exon_seq);
			@ref_nts = reverse(@ref_nts);	

			$len_before+=length($exon_seq);               
            
           	my @nucleotides = split(//,$exon_seq);
          	@nucleotides = reverse(@nucleotides);   #since we want to go from 3' to 5'

            my %to_consider=();
   	        my @to_print;

            if($trans_strand{$transcript} eq '+')
   	        {
      	        for (my $k=$exon_end;$k>=$exon_start;$k--)
               	{
				 	$to_consider{$k} = $nucleotides[$exon_end - $k];
                   	if(defined $to_delete{$k})
                   	{
                       	$to_consider{$k}="";
						push @ref_alignment,$ref_nts[$exon_end - $k];
						push @tar_alignment,"\-"; 
                   	}
                   	if(defined $to_insert{$k})
                   	{
                       	$to_consider{$k}.=$to_insert{$k};
						for(my $l=0;$l<length($to_insert{$k});$l++)
						{
							push @ref_alignment,"\-";
						}
                        my @aux_ins = split(//,$to_insert{$k});
                        for(my $l=$#aux_ins;$l>=0;$l--)
                        {
	                        push @tar_alignment,$aux_ins[$l];
                        }
                   	}
					if(!defined $to_delete{$k})
					{
						push @ref_alignment,$ref_nts[$exon_end - $k];
						push @tar_alignment,$nucleotides[$exon_end - $k];
					}
                   	push @to_print,$to_consider{$k};
               	}
			}
	        else
            {
	            for (my $k=$exon_end;$k<=$exon_start;$k++)
	           	{
	               	$to_consider{$k} = $nucleotides[$k-$exon_end];
	                if(!defined $to_delete{$k})
	                {
	                    push @ref_alignment,$ref_nts[$k - $exon_end];
	                    push @tar_alignment,$nucleotides[$k - $exon_end];
	                }
	
	               	if(defined $to_delete{$k})
	               	{
                       	$to_consider{$k}="";  # 2012-10-05 | CF | bugfix; do not consider this base later when reconstructing modified cDNA
	                    push @ref_alignment,$ref_nts[$k - $exon_end];
	                  	push @tar_alignment,"\-";
	               	}
		        	if(defined $to_insert{$k})
	               	{
	    				$to_consider{$k}=$to_insert{$k} . $to_consider{$k};
	                    for(my $l=0;$l<length($to_insert{$k});$l++)
	                    {
	                    	push @ref_alignment,"\-";
						}
						my @aux_ins = split(//,$to_insert{$k});
						for(my $l=$#aux_ins;$l>=0;$l--)
						{
							push @tar_alignment,$aux_ins[$l];
						}
					}
					push @to_print,$to_consider{$k};
				}
           	}
			my $new_exon="";
			for(my $k=$#to_print;$k>=0;$k--)
			{
				$new_exon.=$to_print[$k];
			}
			$exons[$j]=~s/$exon_seq/$new_exon/;

			my $ref_ali_exon=$exon_prefix . "\t";
			my $len_prefix = length($exon_prefix);

			my $tar_ali_exon="";
			for(my $k=0;$k<$len_prefix;$k++)
			{
				$tar_ali_exon.=" ";
			}
			$tar_ali_exon.="\t";
	
			my $count_spacing = 0;
			for (my $k=$#ref_alignment;$k>=0;$k--)
			{
				if($count_spacing == 80)
				{
					$ref_ali_exon.="\n";
					for(my $l=0;$l<$len_prefix;$l++)
					{
						$ref_ali_exon.=" ";
					}
					$ref_ali_exon.="\t";
					$count_spacing=0;
				}	
				$ref_ali_exon.=$ref_alignment[$k];
				$count_spacing++;
			}


			$count_spacing = 0;
  	        for (my $k=$#tar_alignment;$k>=0;$k--)
          	{ 
				if($count_spacing == 80)
				{
					$tar_ali_exon.="\n";
					for(my $l=0;$l<$len_prefix;$l++)
					{ 
						$tar_ali_exon.=" ";
					}
					$tar_ali_exon.="\t";
					$count_spacing=0;
				}
				$tar_ali_exon.=$tar_alignment[$k];
				$count_spacing++;
			}

			unshift @ref_ali_exons,$ref_ali_exon;
			unshift @tar_ali_exons,$tar_ali_exon;
		}
		
		for (my $j=0;$j<@ref_ali_exons;$j++)
		{
			my @ref_segments = split(/\n/,$ref_ali_exons[$j]);
			my @tar_segments = split(/\n/,$tar_ali_exons[$j]);
			for(my $k=0;$k<@ref_segments;$k++)
			{
				print ALIGNMENT "$ref_segments[$k]\n$tar_segments[$k]\n\n";
				$alignments ++;
			}
		}

		my $modified_cDNA="";
      	@exons = reverse(@exons);
       	for(my $j=0;$j<@exons;$j++)
       	{
           	print MOD_PFA $exons[$j],"\n";
           	$mod_pfas ++;
            $exons[$j]=~/^\S+\s+(\S*)/;
			my $aux_seq = $1;
            print MOD_cDNA $aux_seq;
			$modified_cDNA.=$aux_seq;
      	}
        print MOD_cDNA "\n"; # 2012-10-07 | CF | add newline character after sequence
       	$mod_cdnas ++;
        	
		$len_after = length($modified_cDNA);

		my @assessed_pep = evaluate_cDNA($modified_cDNA,$transcript);

		my $mod_pep_seq= $assessed_pep[0];
		my $mod_pep_status = $assessed_pep[1];
		my $mod_pep_perc = $assessed_pep[2];

		if(!defined $trans2snp{$transcript} && 
		!defined $trans2ins{$transcript} && 
		!defined $trans2del{$transcript})
		{
			$mod_pep_status = "ORF_INTACT";
		}

		# 2012-10-11 | CF | determine indels impacting splice junctions
		my (%sj_dels, %sj_ins);
		if (exists $sj{$transcript})
		{
			my @tsjs = split(" ", $sj{$transcript});
			my ($junctions_acc, $junctions_don) = ("", "");
			foreach my $tsj (@tsjs)
			{
				my @f = split(";", $tsj);
				my ($chr, $start, $end) = $f[0] =~ /([^:]+):(\d+)..(\d+)/;
				if ($f[2] eq 'donor')
				{
					$junctions_don .= "($start-$end) \n";
				}
				elsif ($f[2] eq 'acceptor')
				{
					$junctions_acc .= "($start-$end) \n";
				}
				else
				{
					die "[apply-insertions-deletions.pl] ERROR: unknown splice site type: $tsj\n";
				}
			}
#			print "DONORS: $junctions_don\n";
#			print "ACCEPTORS: $junctions_acc\n";
			
			# insertions impacting splice acceptor sites
			{
				my $ods = get_overlap_ins($key, $junctions_acc);
				foreach my $od (split(/\s+/,$ods))
				{
#					print "INS HIT ACCEPTOR: $key:$od\n";
					$ins2gvf_variant{"$od\_\_splice\_acceptor\_variant"} .= $transcript . "\,";
					$seen_ins{"$key:$od"} = 1;
					$sj_ins{"$key:$od"} = 1;
					$trans2splice{$transcript} = 1;
				}				
			}
			# insertions impacting splice donor sites
			{
				my $ods = get_overlap_ins($key, $junctions_don);
				foreach my $od (split(/\s+/,$ods))
				{
#					print "INS HIT DONOR: $key:$od\n";
					$ins2gvf_variant{"$od\_\_splice\_donor\_variant"} .= $transcript . "\,";
					$seen_ins{"$key:$od"} = 1;
					$sj_ins{"$key:$od"} = 1;
					$trans2splice{$transcript} = 1;
				}				
			}
			# deletions impacting splice acceptor sites
			{
				my $ods = get_overlap_del($key, $junctions_acc);
				foreach my $od (split(/\s+/,$ods))
				{
#					print "DEL HIT ACCEPTOR: $key:$od\n";
					$del2gvf_variant{"$od\_\_splice\_acceptor\_variant"} .= $transcript . "\,";
					$seen_del{"$key:$od"} = 1;
					$sj_dels{"$key:$od"} = 1;
					$trans2splice{$transcript} = 1;
				}				
			}
			# deletions impacting splice donor sites
			{
				my $ods = get_overlap_del($key, $junctions_don);
				foreach my $od (split(/\s+/,$ods))
				{
#					print "DEL HIT DONOR: $key:$od\n";
					$del2gvf_variant{"$od\_\_splice\_donor\_variant"} .= $transcript . "\,";
					$seen_del{"$key:$od"} = 1;
					$sj_dels{"$key:$od"} = 1;
					$trans2splice{$transcript} = 1;
				}
				
			}
		}
		
		for(my $j=0;$j<@ins;$j++)
		{
        	# 2012-10-11 | CF | skip insertions already classified as splice junction variants
        	next if (exists $sj_ins{"$key:$ins[$j]"});

			# 2012-10-05 | CF | classify insertion into frameshift or inframe based on length, not based on ORF status; the ORF status remains unchanged though
			$to_ins{"$key\:$ins[$j]"}=~/^(\d+)\-(\S+)$/;
			my $len = length($2);
			if ($len % 3 == 0)
			{
				$ins2gvf_variant{"$ins[$j]\_\_inframe\_variant"}.=$transcript . "\,";
			}
			else
			{
				$ins2gvf_variant{"$ins[$j]\_\_frameshift\_variant"}.=$transcript . "\,";				
			}
		}

        for(my $j=0;$j<@dels;$j++)
        {
        	# 2012-10-11 | CF | skip deletions already classified as splice junction variants
        	next if (exists $sj_dels{"$key:$dels[$j]"});
        	
			# 2012-10-05 | CF | classify deletion into frameshift or inframe based on length, not based on ORF status; the ORF status remains unchanged though
        	my ($start_del, $end_del) = $dels[$j]=~/(\d+)\.\.(\d+)/;
        	my $len = abs($end_del-$start_del)+1;
			if ($len % 3 == 0)
			{
            	$del2gvf_variant{"$dels[$j]\_\_inframe\_variant"}.=$transcript . "\,";
			}
			else
			{
            	$del2gvf_variant{"$dels[$j]\_\_frameshift\_variant"}.=$transcript . "\,";
			}
        }

		# 2012-06-19 | CF | ORFs impacted by splice site variants are considered disrupted
		$mod_pep_status = 'ORF_DISRUPTED'
			if ($mod_pep_status eq 'ORF_PRESERVED' and $trans2splice{$transcript});

		$transcript2fate{$transcript} =  $mod_pep_status;

		print GENE_SUMMARY "$len_before\t$len_after\t$mod_pep_status\t";

		#if there is a case of ne NA and eq ORF_INTACT means that the reference protein comes with an internal stop codon
		if($mod_pep_perc ne 'N.A.' && $mod_pep_status ne 'ORF_INTACT')
		{
			print GENE_SUMMARY $mod_pep_perc;
			$dist_disrupted[floor($mod_pep_perc/10)]++;
			$transcript2fate{$transcript}.="\($mod_pep_perc\)";
		}
		$count_orf_status{$mod_pep_status}++;

		print GENE_SUMMARY "\n";
		$summaries ++;
        print MOD_PEP "\>$transcript\n";
        print MOD_PEP $mod_pep_seq,"\n";
        $mod_peps ++;
		#print peptide sequence in modified fasta file for proteins
		#PRINT status in gene summary
	}

	for my $key2 (keys %ins2gvf_variant)
	{
		my ($start_ins,$end_ins,$seq_ins,$variant);
		if($key2=~/(\d+)\_\_(\S+)/)
		{
			$start_ins = $1;
			$variant = $2;
		}
		if($to_ins{"$key\:$start_ins"}=~/^(\d+)\-(\S+)$/)
		{
			$end_ins = $1;
			$seq_ins = $2;
		}
		print GVF "$key\tCooVar\tinsertion\t$start_ins\t$end_ins\t\.\t\+\t\.\t";
		print GVF "ID\=ins\_$ins_id\;Variant\_seq\=$seq_ins\;Reference\_seq\=\~\;Variant\_type\=$variant;Variant\_effect\=$variant"; # 2011-11-21 | CF | added variant_type for easier track grouping in Gbrowse
		chop($ins2gvf_variant{$key2});
		print GVF " 0 mRNA $ins2gvf_variant{$key2}\n";
		$gvfs ++;

		my @trans = split(/\,/,$ins2gvf_variant{$key2});
		for(my $j=0;$j<@trans;$j++)
		{
			$trans2ins_id{$trans[$j]}.= "ins\_$ins_id\,"; 
		}

		$ins_id++;
	}

	for(my $j=0;$chrom_ins{$key} and $j<@{$chrom_ins{$key}};$j++)
	{
		my $seen_key = "$key:".$chrom_ins{$key}->[$j];
		$seen_key =~ s/\.\..*//;
		next if (defined $seen_ins{$seen_key});
#		print "NOT SEEN: $seen_key\n";
		#then it is a silent mutation
		$chrom_ins{$key}->[$j]=~/(\d+)\.\.(\d+)\-(\S+)/;

		my $start_ins = $1;
		my $end_ins = $2;
		my $seq_ins = $3;

		print GVF "$key\tCooVar\tinsertion\t$start_ins\t$end_ins\t\.\t\+\t\.\t";
		print GVF "ID\=ins\_$ins_id\;Variant\_seq\=$seq_ins\;Reference\_seq\=\~\;Variant\_type\=silent\_mutation;Variant\_effect\=silent\_mutation\n"; # 2011-11-21 | CF | added variant_type for easier track grouping in Gbrowse
		$ins_id++;
		$gvfs ++;

	}

	for my $key2 (keys %del2gvf_variant)
	{
		$key2=~/(\d+)\.\.(\d+)\_\_(\S+)/;
		my $start_del =	$1;
		my $end_del = $2;
		my $variant = $3;

		print GVF "$key\tCooVar\tdeletion\t$start_del\t$end_del\t\.\t\+\t\.\t";
		print GVF "ID\=del\_$del_id\;Variant\_seq\=\-\;Reference\_seq\=\~\;Variant\_type\=$variant;Variant\_effect\=$variant";  # 2011-11-21 | CF | added variant_type for easier track grouping in Gbrowse
		chop($del2gvf_variant{$key2});
		print GVF " 0 mRNA $del2gvf_variant{$key2}\n";
		$gvfs ++;
		my @trans = split(/\,/,$del2gvf_variant{$key2});
		for(my $j=0;$j<@trans;$j++)
		{
			$trans2del_id{$trans[$j]}.= "del\_$del_id\,";
		}
		$del_id++;
	}

	for(my $j=0;$chrom_del{$key} and $j<@{$chrom_del{$key}};$j++)
	{
		next if (defined $seen_del{"$key:".$chrom_del{$key}->[$j]});
		#then it is a silent mutation
		$chrom_del{$key}->[$j]=~/(\d+)\.\.(\d+)/;

		my $start_del = $1;
		my $end_del = $2;

		print GVF "$key\tCooVar\tdeletion\t$start_del\t$end_del\t\.\t\+\t\.\t";
		print GVF "ID\=del\_$del_id\;Variant\_seq\=\-\;Reference\_seq\=\~\;Variant\_type\=silent\_mutation;Variant\_effect\=silent\_mutation\n"; # 2011-11-21 | CF | added variant_type for easier track grouping in Gbrowse
		$del_id++;
		$gvfs ++;
	}
}
close(MOD_cDNA);
close(MOD_PFA);
close(MOD_PEP);
close(GENE_SUMMARY);
close(GVF);
close(ALIGNMENT);

print "[apply-insertions-deletions.pl] $gvfs lines added to $gvf_gvs\n";
print "[apply-insertions-deletions.pl] $mod_pfas lines written to $out_dir/transcripts/variant_cdna.exons\n";
print "[apply-insertions-deletions.pl] $mod_cdnas lines written to $out_dir/transcripts/variant_cdna.fasta\n";
print "[apply-insertions-deletions.pl] $mod_peps lines written to $out_dir/transcripts/variant_peptides.fasta\n";
print "[apply-insertions-deletions.pl] $summaries lines written to $out_dir/intermediate-files/variant.summary\n";
print "[apply-insertions-deletions.pl] $alignments lines written to $out_dir/transcripts/reference_variant.alignment\n";

open(FINAL_GFF3,">$out_dir/transcripts.gff3") || die "[apply-insertions-deletions.pl] $!: $out_dir/transcripts.gff3\n";

while(<TRANS_GFF3>)
{
        chomp($_);
        if($_=~/^\#/ || $_=~/\tCooVar\tCDS\t/)
        {
        	print FINAL_GFF3 $_,"\n";
            next;
        }
        else
        {
			my @line = split(/\t/,$_);
            $line[8]=~/ID\=([^;\s]+)/;
            print "[apply-insertions-deletions.pl] ERROR: could not parse ID from following line:\n$_\n" if (!defined $1);
            my $trans = $1;
            print FINAL_GFF3 $_,"$trans2ins_id{$trans}$trans2del_id{$trans}$transcript2fate{$trans}\n";
        }
}
close(TRANS_GFF3);
close(FINAL_GFF3);

for my $key (keys %count_orf_status)
{
        print STAT "$key\t$count_orf_status{$key}\n";
}

print STAT "\n\n";

print STAT "Distribution of Early Stop Codons on Peptides\n";
print STAT "--------------------------------------------\n";

for (my $i=0;$i<@dist_disrupted;$i++)
{
        next if (!($dist_disrupted[$i]>=1));
        print STAT "\[" , $i*10 , "\," , ($i+1)*10 , "\)\t$dist_disrupted[$i]\n";
}
print STAT "\n";

my %count_overall_impact=();

for my $key (keys %all_transcripts)
{
        my $overall="";
        if(defined $trans2snp{$key})
        {
                $overall.="SNP";
        }
        if(defined $trans2ins{$key})
        {
                $overall.=" INS";
        }
        if(defined $trans2del{$key})
        {
                $overall.=" DEL";
        }
	if($overall eq "")
	{
		$overall = "NONE";
	}
        $count_overall_impact{$overall}++;
}

print STAT "OVERALL SUMMARY IMPACT OF GVs ON TRANSCRIPTS\n";
print STAT "--------------------------------------------\n";
for my $key (sort keys %count_overall_impact)
{
        print STAT "$key\t$count_overall_impact{$key}\n";
}

close(STAT);
print "[apply-insertions-deletions.pl] Variant statistics written to $out_dir/transcripts/variant.stat\n";

print DIST_INS "Length Distribution of All Insertions\n";
print DIST_INS "-------------------------------------\n";

for (my $i=0;$i<@dist_ins;$i++)
{
	next if (!($dist_ins[$i]>=1));
	print DIST_INS "$i\t$dist_ins[$i]\n";
}
print DIST_INS "\n";

print DIST_INS "Length Distribution of Insertions impacting transcripts\n";
print DIST_INS "-------------------------------------------------------\n";


for (my $i=0;$i<@dist_ins_exons;$i++)
{ 
        next if (!($dist_ins_exons[$i]>=1));
        print DIST_INS "$i\t$dist_ins_exons[$i]\n";
}


print DIST_DEL "Length Distribution of All Deletions\n";
print DIST_DEL "------------------------------------\n";

for (my $i=0;$i<@dist_del;$i++)
{ 
        next if (!($dist_del[$i]>=1));
        print DIST_DEL "$i\t$dist_del[$i]\n";
}
print DIST_DEL "\n";

print DIST_DEL "Length Distribution of Deletions impacting transcripts\n";
print DIST_DEL "-------------------------------------------------------\n";

for (my $i=0;$i<@dist_del_exons;$i++)
{ 
        next if (!($dist_del_exons[$i]>=1));
        print DIST_DEL "$i\t$dist_del_exons[$i]\n";
}

close(DIST_INS);
close(DIST_DEL);

print "[apply-insertions-deletions.pl] Insertion distribution written to $out_dir/insertions/distribution_insertions.out\n";
print "[apply-insertions-deletions.pl] Deletion distribution written to $out_dir/deletions/distribution_deletions.out\n";


print "[apply-insertions-deletions.pl] Done at ";
print localtime()."\n";
