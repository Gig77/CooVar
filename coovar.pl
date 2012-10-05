#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
use Bio::DB::Fasta;

my $exon_file;
my $dna_file;
my $snp_file;
my $ins_file;
my $del_file;
my $tab_file;
my $vcf_file;
my $out_dir;
my $circos='';

sub getTimestamp
{
	my ($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time);
	$Year += 1900;
	$Month++;
	$Month = "0$Month" if ($Month < 10);	
	$Day = "0$Day" if ($Day < 10);	
	$Hour = "0$Hour" if ($Hour < 10);	
	$Minute = "0$Minute" if ($Minute < 10);	
	$Second = "0$Second" if ($Second < 10);
		
	return "$Year$Month$Day"."_$Hour$Minute$Second";
}

sub execute
{
	my $cmd = shift;
	system($cmd) == 0 or die "[coovar.pl] ERROR executing command: $cmd\n.";
}

my $params = join(" ", @ARGV);
GetOptions ("exon|e=s" => \$exon_file,
		"ref|r=s"  => \$dna_file,
		"tab|t=s"  => \$tab_file,
		"vcf|v=s"  => \$vcf_file,
		"out|o=s"  => \$out_dir,
		"circos"   => \$circos
		);	

unless ($exon_file && $dna_file && ($vcf_file || $tab_file))
{
    print "USAGE: ./coovar.pl -e EXONS_GFF -r REFERENCE_FASTA (-t GVS_TAB_FORMAT | -v GVS_VCF_FORMAT) [-o OUTPUT_DIRECTORY] [--circos]\n";
    print "Program parameter details provided in file README.\n";
    exit;
}

$out_dir = cwd()."/output_".getTimestamp() if (!$out_dir);

my ($prog_dir) = __FILE__ =~ /(.*)\//;
$prog_dir = "." if (!$prog_dir);

print "[coovar.pl] Start executing script on ";
system("date");

print "[coovar.pl] Program directory: $prog_dir\n";
print "[coovar.pl] Command line: ".__FILE__." $params\n";
print "[coovar.pl]   REFERENCE: $dna_file\n";
print "[coovar.pl]   GVS_TAB_FORMAT: $tab_file\n";
print "[coovar.pl]   GVS_VCF_FORMAT: $vcf_file\n";
print "[coovar.pl]   OUTPUT DIRECTORY: $out_dir\n";
print "[coovar.pl]   CIRCOS FLAG: $circos\n";

die "[coovar.pl] ERROR: EXONS_GFF file $exon_file not found.\n" unless (-e $exon_file);
die "[coovar.pl] ERROR: REFERENCE file $dna_file not found.\n" unless (-e $dna_file);
die "[coovar.pl] ERROR: GVS_TAB_FORMAT file $tab_file not found.\n" unless (!$tab_file or -e $tab_file);
die "[coovar.pl] ERROR: GVS_VCF_FORMAT file $vcf_file not found.\n" unless (!$vcf_file or -e $vcf_file);

execute("mkdir $out_dir") unless (-d $out_dir);
execute("mkdir $out_dir/CV_Intermediate_Files") unless (-d "$out_dir/CV_Intermediate_Files");

# index reference FASTA
print "[coovar.pl] Indexing FASTA file $dna_file on ";
system("date");
my $db = Bio::DB::Fasta->new($dna_file);
die ("[coovar.pl]   ERROR: Could not index FASTA file. Do you have write permissions to the directory containing the FASTA file?\n")
	if (!$db or !-e "$dna_file.index");
print "[coovar.pl]   Done indexing FASTA file on ";
system("date");

print "[coovar.pl] Extracting contig information from FASTA on ";
system("date");
open(CHR, ">$out_dir/CV_Intermediate_Files/contigs.summary")
	or die("[coovar.pl]   ERROR: Could not write contig information to file $out_dir/CV_Intermediate_Files/contigs.summary\n");
my $stream = $db->get_PrimarySeq_stream();
while(my $seq = $stream->next_seq())
{
	print CHR $seq->id."\t".$seq->length."\n";
}
close(CHR);

print "[coovar.pl] Parsing GV files on ";
system("date");

if(defined $tab_file)
{
	execute("perl $prog_dir/scripts/TAB2CV.pl $tab_file $out_dir/CV_Intermediate_Files");
	$snp_file = "$out_dir/CV_Intermediate_Files/".basename($tab_file) . '.CV_SNP';
	$ins_file = "$out_dir/CV_Intermediate_Files/".basename($tab_file) . '.CV_INS';
	$del_file = "$out_dir/CV_Intermediate_Files/".basename($tab_file) . '.CV_DEL';
}
if(defined $vcf_file)
{
	execute("perl $prog_dir/scripts/VCF2CV.pl $vcf_file $out_dir/CV_Intermediate_Files");
    $snp_file = "$out_dir/CV_Intermediate_Files/".basename($vcf_file) .	'.CV_SNP';
    $ins_file = "$out_dir/CV_Intermediate_Files/".basename($vcf_file) .	'.CV_INS';
    $del_file = "$out_dir/CV_Intermediate_Files/".basename($vcf_file) .	'.CV_DEL';
}

execute("mkdir $out_dir/CV_SNPs") unless (-d "$out_dir/CV_SNPs");
execute("mkdir $out_dir/CV_Insertions") unless (-d "$out_dir/CV_Insertions");
execute("mkdir $out_dir/CV_Deletions") unless (-d "$out_dir/CV_Deletions");
execute("mkdir $out_dir/CV_Transcripts") unless (-d "$out_dir/CV_Transcripts");

print "[coovar.pl] Checking consistency of GV files on ";
system("date");
execute("perl $prog_dir/scripts/revise-GV-files.pl $snp_file $ins_file $del_file $out_dir");

my $ex_snps = "$out_dir/CV_SNPs/excluded_" . basename($snp_file);
my $ex_ins = "$out_dir/CV_Insertions/excluded_" . basename($ins_file);
my $ex_del = "$out_dir/CV_Deletions/excluded_" . basename($del_file);

$snp_file = "$out_dir/CV_SNPs/kept_" . basename($snp_file);
$ins_file = "$out_dir/CV_Insertions/kept_" . basename($ins_file);
$del_file = "$out_dir/CV_Deletions/kept_" . basename($del_file);

print "[coovar.pl] Extracting cDNA sequence from reference on ";
system("date");
execute("perl $prog_dir/scripts/extract-cDNA.pl $exon_file $dna_file $out_dir");

print "[coovar.pl] Applying SNPs to cDNA sequences on ";
system("date");
execute("perl $prog_dir/scripts/apply-SNPs.pl $out_dir/CV_Transcripts/reference_cDNA.exons $snp_file $out_dir");

print "[coovar.pl] Categorizing SNPs on ";
system("date");
execute("perl $prog_dir/scripts/categorize-SNPs.pl $out_dir/CV_Transcripts/reference_cDNA.exons $out_dir/CV_Intermediate_Files/target_cDNA_SNPs.exons $prog_dir/lib/grantham_matrix  $snp_file $out_dir/CV_Intermediate_Files/splice_junctions.tmp $out_dir/CV_Intermediate_Files/transcripts.gff3.tmp $out_dir/CV_Intermediate_Files/snps_frequency.table $out_dir");

print "[coovar.pl] Generating stats for SNP categorization on ";
system("date");
execute("perl $prog_dir/scripts/get-stats-SNPs.pl $out_dir/CV_Intermediate_Files/categorized_snp_coords.list $out_dir");

print "[coovar.pl] Analyzing insertions and deletions on ";
system("date");
execute("perl $prog_dir/scripts/apply-insertions-deletions.pl $out_dir/CV_Transcripts/reference_cDNA.exons $out_dir/CV_Intermediate_Files/target_cDNA_SNPs.exons $out_dir/CV_Intermediate_Files/categorized_snp_coords.list $ins_file $del_file $out_dir/CV_categorized_GVs.gvf $out_dir/CV_Intermediate_Files/transcripts_snps_applied.gff3.tmp $out_dir");

if($circos == 1)
{
	print "[coovar.pl] Generating CIRCOS files for SNPs, insertions, deletions and coding regions on ";
	system("date");
	execute("perl $prog_dir/scripts/parse2circos.pl $snp_file $ins_file $del_file $exon_file $dna_file $out_dir"); 
}

print "[coovar.pl] Categorized GVs can be found in $out_dir/CV_categorized_GVs.gvf\n";
print "[coovar.pl] Categorized transcripts can be found in $out_dir/CV_transcripts.gff3\n";

print "[coovar.pl] Done at ";
system("date");
