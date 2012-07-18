#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;

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
	system($cmd) == 0 or die "[variant-analyzer.pl] ERROR executing command: $cmd\n.";
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
    print "USAGE: ./variant-analyzer.pl -e EXONS_GFF -r REFERENCE_FASTA (-t GVS_TAB_FORMAT | -v GVS_VCF_FORMAT) [-o OUTPUT_DIRECTORY] [--circos]\n";
    print "Program parameter details provided in file VA_README.\n";
    exit;
}

$out_dir = cwd()."/output_".getTimestamp() if (!$out_dir);

my ($prog_dir) = __FILE__ =~ /(.*)\//;
$prog_dir = "." if (!$prog_dir);

print "[variant-analyzer.pl] Start executing script on ";
system("date");

print "[variant-analyzer.pl] Program directory: $prog_dir\n";
print "[variant-analyzer.pl] Command line: ".__FILE__." $params\n";
print "[variant-analyzer.pl]   REFERENCE: $dna_file\n";
print "[variant-analyzer.pl]   GVS_TAB_FORMAT: $tab_file\n";
print "[variant-analyzer.pl]   GVS_VCF_FORMAT: $vcf_file\n";
print "[variant-analyzer.pl]   OUTPUT DIRECTORY: $out_dir\n";
print "[variant-analyzer.pl]   CIRCOS FLAG: $circos\n";

die "[variant-analyzer.pl] ERROR: EXONS_GFF file $exon_file not found.\n" unless (-e $exon_file);
die "[variant-analyzer.pl] ERROR: REFERENCE file $dna_file not found.\n" unless (-e $dna_file);
die "[variant-analyzer.pl] ERROR: GVS_TAB_FORMAT file $tab_file not found.\n" unless (!$tab_file or -e $tab_file);
die "[variant-analyzer.pl] ERROR: GVS_VCF_FORMAT file $vcf_file not found.\n" unless (!$vcf_file or -e $vcf_file);

execute("mkdir $out_dir") unless (-d $out_dir);
execute("mkdir $out_dir/VA_Intermediate_Files") unless (-d "$out_dir/VA_Intermediate_Files");

print "[variant-analyzer.pl] Parsing GV files on ";
system("date");

if(defined $tab_file)
{
	execute("perl $prog_dir/scripts/TAB2VA.pl $tab_file $out_dir/VA_Intermediate_Files");
	$snp_file = "$out_dir/VA_Intermediate_Files/".basename($tab_file) . '.VA_SNP';
	$ins_file = "$out_dir/VA_Intermediate_Files/".basename($tab_file) . '.VA_INS';
	$del_file = "$out_dir/VA_Intermediate_Files/".basename($tab_file) . '.VA_DEL';
}
if(defined $vcf_file)
{
	execute("perl $prog_dir/scripts/VCF2VA.pl $vcf_file $out_dir/VA_Intermediate_Files");
    $snp_file = "$out_dir/VA_Intermediate_Files/".basename($vcf_file) .	'.VA_SNP';
    $ins_file = "$out_dir/VA_Intermediate_Files/".basename($vcf_file) .	'.VA_INS';
    $del_file = "$out_dir/VA_Intermediate_Files/".basename($vcf_file) .	'.VA_DEL';
}

execute("mkdir $out_dir/VA_SNPs") unless (-d "$out_dir/VA_SNPs");
execute("mkdir $out_dir/VA_Insertions") unless (-d "$out_dir/VA_Insertions");
execute("mkdir $out_dir/VA_Deletions") unless (-d "$out_dir/VA_Deletions");
execute("mkdir $out_dir/VA_Transcripts") unless (-d "$out_dir/VA_Transcripts");

print "[variant-analyzer.pl] Checking consistency of GV files on ";
system("date");
execute("perl $prog_dir/scripts/revise-GV-files.pl $snp_file $ins_file $del_file $out_dir");

my $ex_snps = "$out_dir/VA_SNPs/excluded_" . basename($snp_file);
my $ex_ins = "$out_dir/VA_Insertions/excluded_" . basename($ins_file);
my $ex_del = "$out_dir/VA_Deletions/excluded_" . basename($del_file);

$snp_file = "$out_dir/VA_SNPs/kept_" . basename($snp_file);
$ins_file = "$out_dir/VA_Insertions/kept_" . basename($ins_file);
$del_file = "$out_dir/VA_Deletions/kept_" . basename($del_file);

print "[variant-analyzer.pl] Extracting cDNA sequence from reference on ";
system("date");
execute("perl $prog_dir/scripts/extract-cDNA.pl $exon_file $dna_file $out_dir");

print "[variant-analyzer.pl] Applying SNPs to cDNA sequences on ";
system("date");
execute("perl $prog_dir/scripts/apply-SNPs.pl $out_dir/VA_Transcripts/reference_cDNA.exons $snp_file $out_dir");

print "[variant-analyzer.pl] Categorizing SNPs on ";
system("date");
execute("perl $prog_dir/scripts/categorize-SNPs.pl $out_dir/VA_Transcripts/reference_cDNA.exons $out_dir/VA_Intermediate_Files/target_cDNA_SNPs.exons $prog_dir/lib/grantham_matrix  $snp_file $out_dir/VA_Intermediate_Files/splice_junctions.tmp $out_dir/VA_Intermediate_Files/transcripts.gff3.tmp $out_dir/VA_Intermediate_Files/snps_frequency.table $out_dir");

print "[variant-analyzer.pl] Generating stats for SNP categorization on ";
system("date");
execute("perl $prog_dir/scripts/get-stats-SNPs.pl $out_dir/VA_Intermediate_Files/categorized_snp_coords.list $out_dir");

print "[variant-analyzer.pl] Analyzing insertions and deletions on ";
system("date");
execute("perl $prog_dir/scripts/apply-insertions-deletions.pl $out_dir/VA_Transcripts/reference_cDNA.exons $out_dir/VA_Intermediate_Files/target_cDNA_SNPs.exons $out_dir/VA_Intermediate_Files/categorized_snp_coords.list $ins_file $del_file $out_dir/VA_categorized_GVs.gvf $out_dir/VA_Intermediate_Files/transcripts_snps_applied.gff3.tmp $out_dir");

if($circos == 1)
{
	print "[variant-analyzer.pl] Generating CIRCOS files for SNPs, insertions, deletions and coding regions on ";
	system("date");
	execute("perl $prog_dir/scripts/parse2circos.pl $snp_file $ins_file $del_file $exon_file $dna_file $out_dir"); 
}

print "[variant-analyzer.pl] Done at ";
system("date");
