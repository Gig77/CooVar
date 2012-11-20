#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
use Bio::DB::Fasta;
use Config;

our $VERSION = '0.07';

my $exon_file;
my $dna_file;
my $snp_file;
my $ins_file;
my $del_file;
my $tab_file;
my $vcf_file;
my $out_dir;
my $circos='';
my $feature_source = '';
my $feature_type = 'CDS';
my $no_contig_sum = 0;

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
		"circos"   => \$circos,
		"feature_source=s" => \$feature_source,
		"feature_type=s"   => \$feature_type,
		"no_contig_sum" => \$no_contig_sum  # undocumented parameter for test purposes
		);	

unless ($exon_file && $dna_file && ($vcf_file || $tab_file))
{
    print "USAGE: ./coovar.pl -e EXONS_GFF -r REFERENCE_FASTA (-t GVS_TAB_FORMAT | -v GVS_VCF_FORMAT) [-o OUTPUT_DIRECTORY] [--circos] [--feature_source] [--feature_type]\n";
    print "Program parameter details provided in file README.\n";
    exit;
}

$out_dir = cwd()."/output_".getTimestamp() if (!$out_dir);
$out_dir =~ s/\/$//; # remove trailing '/'

my ($prog_dir) = __FILE__ =~ /(.*)\//;
$prog_dir = "." if (!$prog_dir);

print "[coovar.pl] Start executing script on ";
print localtime()."\n";

print "[coovar.pl] Operating system: $Config{'osname'} $Config{'archname'}\n";
print "[coovar.pl] Program directory: $prog_dir\n";
print "[coovar.pl] Program version: $VERSION\n";
print "[coovar.pl] Command line: ".__FILE__." $params\n";
print "[coovar.pl]   REFERENCE: $dna_file\n";
print "[coovar.pl]   CODING EXONS: $exon_file\n";
print "[coovar.pl]   GVS_TAB_FORMAT: $tab_file\n";
print "[coovar.pl]   GVS_VCF_FORMAT: $vcf_file\n";
print "[coovar.pl]   OUTPUT DIRECTORY: $out_dir\n";
print "[coovar.pl]   CIRCOS FLAG: $circos\n";
print "[coovar.pl]   FEATURE SOURCE: $feature_source\n";
print "[coovar.pl]   FEATURE TYPE: $feature_type\n";

die "[coovar.pl] ERROR: EXONS_GFF file $exon_file not found.\n" unless (-e $exon_file);
die "[coovar.pl] ERROR: REFERENCE file $dna_file not found.\n" unless (-e $dna_file);
die "[coovar.pl] ERROR: GVS_TAB_FORMAT file $tab_file not found.\n" unless (!$tab_file or -e $tab_file);
die "[coovar.pl] ERROR: GVS_VCF_FORMAT file $vcf_file not found.\n" unless (!$vcf_file or -e $vcf_file);

mkdir("$out_dir") unless (-d $out_dir);
mkdir("$out_dir/intermediate-files") unless (-d "$out_dir/intermediate-files");

# index reference FASTA
print "[coovar.pl] Indexing FASTA file $dna_file on ";
print localtime()."\n";
my $db = Bio::DB::Fasta->new($dna_file);
die ("[coovar.pl]   ERROR: Could not index FASTA file. Do you have write permissions to the directory containing the FASTA file?\n")
	if (!$db or !-e "$dna_file.index");
print "[coovar.pl]   Done indexing FASTA file on ";
print localtime()."\n";

if (!$no_contig_sum or !-e "$out_dir/intermediate-files/contigs.summary")
{
	print "[coovar.pl] Extracting contig information from FASTA on ";
	print localtime()."\n";
	open(CHR, ">$out_dir/intermediate-files/contigs.summary")
		or die("[coovar.pl]   ERROR: Could not write contig information to file $out_dir/intermediate-files/contigs.summary\n");
	my $stream = $db->get_PrimarySeq_stream();
	while(my $seq = $stream->next_seq())
	{
		print CHR $seq->id."\t".$seq->length."\n";
	}
	close(CHR);	
}
$db = undef; # 2012-10-22 | CF | free file handles; otherwise sequence fetch would fail after indexing in extract-cdna.pl


print "[coovar.pl] Parsing GV files on ";
print localtime()."\n";

if(defined $tab_file)
{
	execute("perl $prog_dir/scripts/tab2cv.pl $tab_file $out_dir/intermediate-files");
	$snp_file = "$out_dir/intermediate-files/".basename($tab_file) . '.snp';
	$ins_file = "$out_dir/intermediate-files/".basename($tab_file) . '.ins';
	$del_file = "$out_dir/intermediate-files/".basename($tab_file) . '.del';
}
if(defined $vcf_file)
{
	execute("perl $prog_dir/scripts/vcf2cv.pl $vcf_file $out_dir/intermediate-files");
    $snp_file = "$out_dir/intermediate-files/".basename($vcf_file) .	'.snp';
    $ins_file = "$out_dir/intermediate-files/".basename($vcf_file) .	'.ins';
    $del_file = "$out_dir/intermediate-files/".basename($vcf_file) .	'.del';
}

mkdir("$out_dir/snps") unless (-d "$out_dir/snps");
mkdir("$out_dir/insertions") unless (-d "$out_dir/insertions");
mkdir("$out_dir/deletions") unless (-d "$out_dir/deletions");
mkdir("$out_dir/transcripts") unless (-d "$out_dir/transcripts");

print "[coovar.pl] Checking consistency of GV files on ";
print localtime()."\n";
execute("perl $prog_dir/scripts/revise-gv-files.pl $snp_file $ins_file $del_file $out_dir");

my $ex_snps = "$out_dir/snps/excluded_" . basename($snp_file);
my $ex_ins = "$out_dir/insertions/excluded_" . basename($ins_file);
my $ex_del = "$out_dir/deletions/excluded_" . basename($del_file);

$snp_file = "$out_dir/snps/kept_" . basename($snp_file);
$ins_file = "$out_dir/insertions/kept_" . basename($ins_file);
$del_file = "$out_dir/deletions/kept_" . basename($del_file);

print "[coovar.pl] Extracting cDNA sequence from reference on ";
print localtime()."\n";
execute("perl $prog_dir/scripts/extract-cdna.pl $exon_file $dna_file $out_dir $feature_type $feature_source");

print "[coovar.pl] Applying SNPs to cDNA sequences on ";
print localtime()."\n";
execute("perl $prog_dir/scripts/apply-snps.pl $out_dir/transcripts/reference_cdna.exons $snp_file $out_dir");

print "[coovar.pl] Categorizing SNPs on ";
print localtime()."\n";
execute("perl $prog_dir/scripts/categorize-snps.pl $out_dir/transcripts/reference_cdna.exons $out_dir/intermediate-files/target_cdna_snps.exons $prog_dir/lib/grantham_matrix  $snp_file $out_dir/intermediate-files/splice_junctions.tmp $out_dir/intermediate-files/transcripts.gff3.tmp $out_dir/intermediate-files/snps_frequency.table $out_dir");

print "[coovar.pl] Generating stats for SNP categorization on ";
print localtime()."\n";
execute("perl $prog_dir/scripts/get-stats-snps.pl $out_dir/intermediate-files/categorized_snp_coords.list $out_dir");

print "[coovar.pl] Analyzing insertions and deletions on ";
print localtime()."\n";
execute("perl $prog_dir/scripts/apply-insertions-deletions.pl $out_dir/transcripts/reference_cdna.exons $out_dir/intermediate-files/target_cdna_snps.exons $out_dir/intermediate-files/categorized_snp_coords.list $ins_file $del_file $out_dir/categorized-gvs.gvf $out_dir/intermediate-files/transcripts_snps_applied.gff3.tmp $out_dir/intermediate-files/splice_junctions.tmp $out_dir");

if($circos == 1)
{
	print "[coovar.pl] Generating CIRCOS files for SNPs, insertions, deletions and coding regions on ";
	print localtime()."\n";
	execute("perl $prog_dir/scripts/parse2circos.pl $snp_file $ins_file $del_file $out_dir/intermediate-files/sorted_gff3.tmp $dna_file $out_dir"); 
}

print "[coovar.pl] Categorized GVs can be found in $out_dir/categorized-gvs.gvf\n";
print "[coovar.pl] Categorized transcripts can be found in $out_dir/transcripts.gff3\n";

print "[coovar.pl] Done at ";
print localtime()."\n";
