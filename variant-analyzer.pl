#!/usr/bin/perl

use strict;
use Getopt::Long;

my $exon_file;
my $dna_file;
my $snp_file;
my $ins_file;
my $del_file;
my $tab_file;
my $vcf_file;

my $circos='';

GetOptions ("exon|e=s" => \$exon_file,
		"ref|r=s"  => \$dna_file,
		"tab|t=s"  => \$tab_file,
		"vcf|v=s"  => \$vcf_file,
		"circos"   => \$circos
		);	

unless ($exon_file && $dna_file && ($vcf_file || $tab_file))
{
    print "To run\:\n \.\/variant-analyzer.pl -e EXONS -r REFERENCE \(-t GVS_TAB_FORMAT \| -v GVS_VCF_FORMAT\) [\-\-circos]\n";
    exit;
}

print "Parsing GV files on ";
system(date);
if(defined $tab_file)
{
	system("./scripts/TAB2VA.pl $tab_file");
	$snp_file = $tab_file . '.VA_SNP';
	$ins_file = $tab_file . '.VA_INS';
	$del_file = $tab_file . '.VA_DEL';
}

if(defined $vcf_file)
{
	system("./scripts/VCF2VA.pl $vcf_file");
        $snp_file = $vcf_file .	'.VA_SNP';
        $ins_file = $vcf_file .	'.VA_INS';
        $del_file = $vcf_file .	'.VA_DEL';

}
print "Done at ";
system(date);

unless (-d "./VA_Intermediate_Files")
{
        system("mkdir VA_Intermediate_Files");
}

print "Checking the consistency of GV files on ";
system(date);
system("./scripts/revise-GV-files.pl $snp_file $ins_file $del_file");
print "Done at ";
system(date);

my $ex_snps = "excluded_" . $snp_file;
my $ex_ins = "excluded_" . $ins_file;
my $ex_del = "excluded_" . $del_file;

system("mv $snp_file $ins_file $del_file VA_Intermediate_Files");

$snp_file = "kept_" . $snp_file;
$ins_file = "kept_" . $ins_file;
$del_file = "kept_" . $del_file;

unless (-d "./VA_SNPs")
{
        system("mkdir VA_SNPs");
}

unless (-d "./VA_Insertions")
{
        system("mkdir VA_Insertions");
}

unless (-d "./VA_Deletions")
{
        system("mkdir VA_Deletions");
}
unless (-d "./VA_Transcripts")
{
        system("mkdir VA_Transcripts");
}

system("mv $snp_file $ex_snps VA_SNPs");
system("mv $ins_file $ex_ins VA_Insertions");
system("mv $del_file $ex_del VA_Deletions");

print "Extracting cDNA sequence from reference on ";
system(date);
system("./scripts/extract-cDNA.pl $exon_file $dna_file");
system("mv reference_cDNA.exons reference_cDNA.fasta reference_peptides.fasta VA_Transcripts");
print "Done at ";
system(date);

print "Applying SNPs to cDNA sequences on ";
system(date);
system("./scripts/apply-SNPs.pl VA_Transcripts/reference_cDNA.exons VA_SNPs/$snp_file");
print "Done at ";
system(date);

print "Categorizing SNPs on ";
system(date); 
system("./scripts/categorize-SNPs.pl VA_Transcripts/reference_cDNA.exons target_cDNA_SNPs.exons ./lib/grantham_matrix  VA_SNPs/$snp_file splice_junctions.tmp transcripts.gff3.tmp > snps_frequency.table");
print "Done at ";
system(date);

print "Generating stats for SNP categorization on ";
system(date);
system("./scripts/get-stats-SNPs.pl categorized_snp_coords.list");
system("mv categorized_snp_coords.list snps_frequency.table codon_bias.stat ambiguous_SNPs.list categorized_SNPs.stat VA_SNPs");
print "Done at ";
system(date);

print "Analyzing insertions and deletions on ";
system(date);
system("./scripts/apply-insertions-deletions.pl VA_Transcripts/reference_cDNA.exons target_cDNA_SNPs.exons VA_SNPs/categorized_snp_coords.list VA_Insertions/$ins_file VA_Deletions/$del_file VA_categorized_GVs.gvf transcripts_snps_applied.gff3.tmp");
system("mv variant_peptides.fasta variant_cDNA.exons variant_cDNA.fasta reference_variant.alignment variant.stat VA_Transcripts");
system("mv variant.summary ./VA_Intermediate_Files");
system("mv distribution_insertions.out ./VA_Insertions");
system("mv distribution_deletions.out ./VA_Deletions");

system("mv transcripts.gff3.tmp target_cDNA_SNPs.exons target_cDNA_SNPs.fasta splice_junctions.tmp transcripts_snps_applied.gff3.tmp VA_SNPs/snps_frequency.table VA_SNPs/categorized_snp_coords.list VA_Intermediate_Files");

print "Analysis of insertions and deletions done on ";
system("date");

if($circos == 1)
{
	print "Generating CIRCOS files for SNPs, insertions, deletions and coding regions on ";
	system(date);
	system("./scripts/parse2circos.pl VA_SNPs/$snp_file VA_Insertions/$ins_file VA_Deletions/$del_file $exon_file $dna_file");
	my $exon_circos_out = $exon_file . '.circos';
	system("mv $exon_circos_out ./VA_Transcripts");
	print "Generation of CIRCOS input files done on ";
	system(date);
}
