#-------------------------------------------------------------------------------
# To run CooVar tests, please download and point to the following files:
#
#   - Homo_sapiens.GRCh37.68.dna.allchr.fa: human reference genome
#   - WS210_elegans.fasta: C. elegans reference genome
#   - WS210_elegans.gff3: C. elegans gene models
#   - gvs_hawaiian_strain.list: list of GVs, tab format
#
# You can obtain these files as part of the test data sets available
# for download at http://genome.sfu.ca/projects/coovar/datasets:
#
#   - HG00732-200-37-ASM.tar.gz: includes human reference genome (GRCh37.68)
#   - Celegans_Hawaiian_CB4856.tar.gz: includes worm reference genome (WS210) 
#
#-------------------------------------------------------------------------------

HUMAN_REF=/home/cfa24/genomes/Homo_sapiens/Homo_sapiens.GRCh37.68.dna.allchr.fa
WORM_REF=/home/cfa24/va_test/Celegans_Hawaiian_CB4856/WS210_elegans.fasta
WORM_GFF=/home/cfa24/va_test/Celegans_Hawaiian_CB4856/WS210_elegans.gff3
WORM_GVS=/home/cfa24/va_test/Celegans_Hawaiian_CB4856/gvs_hawaiian_strain.list

KEEP_OLD_OUTPUT=0

sh fully_deleted/run.sh $KEEP_OLD_OUTPUT $WORM_REF
sh already_disrupted/run.sh $KEEP_OLD_OUTPUT $HUMAN_REF
sh missense_snp/run.sh $KEEP_OLD_OUTPUT $HUMAN_REF
sh synonymous_snp/run.sh $KEEP_OLD_OUTPUT $HUMAN_REF
sh stop_lost/run.sh $KEEP_OLD_OUTPUT $HUMAN_REF
sh stop_gained/run.sh $KEEP_OLD_OUTPUT $HUMAN_REF
sh splice_junction_insertion/run.sh $KEEP_OLD_OUTPUT $HUMAN_REF
sh splice_junction_deletion/run.sh $KEEP_OLD_OUTPUT $HUMAN_REF
sh complex_deletion/run.sh $KEEP_OLD_OUTPUT $HUMAN_REF
sh inframe_insertion/run.sh $KEEP_OLD_OUTPUT $HUMAN_REF
sh celegans_CB4856/run.sh $KEEP_OLD_OUTPUT $WORM_REF $WORM_GFF $WORM_GVS


