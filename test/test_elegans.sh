#----------------------------------------
# CooVar test run on a C. elegans dataset
# this dataset can be obtained from http://genome.sfu.ca/projects/coovar/elegans.tar.gz
# please change INPUTDIR and OUTPUTDIR after extracting the tarball and before running this script
# this script is meant to be run from your CooVar installation directory, using 'sh test/test_elegans.sh'
#----------------------------------------
 
INPUTDIR=~/va_test/ELEGANS_INPUT
OUTPUTDIR=~/va_test/ELEGANS_OUTPUT

perl coovar.pl \
	-e $INPUTDIR/WS210_elegans.gff3 \
	-r $INPUTDIR/WS210_elegans.fasta \
	-t $INPUTDIR/GVs_hawaiian_strain.list \
	-o $OUTPUTDIR \
	--circos