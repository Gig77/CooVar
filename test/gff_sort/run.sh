printf "%-70s" "GFF filtering and sorting..."

# change current directory to script directory
cd $(dirname $(readlink -f $0))

# remove and recreate output directory
if [ $1 == 0 ]; then
  if [ -d "out" ]; then
    rm -rf out/
  fi
  mkdir out
fi

# run CooVar
perl ../../coovar.pl \
  -e unsorted_transcripts.gff3 \
  -t variants.list \
  -r $2 \
  -o out \
  --no_contig_sum \
  --feature_type exon \
  --feature_source mysource \
  --circos \
  | tee out/coovar.log | grep -iP "(error|warning)"

# check output
s=`cat out/intermediate-files/sorted_gff3.tmp | wc -l`
if [ $s != 8 ] ; then
    echo "[FAILED] filtering: $s"
    exit 1
fi

s=`head -1 out/intermediate-files/sorted_gff3.tmp | grep 17486937 | wc -l`
if [ $s != 1 ] ; then
    echo "[FAILED] sorting: $s"
    exit 1
fi

s=`tail -1 out/intermediate-files/sorted_gff3.tmp | grep 17714656 | wc -l`
if [ $s != 1 ] ; then
    echo "[FAILED] sorting: $s"
    exit 1
fi

echo "[PASSED]"
exit 0
