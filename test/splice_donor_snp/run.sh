printf "%-70s" "Splice donor SNP with GFF3 input..."

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
  -e test.gff3 \
  -t test.table \
  -r $2 \
  -o out \
  --no_contig_sum \
  | tee out/coovar.log | grep -iP "(error|warning)"

# check output
n=`grep splice_donor_variant out/categorized-gvs.gvf | wc -l`
if [ $n != 1 ] ; then
    echo "[FAILED]"
    exit 1
fi

echo "[PASSED]"
exit 0
