printf "%-70s" "Stop lost..."

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
  -e Homo_sapiens.GRCh37.68.chrom.stop_lost.gtf \
  -v HG00731-200-37-ASM.stop_lost.vcf \
  -r $2 \
  -o out \
  --no_contig_sum \
  | tee out/coovar.log | grep -iP "(error|warning)"

# check output
n=`grep stop_lost out/categorized-gvs.gvf | wc -l`
if [ $n != 1 ] ; then
    echo "[FAILED]"
    exit 1
fi

echo "[PASSED]"
exit 0
