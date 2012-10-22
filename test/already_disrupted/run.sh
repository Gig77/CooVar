printf "%-70s" "ORF already disrupted in reference..."

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
  -e ENST00000426379.gtf \
  -v ENST00000426379.vcf \
  -r $2 \
  -o out \
  --no_contig_sum \
  | tee out/coovar.log | grep -iP "(error|warning)"

# check output
n=`grep ORF_PRESERVED out/transcripts.gff3 | wc -l`
if [[ $n != 1 ]] ; then
    echo "[FAILED]"
    exit 1
fi

echo "[PASSED]"
exit 0
