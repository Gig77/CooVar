printf "%-70s" "Fully deleted..."

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
  -e transcript.gff3 \
  -t deletion.tab \
  -r $2 \
  -o out \
  --no_contig_sum \
  | tee out/coovar.log | grep -iP "(error|warning)"

# check output
n=`grep FULLY_DELETED out/transcripts.gff3 | wc -l`
t=`grep 'Gene:Y38E10A.14' out/transcripts.gff3 | wc -l`

if [ $n != 1 ] || [ $t != 1 ] ; then
    echo "[FAILED]"
    exit 1
fi

echo "[PASSED]"
exit 0

