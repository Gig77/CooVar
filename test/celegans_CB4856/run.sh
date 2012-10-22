printf "%-70s" "C. elegans hawaiian strain (CB4856)..."

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
  -e $3 \
  -t $4 \
  -r $2 \
  -o out \
  --no_contig_sum \
  | tee out/coovar.log | grep -iP "(error|warning)"

# check output
s=`grep synonymous out/categorized-gvs.gvf | wc -l`
if [ $s != 8581 ] ; then
    echo "[FAILED] synonymous: $s"
    exit 1
fi

m=`grep missense out/categorized-gvs.gvf | wc -l`
if [ $m != 9804 ] ; then
    echo "[FAILED] missense: $m"
    exit 1
fi

sl=`grep stop_lost out/categorized-gvs.gvf | wc -l`
if [ $sl != 40 ] ; then
    echo "[FAILED] stop_lost: $sl"
    exit 1
fi

sg=`grep stop_gained out/categorized-gvs.gvf | wc -l`
if [ $sg != 152 ] ; then
    echo "[FAILED] stop_gained: $sg"
    exit 1
fi

sd=`grep splice_donor out/categorized-gvs.gvf | wc -l`
if [ $sd != 55 ] ; then
    echo "[FAILED] splice_donor: $sd"
    exit 1
fi

sa=`grep splice_acceptor out/categorized-gvs.gvf | wc -l`
if [ $sa != 41 ] ; then
    echo "[FAILED] splice_acceptor: $sa"
    exit 1
fi

fs=`grep frameshift out/categorized-gvs.gvf | wc -l`
if [ $fs != 263 ] ; then
    echo "[FAILED] frameshift: $fs"
    exit 1
fi

if=`grep inframe out/categorized-gvs.gvf | wc -l`
if [ $if != 19 ] ; then
    echo "[FAILED] inframe: $if"
    exit 1
fi

oi=`grep ORF_INTACT out/transcripts.gff3 | wc -l`
if [ $oi != 15293 ] ; then
    echo "[FAILED] ORF_INTACT: $oi"
    exit 1
fi

od=`grep ORF_DISRUPTED out/transcripts.gff3 | wc -l`
if [ $od != 517 ] ; then
    echo "[FAILED] ORF_DISRUPTED: $od"
    exit 1
fi

op=`grep ORF_PRESERVED out/transcripts.gff3 | wc -l`
if [ $op != 8446 ] ; then
    echo "[FAILED] ORF_PRESERVED: $op"
    exit 1
fi

echo "[PASSED]"
exit 0
