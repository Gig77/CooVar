printf "%-70s" "Missense SNP..."

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
  -e ENST00000295641.gtf \
  -v HG00732-200-37-ASM.missense.vcf \
  -r $2 \
  -o out \
  --no_contig_sum \
  | tee out/coovar.log | grep -iP "(error|warning)"

# check output
n=`grep non_conservative_missense_codon out/categorized-gvs.gvf | wc -l`
m=`grep 'tct>tTt_S>F_aa752_codon_loc2_RADICAL(155),tct>tTt_S>F_aa741_codon_loc2_RADICAL(155)' out/categorized-gvs.gvf | wc -l`

if [ $n != 1 ] || [ $m != 1 ] ; then
    echo "[FAILED]"
    exit 1
fi

echo "[PASSED]"
exit 0
