#!/bin/bash
#  Copyright (c) 2013 University of Pennsylvania
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.


if [ $# -lt 11 ]
then
  echo "USAGE: in_bam genome_fas min_q min_coverage seq_err \
 hypothesis max_p max_fdr out_file filter_ends[0|1] not_strand_specific[0|1]" >&2
  exit 1
fi

in_bam=$1
genome_fas=$2
min_q=$3
min_coverage=$4
seq_err=$5
hypothesis=$6
max_p=$7
max_fdr=$8
outfn=$9
shift
filterends=$9
shift
not_strand_specific=$9

if [[ $hypothesis != "H1" && $hypothesis != "H4" ]]
then
  echo "hypothesis must be 'H1' or 'H4'" >&2
  exit 1
fi

if [[ $not_strand_specific = "1" ]]
then
    not_strand_specific=--noss
else
    not_strand_specific=
fi

if [[ -z $TMPDIR ]]; then
    export TMPDIR=/tmp
fi

mmbed=${TMPDIR}/$$.mismatches.bed

echo "getting mismatches..." >&2
./rnapileup $not_strand_specific $in_bam $genome_fas | \
  ./filter_pileup /dev/stdin $min_q $filterends | \
  awk '$4 >= '$min_coverage | \
  ./rnapileup2mismatchbed - > $mmbed

# it's easier to process consecutive lines for each locus
# but order isn't guaranteed wrt strand in mismatches.bed
# so split it into + and - files for processing

echo "converting fwd strand to nuc freq table..." >&2
awk '$6 == "+"' $mmbed | bash mismatchbed2table.sh | \
  awk '$9 > 0' > \
    ${mmbed}.txt.fwd &

echo "converting rev strand to nuc freq table..." >&2
awk '$6 == "-"' $mmbed | bash mismatchbed2table.sh | \
  awk '$9 > 0' > \
    ${mmbed}.txt.rev &

wait

echo "sorting nuc freq table..." >&2
sort -m -k1,1 -k2n,2n -k3 ${mmbed}.txt.fwd ${mmbed}.txt.rev > \
  ${mmbed}.txt

echo "testing for statistical significance..." >&2
Rscript detect_mods.R ${mmbed}.txt \
  $seq_err $hypothesis $max_p $max_fdr > \
  $outfn

rm $mmbed ${mmbed}.txt.fwd ${mmbed}.txt.rev ${mmbed}.txt
