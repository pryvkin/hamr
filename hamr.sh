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

PROGRAM="HAMR"
VERSION="1.2.0"

function print_usage() {
  echo "USAGE: $0 [OPTIONS] reads.bam genome.fasta output_prefix" >&2
  echo "" >&2
  echo "    OPTIONS:" >&2
  echo "     Sequencing data options:" >&2
  ./hamr_cmd rnapileup --list-options
  echo "      --no-check-sorted      Don't check if BAM is sorted" >&2
  echo "" >&2
  echo "     Modification detection options:" >&2
  ./hamr_detect_mods.R --list-options
  echo "" >&2
  echo "      --help                 Print this message and exit" >&2
  echo "      --version              Print version number and exit" >&2
}

opts=()
args=()

req_args=3

# parse command line arguments
while test $# -gt 0
do
    if [[ "${1:0:1}" == "-" ]]; then
	opts=("${opts[@]}" "$1")
    else
	args=("${args[@]}" "$1")
    fi
    shift
done

# parse switches and value options
for i in "${opts[@]}"; do
    case $i in
	--help) print_usage; exit 0;;
	--version) echo "${PROGRAM} ${VERSION}"; exit 0;;
	--no-check-sorted) no_check_sorted=1;;
    esac
done

# ensure we have enough positional arguments
nargs=${#args[@]}
if [[ "${nargs}" -lt "${req_args}" ]]; then
    print_usage; exit 1;
fi

in_bam="${args[0]}"
genome_fas="${args[1]}"
outpre="${args[2]}"

seq_err=$5
hypothesis=$6
max_p=$7
max_fdr=$8

if [[ ! -e `dirname $outpre` ]]; then
    echo "Creating output directory \"`dirname $outpre`\" ..."
    mkdir -p `dirname $outpre`
fi

if [[ -z $no_check_sorted ]]; then
    echo "Checking if BAM file is sorted..." >&2
    samtools view "${in_bam}" | sort -k3,3 -k4n,4n --check=quiet
    if [[ $? -ne 0 ]]; then
	echo "ERROR: BAM file not sorted. If this is incorrect," >&2
	echo "         proceed anyway with --no-check-sorted" >&2
	exit 1
    fi
else
    echo "Skipping check of BAM sortedness..." >&2
fi

# Generate RNA pileup
echo "Computing RNA pileup..." >&2
./hamr_cmd rnapileup ${opts[@]} "${in_bam}" "${genome_fas}" \
 > ${outpre}.rnapileup

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to generate RNA pileup" >&2
    exit 1
else
    echo "Succesfully generated RNA pileup" >&2
fi

# Convert RNA pileup to BED file
echo "Converting pileup to BED" >&2
./hamr_cmd rnapileup2mismatchbed ${outpre}.rnapileup \
 > ${outpre}_mismatches.bed


# it's easier to process consecutive lines for each locus
# but order isn't guaranteed wrt strand in mismatches.bed
# so split it into + and - files for processing

# NOTE: doing this in parallel, thus requiring 2 CPUs

echo "Converting fwd strand to nuc freq table..." >&2
awk '$6 == "+"' ${outpre}_mismatches.bed | \
  ./hamr_mismatchbed2table.sh | \
  awk '$9 > 0' > \
    ${outpre}_mismatches_fwd.txt &

echo "Converting rev strand to nuc freq table..." >&2
awk '$6 == "-"' ${outpre}_mismatches.bed | \
  ./hamr_mismatchbed2table.sh | \
  awk '$9 > 0' > \
    ${outpre}_mismatches_rev.txt &

wait

echo "Sorting nuc freq table..." >&2
sort -m -k1,1 -k2n,2n -k3 \
  ${outpre}_mismatches_fwd.txt ${outpre}_mismatches_rev.txt \
  > ${outpre}_mismatches_sorted.txt

# Detect modifications using statistical testing
echo "Testing for statistical significance..." >&2
./hamr_detect_mods.R ${opts[@]} ${outpre}_mismatches_sorted.txt \
  > ${outpre}_mods.txt

if [[ $? -ne 0 ]]; then
    echo "ERROR: statistical testing failed" >&2
    exit 1
else
    echo "Statistical testing successful" >&2
fi

echo "Analysis complete." >&2
