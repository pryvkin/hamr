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

# takes mismatch bed file on stdin and counts nucleotide
#   frequencies at each site, resulting in a 
#   nucleotide frequency table (zero-based coords)

awk '
function init() {
  counts["A"] = 0;
  counts["C"] = 0;
  counts["G"] = 0;
  counts["T"] = 0;
  counts["NR"] = 0;
}
function output() {
  print prev_chr, prev_bp, prev_str, prev_from_nuc,
     counts["A"], counts["C"], counts["G"],counts["T"],counts["NR"]
}
BEGIN {
  FS="\t"
  OFS="\t"
  init()
}
{
  loc = $1";"$2";"$6

  if (NR > 1 && (loc != prev_loc)) {
    output()
    init()
  }

  split($4,a,">")
  from_nuc = a[1]
  to_nuc = a[2]

  split($5,b,";")
  count = b[1]

  if (to_nuc == ".")
    to_nuc = from_nuc
  else
    counts["NR"] += count

  counts[to_nuc] += count

  prev_from_nuc = from_nuc
  prev_chr = $1
  prev_bp = $2
  prev_str = $6
  prev_loc = $1";"$2";"$6
}
END { output() }'

