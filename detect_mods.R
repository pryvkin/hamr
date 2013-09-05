#!/usr/bin/Rscript env
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

if (length(commandArgs(T)) < 5) {
  cat("USAGE: detect_mods.R in_table seq_err hyp max_p max_fdr\n")
  q()
}

infn = commandArgs(T)[1]
seq.err = as.numeric(commandArgs(T)[2])
hypothesis = commandArgs(T)[3]
maxp = as.numeric(commandArgs(T)[4])
maxq = as.numeric(commandArgs(T)[5])

nucs = c('A','C','G','T')

x = read.table(infn, as.is=T, sep="\t",
  col.names=c('chr', 'bp', 'strand',
    'refnuc', 'A', 'C', 'G', 'T', 'nonref'))


x$ref = rowSums(x[,nucs]) - x$nonref

hyps = c('AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT')
hyp.ps = array(NA, dim=c(nrow(x), length(hyps)))
colnames(hyp.ps) = hyps

for(h in hyps) {
  correct.nucs = unique(unlist(strsplit(h,'')))
  err.nucs = nucs[!nucs %in% correct.nucs]
  correct.counts = rowSums(data.frame(x[,correct.nucs]))
  err.counts = rowSums(data.frame(x[,err.nucs]))
  hyp.ps[,h] = pbinom(correct.counts,
                        rowSums(cbind(correct.counts, err.counts)),
                        1-seq.err,
                        lower.tail=T)
}

hyp.union.ps = array(NA, dim=c(nrow(x), 2))
colnames(hyp.union.ps) = paste("H", c(1,4), sep='')

# RR
hyp.union.ps[,'H1'] = sapply(1:nrow(x), function(i) {
  hyp.ps[i,sprintf("%s%s",x$refnuc[i],x$refnuc[i])]
})

# all
hyp.union.ps[,'H4'] = apply(hyp.ps, 1, max)

hyp.union.ps.adj = apply(hyp.union.ps, 2, p.adjust, method='BH')

x$h1.p = hyp.union.ps[,'H1']
x$h1.padj = hyp.union.ps.adj[,'H1']
x$h4.p = hyp.union.ps[,'H4']
x$h4.padj = hyp.union.ps.adj[,'H4']

x$sig = hyp.union.ps[,hypothesis] < maxp &
  hyp.union.ps.adj[,hypothesis] < maxq

write.table(x, file="", row.names=F, col.names=T, quote=F,sep="\t")



