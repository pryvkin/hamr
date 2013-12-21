#!/usr/bin/env Rscript
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

### detect_mods.R
## Detect modifications using nucleotide frequencies
# Input:   - Nucleotide frequency table
#          - Sequencing error rate
#          - Hypothesis type (H1 or H4)
#          - Maximum p-value
#          - Maximum FDR-adjusted p-value
#
# Output:  - Table containing results for each site

logWrite <- function(s) { write(s, stderr()) }

print_usage <- function(list.opts.only=F) {
    if (!list.opts.only) {
        logWrite("USAGE: ./hamr_detect_mods.R [OPTIONS] in.nuc_freq_table")
        logWrite("")
        logWrite("    OPTIONS:")
    }
    logWrite("      --hypothesis=H1|H4     H1: loose, allow SNP/edit-like")
    logWrite("                               H4: strict, only mod-like (default)")
    logWrite("      --max-p=P              Use unadj. p-value cutoff P (1.0)")
    logWrite("      --max-q=Q              Use FDR-controlled cutoff Q (0.05)")
    logWrite("      --seq-error-rate       Assumed rate of seq. errors (0.01)")
}

# parse command line arguments
req.positional.args <- 1
args <- commandArgs(T)
is.opt.arg <- grepl("^-", args)
option.args <- args[is.opt.arg]
positional.args <- args[!is.opt.arg]

option.args <- strsplit(option.args, "=")

hypothesis <- "H4"
maxp <- 1.0
maxq <- 0.05
seq.err <- 0.01

for(option.arg in option.args) {
    opt.key <- option.arg[1]
    opt.value <- option.arg[2]
    if (opt.key == "--list-options") {
        print_usage(list.opts.only=T)
        quit(status=0)
    } else if (opt.key == "--hypothesis") {
        hypothesis <- opt.value
        if ( !(hypothesis %in% c("H1", "H4"))) {
            logWrite(sprintf("ERROR: invalid null hypothesis (%s): must be H1 or H4",
                             opt.value))
            quit(status=1)
        }
    } else if (opt.key == "--max-p") {
        maxp <- as.numeric(opt.value)
        if (maxp < 0 || maxp > 1) {
            logWrite(sprintf("ERROR: invalid p-value threshold (%s): must be real number in [0,1]",
                             opt.value))
            quit(status=1)
        }
    } else if (opt.key == "--seq-error-rate") {
        seq.err <- as.numeric(opt.value)
        if (seq.err < 0 || seq.err > 1) {
            logWrite(sprintf("ERROR: invalid sequencing err. rate (%s): must be real number in [0,1]",
                             opt.value))
            quit(status=1)
        }
    }
}

if (length(positional.args) < req.positional.args) {
    print_usage()
    quit(status=1)
}


infn <- positional.args[1]

logWrite(sprintf("  Using null hypothesis %s", hypothesis))
logWrite(sprintf("  Assuming sequencing error rate is %f", seq.err))
logWrite(sprintf("  Using p-value threshold %f", maxp))
logWrite(sprintf("  Using FDR adj. p-value threshold %f", maxq))

nucs = c('A','C','G','T')

tryCatch( {
    x = read.table(infn, as.is=T, sep="\t",
        col.names=c('chr', 'bp', 'strand',
            'refnuc', 'A', 'C', 'G', 'T', 'nonref'))
}, error = function(cond) {
    logWrite(sprintf("ERROR: Could not open file %s", infn))
    message(cond)
    quit(status=1)
} )

if (nrow(x) < 1) {
    logWrite("WARNING: empty input to detect_mods (no mismatches found?)")
    quit(status=1)
}

# count reference nucleotide observations at each site
x$ref = rowSums(x[,nucs]) - x$nonref

hyps = c('AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT')
hyp.ps = array(NA, dim=c(nrow(x), length(hyps)))
colnames(hyp.ps) = hyps

# compute p-value for each possible genotype
for(h in hyps) {
  correct.nucs = unique(unlist(strsplit(h,'')))
  err.nucs = nucs[!nucs %in% correct.nucs]
  correct.counts = rowSums(data.frame(x[,correct.nucs]))
  err.counts = rowSums(data.frame(x[,err.nucs]))
  # binomial distribution
  hyp.ps[,h] = pbinom(correct.counts,
                        rowSums(cbind(correct.counts, err.counts)),
                        1-seq.err,
                        lower.tail=T)
}

hyp.union.ps = array(NA, dim=c(nrow(x), 2))
colnames(hyp.union.ps) = paste("H", c(1,4), sep='')

# H0_1: null hypothesis is: genotype = homozygous reference
# RR (ref nuc, ref nuc)
hyp.union.ps[,'H1'] = sapply(1:nrow(x), function(i) {
  hyp.ps[i,sprintf("%s%s",x$refnuc[i],x$refnuc[i])]
})

# H0_4: null hypothesis is: genotype = any one or two alelle(s)
# take maximum p-value across all genotype hypotheses
hyp.union.ps[,'H4'] = apply(hyp.ps, 1, max)

# adjust p-values
hyp.union.ps.adj = apply(hyp.union.ps, 2, p.adjust, method='BH')

x$h1.p = hyp.union.ps[,'H1']
x$h1.padj = hyp.union.ps.adj[,'H1']
x$h4.p = hyp.union.ps[,'H4']
x$h4.padj = hyp.union.ps.adj[,'H4']

x$sig = hyp.union.ps[,hypothesis] < maxp &
  hyp.union.ps.adj[,hypothesis] < maxq

write.table(x, file="", row.names=F, col.names=T, quote=F,sep="\t")

