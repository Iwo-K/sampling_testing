#' # Is splitting libraries for sequencing a prbolem?
#' Author: Iwo Kucinski
#' Last updated: `r Sys.time()`
#+ setup, include=FALSE
library(ggplot2)
library(dplyr)
library(gtools)
library(reshape2)
library(data.table)
source("sample_funs.R")
theme_TFpaper = theme_bw(base_size = 18, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_TFpaper)

#' Here I am looking at differences in estimating abundances (e.g. gene expression) across various sampling strategies. For instance, whether a single large sample is equivalent to a pool of two smaller samples. This is relevant to sequencing experiments, where a library can either be sequenced on a one lane or split into smaller aliquots sequences across multiple lanes.
#'
#' _In short: does it matter whether a library is sequenced once (one big run) or we use a pool of smaller runs?_

#' ## Simple sampling simulation
#' Instinctively it seems that one large sample from an infinite population should be comparable to a pool of two samples half the size. We will use a simple simulation to empirically verify this.

#' First generating 100000 random reads using a negative binomial distribution corresponding to pre-defined 10000 sequences.
reads = simreads(nuniq = 10000, n_total = 100000)

#' Distribution of observed abundances
hist(table(reads), breaks = 100)

#' I sample 10000 reads from the population and store each results. I repeat the process 1000 times to get an accurate estimate for the mean observed abundance. We store two such runs x1 and x2 for later comparisons.
x1 = sim_sampling(reads, sample_size = 10000, simno = 1000)
x2 = sim_sampling(reads, sample_size = 10000, simno = 1000)

#' Now we simulate the smaller samples, we perform 2000 simulation each sampling 5000 reads, and pooling them in pairs, to make equal to the larger sample
x = sim_sampling(reads, sample_size = 5000, simno = 2000)
x2_n2 = collapse_sampling(x, n = 2)

#' Plots below shows the differences in observed abundances:
#' 1. Difference between one large sample vs a pool of two samll samples
plot(apply(x1, 1, mean), apply(x2_n2, 1, mean))
#' 2. Difference between two random large samples
plot(apply(x1, 1, mean), apply(x2, 1, mean))

#' It looks like the estimated abundances are very similar, and the observed variance between pooled and large samples is not greater than comparison of two independent large samples.
#'
#' It seems that the empirical test confirms the intuition, but this is a bit simplistic approach.
#' In reality we start with a large library from which we take aliquots, potentially reducing the complexity of the sample. For reference 20µl of a 15nM (about 100ng of 500bp molecules) library contains around 2*10^11 molecules. A typical run on NovaSeq or HiSeq will consume...
#' So each sequencing runs consumes a ... fraction of the library
#' 
#' Additionally the distribution of fragments are highly skewed, with a small number of highly abundant ones and plenty of very low abundant fragments. Which may be relevant in situation where the complexity of the library is reduced by the first run.


#' ## Real world scenario
#' ### Loading real-world read distribution
#' In order to have a distribution of fragments approximating real libraries I load an example of RNA-seq data. These are counts per genes, which tell us the number of copies in the mixture mapping to each genes. We can operate on these gene-level estimates onwards.
realreads = fread('../sc5v3/procdata/TFnet_counts_QC.csv')
realreads = as.data.frame(realreads)
row.names(realreads) = realreads$V1
realreads = realreads[,-1]
realreads = realreads[,1:100]
realreads = rowSums(realreads)
#' Converting the counts into a vector with reads
realreads = rep(names(realreads), times = realreads)
#' Taking 10^8 reads
realreads = sample(realreads, 1e8)
#' Generating counts again
realcounts = count_reads(realreads)

#' We have this many reads (check):
sum(realcounts$count)
#' ### Simulating library splitting

#' Finally we can test a scenario mimicking a real sequencing experiment.
#' 1. We start with a library with 100 mln mln reads (reflecting molecules in a library).
#' 2. We then simulate sampling 10 mln reads from that library.
#' 3. Finally, we do a two-step sampling. First we sample 5mln reads, remove those reads from the original library (removing around 5% of the library molecules) and sample again 5mln reads, followed by pooling of the two small samples.
#' 4. Finally we repeat the process 100 times to get an estimation for the mean abundance (note: each simulation starts with a fresh, non-depleted library).

## #' Small test sample
## counts = count_reads(reads)
## z = sim_splitting(counts, ref = 10000, splits = c(5000,5000), simno =100, remove_reads = TRUE, quiet =FALSE)
## ## lapply(z, head)
## #+ fig.width = 14
## plot_split_mean(z)
## plot_split_var(z)

realsplit = sim_splitting(realcounts, ref = 1e7, splits = c(5e6, 5e6), simno = 100, remove_reads = TRUE)
realsplit = lapply(realsplit, FUN = function(x) log10(x+1))
#+ fig.width = 14
plot_split_mean(realsplit_small)
plot_split_var(realsplit_small)

#' We can also look what happens if we deplete the library more by performing the first equencing run. Above the depletion is about 5%, here we will try 25%.
#'
#' To test this, we reduce the original library size to 20 mln molecules, and perform the same sampling as above
realreads_small = sample(realreads, 2e7)
realcounts_small = count_reads(realreads_small)

realsplit_small = sim_splitting(realcounts_small, ref = 1e7, splits = c(5e6, 5e6), simno = 100, remove_reads = TRUE, quiet = TRUE)
realsplit_small = lapply(realsplit_small, FUN = function(x) log10(x+1))
#+ fig.width = 14
plot_split_mean(realsplit_small)
plot_split_var(realsplit_small)


#' Again it does not look like there is a big difference in observed abundances.
#'
