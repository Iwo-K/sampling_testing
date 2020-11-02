#' # Analysis of samplig of libraries - unfinished
library(ggplot2)
library(dplyr)
library(gtools)
library(reshape2)
library(data.table)


source("sample_funs.R")
#' ## Generating random reads (nbinom)
reads = simreads(nuniq = 20000, n_total = 100000)

#' ## Loading real world scenario
#' Loading some example counts data and pooling to create one library
realreads = fread('../sc5v3/procdata/TFnet_counts_QC.csv')
realreads = as.data.frame(realreads)
row.names(realreads) = realreads$V1
realreads = realreads[,-1]
realreads = realreads[,1:100]
realreads = rowSums(realreads)

realreads = rep(names(realreads), times = realreads)
realcounts = count_reads(realreads)

#' ## Small and simple simulation of reads

#' ## Simualating library splitting

#' Still needs to solve the no replacement things
#' I think it's negligible because the original library is massive and only a small fraction gets sequenced
sample_counts = function(x, size = 1000){
  x = rep(row.names(x), times = x$count)
  x = sample(x, size = size)
  return(x)
}

sim_splitting = function(counts, ref = 10000, splits = c(5000,5000),
                         simno = 100,
                         remove_reads = TRUE){

  if (sum(splits) != ref){
    print("Error sum of splits does not equal the ref")
    stop()
  }
  reads_uniq = row.names(counts)

  dfref = data.frame(row.names = reads_uniq)
  dfsplit = data.frame(row.names = reads_uniq)

  for (i in 1:simno){
    simname = paste0("sim_", as.character(i))
    #' Simulating the reference
    print(sum(counts$count))
    refx = sample_counts(counts, size = ref)
    refx = count_reads(refx, names = reads_uniq)
    colnames(refx) = simname
    dfref = cbind(dfref, refx)

    counts_temp = counts
    #' Simulating the split + pool
    df = data.frame(row.names = reads_uniq)
    for (j in splits){
      print(sum(counts_temp$count))
      #This part needs fixing so that it removes the previous splits!
      ## if (remove_reads){
      ##   counts_temp = counts - refx[,simname]
      ## }
      splitx = sample_counts(counts_temp, size = j)
      splitx = count_reads(splitx, names = reads_uniq)
      df = cbind(df, splitx)
    }
    df = data.frame(row.names = reads_uniq, count = rowSums(df))
    colnames(df) = paste0("sim_", as.character(i))
    dfsplit = cbind(dfsplit, df)
  }
  return(list(ref = dfref, split = dfsplit))
}
counts = count_reads(reads)
z = sim_splitting(counts, ref = 10000, splits = c(5000,5000), simno =100, remove_reads = TRUE)
lapply(z, head)

plot(rowMeans(z$ref), rowMeans(z$split))

plot_split = function(simsplit){
  #Generating ref x ref comparison (1/2 of simulations)
  r = simsplit$ref
  simno = ncol(r)
  r1 = r[,1:simno/2]
  r2 = r[,(simno/2):ncol(r)]
  g1 = ggplot(data.frame(), aes(x = rowMeans(r1), y = rowMeans(r2))) + geom_point(alpha = 0.5)

  s = simsplit$split
  s1 = s[,1:simno/2]
  g2 = ggplot(data.frame(), aes(x = rowMeans(r1), rowMeans(s1))) + geom_point(alpha = 0.5)

  multiplot(g1, g2, cols = 2)

}

plot_split(z)

sum(realcounts$count)

realsplit = sim_splitting(realcounts, ref = 1e7, splits = c(5e6, 5e6), simno = 100, remove_reads = TRUE)
realsplit = lapply(realsplit, FUN = function(x) log10(x+1))
plot_split(realsplit)

#' To finish off plot this against difference fractions of library that are being sequences (should be more variable if we are sequencing a large chunk etc)
