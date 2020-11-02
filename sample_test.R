#' # Analysis of samplig of libraries - unfinished
library(ggplot2)
library(dplyr)
library(gtools)
library(reshape2)
library(data.table)


#' ## Small and simple simulation of reads
source("sample_funs.R")
reads = simreads(nuniq = 20000, n_total = 100000)
length(reads)

counts = count_reads(reads)

simtest = sim_sampling(reads, simno = 20)

coltest = collapse_sampling(simtest, n = 2)

#' Looking at libraries - 1 lane versus split across two lanes, but the original librar is not depleted
x1 = sim_sampling(reads, sample_size = 10000, simno = 1000)
x2 = sim_sampling(reads, sample_size = 10000, simno = 1000)

x = sim_sampling(reads, sample_size = 5000, simno = 2000)
x2_n2 = collapse_sampling(x, n = 2)

plot(apply(x1, 1, mean), apply(x2_n2, 1, mean))
plot(apply(x1, 1, mean), apply(x2, 1, mean))

#' ## The same with reall-world scenarion (UNIFNISHED)
#' Now the same but depleting the original library - how to deal with the simulations though?
#' For now will remove just an equivalent set (not mathing any of the simulations)
reads2 = reads
length(reads2)

x1 = sim_sampling(reads2, sample_size = 10000, simno = 1000)
#' Now removing 10% of reads (just last 100000) - note: this does not match the ones samples above
reads2 = reads[1:(length(reads2)-10000)]

x = sim_sampling(reads2, sample_size = 5000, simno = 2000)
x2_n2 = collapse_sampling(x, n = 2)

x2_n2 = x2_n2[row.names(x1),]
x1 = x1[row.names(x2_n2),]
plot(apply(x1, 1, mean), apply(x2_n2, 1, mean))
plot(apply(x1, 1, mean), apply(x2, 1, mean))


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
realtest = sim_sampling(realreads, simno = 20)

#' ## Simualating library splitting

#' Still needs to solve the no replacement things
#' I think it's negligible because the original library is massive and only a small fraction gets sequenced
sim_splitting = function(reads, ref = 10000, splits = c(5000,5000), simno = 100, remove_reads = TRUE){

  if (sum(splits) != ref){
    print("Error sum of splits does not equal the ref")
    stop()
  }
  reads_uniq = unique(reads)

  reads_freq = count_reads(reads, names = reads_uniq)
  print(head(reads_freq))
  reads_temp = rep(row.names(reads_freq), times = reads_freq$count)
  print(head(reads_temp))

  dfref = data.frame(row.names = reads_uniq)
  dfsplit = data.frame(row.names = reads_uniq)

  for (i in 1:simno){
    #' Setting a copy so that reads can be removed
    reads_temp = reads
    #' Simulating the reference
    x = sample(reads, size = ref)
    x = count_reads(x, names = reads_uniq)
    colnames(x) = paste0("sim_", as.character(i))
    dfref = cbind(dfref, x)

    #' Simulating the split + pool
    df = data.frame(row.names = reads_uniq)
    for (j in splits){
      if (remove_reads){

      }
      x = sample(reads, size = j)
      x = count_reads(x, names = reads_uniq)
      df = cbind(df, x)
    }
    df = data.frame(row.names = reads_uniq, count = rowSums(df))
    colnames(df) = paste0("sim_", as.character(i))
    dfsplit = cbind(dfsplit, df)
  }
  return(list(ref = dfref, split = dfsplit))
}
z = sim_splitting(reads, ref = 10000, splits = c(5000,5000))
plot(rowMeans(z$ref), rowMeans(z$split))

set.seed(123)
a = sample(reads, 1000)
a2 = count_reads(a, names = unique(reads))

gz = count_reads(reads)
set.seed(123)
b = sample(row.names(gz), size = 1000, prob = gz$count/sum(gz$count))
b = count_reads(b, names=  unique(reads))

head(a2)
head(b)
sum(a2$count)
sum(b$count)


gz = function(){
  start = Sys.time()
  a = sample(1:1000, size = 100,  prob = 1:1000/sum(1:1000) )

  end = Sys.time()
  print(end-start)
  }
gz()

gz = function(){
  start = Sys.time()
  a = sample(1:1000, size = 100,  prob = 1:1000/sum(1:1000) )

  end = Sys.time()
  print(end-start)
}
gz()




sim_splitting = function(reads, refsub = 10000, subs = c(5000, 1000), simno = 10){
  unique_reads = unique(reads)

  #Calculating ref
  df_ref = data.frame()
  for (i in 1:simno){
    x = sim_sampling(reads = reads, sub = refsub)
    x$simno = i
    df_ref = rbind(df_ref, x)
  }
  df_ref = group_by(df_ref, reads, sub)
  df_ref = summarise(df_ref, mean_count = mean(count))

  for (sub in subs){
    for (i in simno){
      mult = refsub/sub
      if (refsub %/% sub != 0){
        print("needs to be divisible")
        break()
      }
      for (i in 1:mult){
      x = sim_sampling(reads = reads, sub)
      x$simno = i
      }
    }
  }
  return(df_ref)
}
test2 = sim_splitting(test)




                                        #Test this out:

sample(1:22,1000,replace=TRUE,prob=c(
                                0,1,0,3,7,14,30,24,5,3,3,2,4,3,1,2,3,2,2,2,1,0
                              )

       ## library(data.table)
       ## df[, mean(z), by= (seq(nrow(df)) - 1) %% n]
