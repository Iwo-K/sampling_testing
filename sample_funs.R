##' Simulate reads using rnbinom distribution
##'
##' .. content for \details{} ..
##' @title 
##' @param nuniq 
##' @param n_total 
##' @return a vector with reads (each type of read repeated a given number of times)
##' @author Iwo Kucinski
simreads = function(nuniq = 10000, n_total = 100000){
  #' Note: #Max nuniq 4^8

  seqs = permutations(4, 8, c("G", "A", "T", "C"), repeats.allowed = TRUE)
  seqs = apply(seqs, 1, FUN = paste0, collapse = "")
  seqs = seqs[1:nuniq]

  ns = rnbinom(nuniq, size = 4, mu = 20)
  print(hist(ns))
  print(sum(ns))
  plot(hist(ns, breaks = 100))

  seqs = rep(seqs, times = ns)

  seqs = sample(seqs, n_total)
  return(seqs)
}


##' Count reads
##'
##' .. content for \details{} ..
##' @title 
##' @param x vector with reads
##' @param names names of reads to be counted (defaults to unique(x))
##' @return data.frame with read names as row.names and counts in the count column. Order is specified by the names argument
##' @author Iwo Kucinski
count_reads = function(x, names = unique(x)){

  x = factor(x, levels = names)
  x = table(x)

  x = data.frame(row.names = names, count = as.vector(x))
  return(x)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param x vector with reads
##' @param names  names of reads to be counted (defaults to unique(x))
##' @param sample_size 
##' @param simno 
##' @return data.frame with reads as rows ad their counts in each column, each column is 1 simulation
##' @author Iwo Kucinski
sim_sampling = function(x, names = unique(x), sample_size = 10000, simno = 10){
  df = data.frame(row.names = unique(x))
  for (i in 1:simno){
    y = sample(reads, size = sample_size)
    y = count_reads(y, names = unique(x))
    colnames(y) = paste0("sim_", as.character(i))
    df = cbind(df, y)

  }
  return(df)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param x data frame with simulated sampling
##' @param n 
##' @return data.frame summarising each n simulation (sum)
##' @author Iwo Kucinski
collapse_sampling = function(x, n){

  df = data.frame(row.names = row.names(x))
  if (ncol(x) %% n != 0){
    print("Choose n that the number of column is divisible by")
    break()
  }
  if (n == ncol(x)) {
    return(x)
  }
  if (n == 1) {
    df$sims = apply(x, 1, sum)
  }
  else{
    group = lapply(as.list(0:(ncol(x)/n-1)), FUN = function(x) 1:n + x*n)
    print(group)

    counter = 1
    for (i in group){
      y = apply(x[,i], 1, sum)
      df[,paste0("sims", counter)] = y
      counter = counter+1
    }
  }
  return(df)

}
