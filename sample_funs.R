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

    counter = 1
    for (i in group){
      y = apply(x[,i], 1, sum)
      df[,paste0("sims", counter)] = y
      counter = counter+1
    }
  }
  return(df)

}


#' Plot multiple ggplots function
#'
#' Function which takes multiple ggpltos and plots them in a single window. Ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects).
#' #' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' Source: Cookbook R website: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' @param cols:   Number of columns in layout
#' @param layout: A matrix specifying the layout. If present, 'cols' is ignored.
#' @export
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' For splitting simulatio

##' Sample from a counttable.. content for \description{} (no empty lines) ..
##'
##' Sample from a count data.frame
##' @title 
##' @param x a data.frame with row.names and corresponding counts stored in a single column ("count")
##' @param size 
##' @param cols: Number of columns in layout
##' @param layout: A matrix specifying the layout. If present, 'cols' is ignored.
##' @return 
##' @author Iwo Kucinski
sample_counts = function(x, size = 1000){
  x = rep(row.names(x), times = x$count)
  x = sample(x, size = size)
  return(x)
}
##' Simulation splitting libraries for sequencing
##'
##' Takes in a counts dataframe, samples ref number of reads - simno number of times and then performs the same for splits but pooling each one of them. remove_reads allows removing reads from the library after each split. Outputs a list of dataframes, one for ref and one for splits (pooled together).
##' @title 
##' @param counts data.frame with names of reads in row.names and counts in count column
##' @param ref 
##' @param splits 
##' @param simno 
##' @param remove_reads 
##' @param cols: Number of columns in layout
##' @param layout: A matrix specifying the layout. If present, 'cols' is ignored.
##' @return 
##' @author Iwo Kucinski
sim_splitting = function(counts, ref = 10000, splits = c(5000,5000),
                         simno = 100,
                         remove_reads = TRUE,
                         quiet = TRUE){

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
    if (!quiet) print(sum(counts$count))
    refx = sample_counts(counts, size = ref)
    refx = count_reads(refx, names = reads_uniq)
    colnames(refx) = simname
    dfref = cbind(dfref, refx)

    counts_temp = counts
    #' Simulating the split + pool
    df = data.frame(row.names = reads_uniq)
    for (j in splits){
      splitx = sample_counts(counts_temp, size = j)
      splitx = count_reads(splitx, names = reads_uniq)
      df = cbind(df, splitx)

      if (!quiet) print(sum(counts_temp$count))
      if (remove_reads){
        counts_temp = counts - splitx[,"count"]
      }
    }

    df = data.frame(row.names = reads_uniq, count = rowSums(df))
    colnames(df) = paste0("sim_", as.character(i))
    dfsplit = cbind(dfsplit, df)
  }
  return(list(ref = dfref, split = dfsplit))
}

plot_split_mean = function(simsplit){
                                        #Generating ref x ref comparison (1/2 of simulations)
  r = simsplit$ref
  simno = ncol(r)
  r1 = r[,1:simno/2]
  r2 = r[,(simno/2):ncol(r)]
  g1 = ggplot(data.frame(), aes(x = rowMeans(r1), y = rowMeans(r2))) + geom_point(alpha = 0.5) + 
  xlab("Large sample mean log10(counts)") + ylab("Large sample mean log10(counts)")

  s = simsplit$split
  s1 = s[,1:simno/2]
  g2 = ggplot(data.frame(), aes(x = rowMeans(r1), rowMeans(s1))) + geom_point(alpha = 0.5) +
    xlab("Large sample mean log10(counts)") + ylab("1/2 sample mean log10(counts)")

  multiplot(g1, g2, cols = 2)

}

rowVars = function(x) apply(x,1, FUN = var)

plot_split_var = function(simsplit){
                                        #Generating ref x ref comparison (1/2 of simulations)
  r = simsplit$ref
  simno = ncol(r)
  r1 = r[,1:simno/2]
  r2 = r[,(simno/2):ncol(r)]
  g1 = ggplot(data.frame(), aes(x = rowVars(r1), y = rowVars(r2))) + geom_point(alpha = 0.5) + 
  xlab("Large sample variance log10(counts)") + ylab("Large sample variance log10(counts)")

  s = simsplit$split
  s1 = s[,1:simno/2]
  g2 = ggplot(data.frame(), aes(x = rowVars(r1), rowVars(s1))) + geom_point(alpha = 0.5) +
    xlab("Large sample variance log10(counts)") + ylab("1/2 sample variance log10(counts)")

  multiplot(g1, g2, cols = 2)
  }

