# setwd("simulation_1")

library(matrixStats)
library(ape)
library(TreeTools)
library(RColorBrewer)
library(ggplot2)
library(viridis)

# specify figure directory
figdir <- "figures/"

# enumerate the analyses
factor <- c(1.25, 1.5, 1.75, 2)
reps   <- 1000

# make all combinations
grid <- expand.grid(factor = factor, rep = 1:reps, stringsAsFactors = FALSE)

# compute all the summaries
summaries <- do.call(rbind, lapply(1:nrow(grid), function(i) {
  
  # get the analysis
  this_grid   <- grid[i,]
  this_factor <- this_grid$factor
  this_rep    <- this_grid$rep
  
  cat(i, " -- ", nrow(grid),"\n", sep = "")
  
  # check if the output file exists
  this_tsv <- paste0("sims/factor_", this_factor, "/sim_", this_rep, "/sim_summary.tsv")
  if ( file.exists(this_tsv) == FALSE ) {
    return(NULL)
  }
  
  # read output file
  results <- read.table(this_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  results <- cbind(data.frame(factor = this_factor, rep = this_rep), results)
  
  return(results)
  
}))

# summaries <- cbind(grid, summaries)

# plot summaries
colors <- gsub("FF", "", turbo(5, begin = 0.1, end = 0.95))
w <- 0.75

layout_mat <- matrix(1:(length(factor)), ncol = length(factor), nrow = 1, byrow = TRUE)

pdf(paste0(figdir, "rejection_summary.pdf"), height = 3, width = 7)
layout(layout_mat)
par(mar=c(0,0,0,0), oma = c(6,4.5,2,0) + 0.1, lend = 2)

for(i in 1:length(factor)) {
  
  # get the factor
  this_f <- factor[i]
  
  # get the relevant samples
  these_summaries <- summaries[summaries$factor == this_f,]

  # get the rejection counts
  these_summaries <- these_summaries[,-c(1:2)]
  
  # add the accepted simulations
  these_summaries$sim.accepted <- 1
  
  # compute the total number of simulations
  total_simulations <- sum(these_summaries)

  # compute sums
  sums <- colSums(these_summaries)
    
  # compute proportions
  props <- sums / total_simulations
  
  # sort
  props <- props[c("sim.accepted", "sim.process_died_before_present", "sim.process_exploded_before_present")]
  names(props) <- c("accepted", "extinct", "exploded")
 
  barplot(props, ylim = c(0,1), border = NA, col = colors, las = 2, yaxt = "n")
  box()
  abline(h = 0.05, lty = 2)
  if ( i == 1 ) {
    axis(side = 2, lwd = 0, lwd.tick = 1, las = 2)
    mtext(side = 2, text = "frequency", line = 2.5)
  }
  mtext(bquote(delta == .(this_f)), side = 3, line = 0.5)
   
}
dev.off()
















