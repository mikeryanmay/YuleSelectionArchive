# setwd("simulation_1")

library(parallel)
library(TESS)
source("scripts/simulate_outbreak.R")

cols <- c("A" = "green", "C" = "blue", "G" = "grey50", "T" = "red")

# diversification parameters
r0      <- 2.5
phi     <- 1 / 7
lambda0 <- r0 * phi
factor  <- c(1.25, 1.5, 1.75, 2)

# mutation rate
genome_size <- 30000
mutation_rate_per_site_per_year <- 0.00084
mutation_rate_per_site_per_day  <- mutation_rate_per_site_per_year / 365
mu <- mutation_rate_per_site_per_day

# the observation times per size
observation_days <- c(19, 14, 12, 10)

# simulation replicates
reps <- 1000

# all combinations
all_sims <- expand.grid(factor = factor, rep = 1:reps)

# do simulation
invisible(mclapply(1:nrow(all_sims), function(i, ...) {
# invisible(lapply(1:nrow(all_sims), function(i, ...) {

  # get this simulation
  this_sim    <- all_sims[i,]
  this_factor <- this_sim$factor
  this_rep    <- this_sim$rep
  
  cat(this_factor, "\t", this_rep, "\n", sep = "")
  
  # get the lambda1
  lambda1 <- lambda0 * this_factor
  
  # make the data dir
  this_dir <- paste0("sims/factor_", this_factor, "/sim_", this_rep)
  
  # get the observation days
  observation_times <- 1:observation_days[factor == this_factor]
  
  # simulate
  sim <- simulateOutbreak(observation_times, lambda0, lambda1, phi, mu, genome_size, nmax = 1200)
  
  # save data
  for(j in 1:length(observation_times)) {
    
    # make directory
    this_day_dir <- paste0(this_dir, "/day_", observation_times[j])
    dir.create(this_day_dir, showWarnings = FALSE, recursive = TRUE)
    
    # get tree
    this_tree <- sim$trees[[j]]

    # save
    if ( nrow(this_tree$edge) == 1 ) {
      file.create(paste0(this_day_dir, "/is_singleton.txt"))
    } else {
      write.nexus(this_tree, file = paste0(this_day_dir, "/tree.nex"))      
    }
    
    # drop the root edge for analysis
    this_tree$root.edge <- NULL
    
    # get alignment
    this_aln <- this_tree$aln
    
    # save
    write.nexus.data(this_aln, file = paste0(this_day_dir, "/seq.nex"))
    
  }
  
  # write rejection summary
  df <- data.frame(sim$process_died_before_present, sim$process_exploded_before_present, sim$process_not_extant, sim$process_didnt_branch)
  write.table(df, file = paste0(this_dir, "/sim_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # write the size summary
  ntips <- sapply(sim$trees, Ntip)
  cat(ntips, file = paste0(this_dir, "/size_summary.tsv"), sep = "\t")
  cat("\n", file = paste0(this_dir, "/size_summary.tsv"), sep = "", append = TRUE)
  
  # make plots
  xlim <- c(0, max(observation_times))
  
  pdf(paste0(this_dir, "/trees.pdf"), height = 4)
  for(j in 1:length(observation_times)) {
    
    # get the sim
    this_sim <- sim$trees[[j]]
    
    if ( nrow(this_sim$edge) > 1 ) {

      # get tip data
      this_aln <- as.character(this_sim$aln[,1])
      tip_data <- toupper(this_aln[this_sim$tip.label,])
      tip_cols <- cols[tip_data]
      
      # plot the tree
      par(mar = c(2,0,0,0))
      plot(this_sim, show.tip.label = FALSE, x.lim = xlim, root.edge = TRUE)
      tiplabels(col = tip_cols, pch = 19, cex = 0.5)
      axis(1, at = c(0, observation_times), labels = c(0, observation_times), las = 2)
      
    } else {
      
      par(mar = c(2,0,0,0))
      plot(NA, xlim = xlim, ylim = c(0,1), xaxt = "n")
      axis(1, at = c(0, observation_times), labels = c(0, observation_times), las = 2)
      
    }
    
  }
  dev.off()
  
  times <- c(0, observation_times)
  tips  <- c(2, sapply(sim$trees, Ntip))
  
  pdf(paste0(this_dir, "/samples_through_time.pdf"), height = 4)
  par(mar=c(4,4,0,0) + 0.1)
  plot(times, tips, type = "s", xlab = "time", ylab = "number of samples")
  dev.off()
  
}, mc.cores = 6, mc.preschedule = FALSE))

