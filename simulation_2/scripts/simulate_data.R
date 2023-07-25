# setwd("simulation_2/")

library(parallel)

# source code
source("scripts/simulate_outbreak.R", chdir = TRUE)

# color scheme
cols <- c("A" = "green", "C" = "blue", "G" = "grey50", "T" = "red")

# simulation settings
tips   <- c(100, 200, 400, 800, 1600)
size   <- 100
# factor <- c(1.5, 2, 2.5, 3)
factor <- c(1.25, 1.5, 1.75, 2, 2.5, 3)
reps   <- 1000

# mutation rate
genome_size <- 30000
mutation_rate_per_site_per_year <- 0.00084
mutation_rate_per_site_per_day  <- mutation_rate_per_site_per_year / 365
gamma <- mutation_rate_per_site_per_day

# diversification rates
r0      <- 2.5
# phi     <- 0.2
phi     <- 1 / 7
lambda0 <- r0 * phi

# overwrite?
overwrite <- FALSE
# overwrite <- TRUE

# all combinations
all_sims <- expand.grid(N = tips, L = size, factor = factor, rep = 1:reps)

invisible(mclapply(1:nrow(all_sims), function(i, ...) {
  
  # get this simulation
  this_sim <- all_sims[i,]
  this_N   <- this_sim$N
  this_L   <- this_sim$L
  this_f   <- this_sim$factor
  this_rep <- this_sim$rep
  
  # create the directory
  dir <- paste0("sims/tips_", this_N, "_f_", this_f, "/rep_", this_rep)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  cat(dir, "\n")
  
  # check if tree already exists
  if ( ("tree.nex" %in% list.files(dir) == FALSE) | overwrite == TRUE ) {
    
    # create the model
    model <- rep("A", this_L)
    
    # compute the factor
    this_factor <- this_f
    
    # simulate the data
    sim <- simulateYuleSelection2(model, lambda0, this_factor, phi, gamma, this_N, model = "multiplicative")
    phy <- ladderize(sim$phy)
    seq <- sim$seq
    seq <- seq[,sim$site, drop = FALSE]
    gain_time <- sim$gain_times[sim$site]
    tot_time  <- sim$time
    
    # make the alignment
    aln <- as.DNAbin(seq)
    
    # simulate neutral sites
    repeat {
      
      # simulate polymorphic sites
      neutral_sites <- as.DNAbin(simSeq(phy, l = genome_size, rate = gamma / 3))
      neutral_sites <- neutral_sites[,apply(neutral_sites, 2, function(x) length(unique(x))) > 1]
      
      # ensure at lest one polymorphic site
      if (ncol(neutral_sites) > 0) {
        break
      }
      
    }
    
    # make one alignment
    aln <- cbind(as.DNAbin(seq), neutral_sites)

    # write the data
    write.nexus(phy, file = paste0(dir, "/tree.nex"))
    write.nexus.data(aln, file = paste0(dir, "/seq.nex"))

    # make a plot
    pdf(paste0(dir,"/tree.pdf"))
    plot(phy, show.tip.label = FALSE, no.margin = TRUE, root.edge = TRUE)
    tiplabels(col = cols[seq[,1]], pch = 19, cex = 0.5)
    dev.off()
    
    # some summary statistics
    num_poly <- ncol(neutral_sites)
    
    # write summaries
    tab <- data.frame(total_time            = tot_time,
                      gain_time             = gain_time,
                      num_sim               = sim$num_sims,
                      num_extinct           = sim$num_extinct,
                      num_no_gains          = sim$num_no_gains,
                      num_reversions        = sim$num_reversions,
                      num_multiple_sites    = sim$num_multiple_sites,
                      num_no_extant_mutants = sim$num_no_extant_mutants,
                      not_enough_samples    = sim$not_enough_samples)
    write.table(tab, file = paste0(dir, "/sim_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    
  }
  
}, mc.cores = 6, mc.preschedule = FALSE))
