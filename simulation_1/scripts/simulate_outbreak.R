library(ape)
library(expm)

# we simulate a lineage with a single site under selection
# we assumed the new allele arises at time t = 0
# we then sample at specified time points after that time
simulateOutbreak <- function(observation_times, lambda0, lambda1, phi, mu, genome_size, nmax = 10000, ladderize = TRUE) {
  
  # condition on making it to the present, and not getting too large
  process_died_before_present <- 0
  process_exploded_before_present <- 0
  process_not_extant <- 0
  process_didnt_branch <- 0
  
  repeat {

    # simulate
    sim <- simulateOutbreakConditional(observation_times, lambda0, lambda1, phi, mu, genome_size, nmax)
    
    if ( inherits(sim, "numeric") ) {
      
      if ( sim == 1 ) { # processes died out
        process_died_before_present <- process_died_before_present + 1
      } else if ( sim == 2 ) { # processes exploded
        process_exploded_before_present <- process_exploded_before_present + 1
      }
      
    } else {
      
      # check that the derived variant persists
      tmp  <- sim$lineages[[length(sim$lineages)]]
      frac <- mean(tmp$selected_state[tmp$status == "extant"] == "A")
      if ( frac == 0 ) {
        process_not_extant <- process_not_extant + 1
      } else {
        break
      }
      
    }
    
  }
  
  # make trees and alignments
  phy <- vector("list", length(observation_times))
  for(i in 1:length(observation_times)) {
    phy[[i]] <- makeTreeAndAlignment(sim$lineages[[i]], sim$genomes[[i]], ladderize)
  }
  
  # bundle stuff together
  sim <- list(trees                           = phy,
              offset                          = sim$offset,
              process_died_before_present     = process_died_before_present,
              process_exploded_before_present = process_exploded_before_present,
              process_not_extant              = process_not_extant,
              process_didnt_branch            = process_didnt_branch)

  return(sim)
    
}

makeTreeAndAlignment <- function(lineages, genomes, ladderize = TRUE) {
  
  if ( nrow(lineages) == 1 ) {

    # make a singleton tree    
    edge        <- matrix(c(2, 1), nrow = 1)
    edge.length <- (lineages$end - lineages$start)
    tip.label   <- "t_1"
    Nnode       <- 1
    phy         <- list(edge        = edge, 
                        edge.length = edge.length,
                        tip.label   = tip.label,
                        root.edge   = edge.length,
                        Nnode       = Nnode)
    class(phy) <- "phylo"
    
    # write.nexus(phy, file = "test.nex")
    # pp <- read.nexus("test.nex")
    
  } else {

    # create the edge matrix
    tmp_edge   <- as.matrix(lineages[-1,1:2])
    old_tips   <- lineages$desc[lineages$status %in% c("sampled", "extant")]  
    old_nodes  <- sort(unique(lineages$anc[-1]))
    new_tips   <- 1:length(old_tips)
    new_nodes  <- 1:length(old_nodes) + length(old_tips)
    new_labels <- c(new_tips, new_nodes)
    names(new_labels) <- c(old_tips, old_nodes)
    
    edge <- cbind(new_labels[as.character(tmp_edge[,1])], new_labels[as.character(tmp_edge[,2])])
    rownames(edge) <- NULL
    
    # create the phylo object
    edge.length <- (lineages$end - lineages$start)[-1]
    tip.label   <- paste0("t_", 1:length(new_tips))
    root.edge   <- (lineages$end - lineages$start)[1]
    Nnode       <- length(tip.label) - 1
    phy         <- list(edge        = edge, 
                        edge.length = edge.length,
                        tip.label   = tip.label,
                        root.edge   = root.edge,
                        Nnode       = Nnode)
    class(phy) <- "phylo"
    
  }
  
  # create the alignment object
  neutral_genomes <- do.call(rbind, genomes)
  neutral_genomes <- neutral_genomes[lineages$status %in% c("sampled", "extant"), , drop = FALSE]
  neutral_genomes_poly <- neutral_genomes[,apply(neutral_genomes, 2, function(x) length(unique(x))) > 1, drop = FALSE]
  aln <- cbind(lineages$selected_state[lineages$status %in% c("sampled", "extant")], neutral_genomes_poly)
  rownames(aln) <- tip.label
  aln <- as.DNAbin(aln)

  # ladderize
  if (ladderize & length(tip.label) > 2) {
    phy <- ladderize(phy)
    aln <- aln[phy$tip.label,]
  }

  # attach alignment to tree
  phy$aln <- aln
  
  return(phy)
    
}

simulateOutbreakConditional <- function(observation_times, lambda0, lambda1, phi, mu, genome_size, nmax) {
  
  # helpers
  nucleotides <- c("A","C","G","T")
  neutral_genome_size <- genome_size - 1

  # make a Q matrix
  Qmat <- matrix(mu / 3, 4, 4)
  diag(Qmat) <- -mu
  colnames(Qmat) <- rownames(Qmat) <- nucleotides
  
  #############################
  # initialize the simulation #
  #############################
  
  # simulate the initial genome
  initial_neutral_genome <- sample(nucleotides, size = neutral_genome_size, replace = TRUE)
  
  # make the initial tables for the two subtrees
  outbreak_init <- data.frame(anc = 1, desc = 2, 
                              start_time = 0, end_time = NA, 
                              selected_state = "A",
                              lambda = lambda1,
                              num_mutations = 0,
                              status = "active", stringsAsFactors = FALSE)
  
  
  # make the lineage table
  lineages <- data.frame(anc = rep(NA, nmax), desc = NA,
                         start_time = NA, end_time = NA,
                         selected_state = NA,
                         lambda = NA,
                         num_mutations = NA,
                         status = NA, stringsAsFactors = FALSE)
  lineages[1,] <- outbreak_init
  dimmax <- nmax
  
  # make the genome list
  genome_list <- vector("list", nmax)
  genome_list[[1]] <- initial_neutral_genome
  
  ###########################################
  # start simulations for observation times #
  ###########################################

  # set the current time
  current_time <- 0
  lineage_index_ctr <- 2
  num_samples <- 0
  
  # make containers for each sample time
  num_observation_times <- length(observation_times)
  
  # make containers for each observation time
  observed_lineages <- vector("list", num_observation_times)
  observed_genomes  <- vector("list", num_observation_times)
  
  # simulate each observation time
  for(i in 1:num_observation_times) {
    
    # get the observation time
    this_observation_time <- observation_times[i]
    
    # cat("simulating time ", this_observation_time, "\n", sep = "")
    
    # simulate forward
    repeat {
      
      # find all the active lineages
      active_lineage_index <- which(lineages$status == "active")
      num_active_lineages  <- length(active_lineage_index)
      
      # make sure there are lineages that events can happen to!
      if ( num_active_lineages == 0 ) {
        return(1)
      }
      
      # make sure we are not too many samples
      if ( num_active_lineages > nmax) {
        return(2)
      }
      
      # get the speciation rates
      active_lineage_speciation_rates <- lineages$lambda[active_lineage_index]
      
      # get the event rate per lineage
      active_lineage_event_rate <- active_lineage_speciation_rates + phi + mu * genome_size
      
      # simulate a waiting time per lineage
      waiting_times <- rexp(num_active_lineages, rate = active_lineage_event_rate)
      min_wait      <- min(waiting_times)
      current_time  <- current_time + min_wait
      
      # terminate if we exceed the time
      if ( current_time > this_observation_time ) {
        current_time <- this_observation_time
        break
      }
      
      # otherwise, determine which lineage is affected
      event_lineage_index <- active_lineage_index[which.min(waiting_times)]
      this_lambda         <- lineages$lambda[event_lineage_index]
      
      # determine the relative probabilities of events for this lineage
      event_probs <- c(this_lambda, phi, mu, mu * neutral_genome_size) / sum(c(this_lambda, phi, mu, mu * neutral_genome_size))
      
      # choose the type of event
      event_type <- sample.int(4, prob = event_probs, size = 1)
      
      # do the event
      if ( event_type == 1 ) { # birth event
        
        if ( lineage_index_ctr + 3 > dimmax ) {
          # expand the lineage matrix
          new_lineages <- data.frame(anc = rep(NA, nmax), desc = NA,
                                     start_time = NA, end_time = NA,
                                     selected_state = NA,
                                     lambda = NA,
                                     num_mutations = NA,
                                     status = NA, stringsAsFactors = FALSE)
          lineages <- rbind(lineages, new_lineages)
          dimmax <- dimmax + nmax
        }
        
        # finalize the affected lineage
        lineages$end_time[event_lineage_index] <- current_time
        lineages$status[event_lineage_index]   <- "inactive"
        
        # make some new lineages
        new_lineage_index <- lineage_index_ctr + 0:1
        lineages[new_lineage_index,] <- data.frame(anc = lineages$desc[event_lineage_index], desc = 1:2 + lineage_index_ctr ,
                                                   start_time = current_time, end_time = NA,
                                                   selected_state = lineages$selected_state[event_lineage_index],
                                                   lambda = lineages$lambda[event_lineage_index],
                                                   num_mutations = lineages$num_mutations[event_lineage_index],
                                                   status = "active")
        
        # copy the genomes
        genome_list[[new_lineage_index[1]]] <- genome_list[[event_lineage_index]]
        genome_list[[new_lineage_index[2]]] <- genome_list[[event_lineage_index]]
        
        # increment lineage counter
        lineage_index_ctr <- lineage_index_ctr + 2
        
      } else if ( event_type == 2 ) { # sample event
        
        # finalize the affected lineage
        lineages$end_time[event_lineage_index] <- current_time
        lineages$status[event_lineage_index]   <- "sampled"
        num_samples <- num_samples + 1
        
      } else if ( event_type == 3 ) { # mutation at selected site
        
        # get the current state
        current_state <- lineages$selected_state[event_lineage_index]
        
        # choose the new state
        new_state <- sample(setdiff(nucleotides, current_state), size = 1)
        
        # update the lineage
        lineages$selected_state[event_lineage_index] <- new_state
        if ( new_state == "A" ) {
          
          # keep track of the number of gains
          num_gains <- num_gains + 1
          
          # update the speciation rate
          lineages$lambda[event_lineage_index] <- lambda1
          
        } else {
          
          # update the speciation rate
          lineages$lambda[event_lineage_index] <- lambda0
          
        }
        
        # update the number of mutations
        lineages$num_mutations[event_lineage_index] <- lineages$num_mutations[event_lineage_index] + 1
        
      } else if ( event_type == 4 ) { # neutral mutation
        
        # get the genome
        this_genome <- genome_list[[event_lineage_index]]
        
        # choose a site
        this_site <- sample.int(neutral_genome_size, size = 1)
        
        # get the current state
        current_state <- this_genome[this_site]
        
        # choose the new state
        new_state <- sample(setdiff(nucleotides, current_state), size = 1)
        
        # update the genome
        this_genome[this_site] <- new_state
        
        # store the genome
        genome_list[[event_lineage_index]] <- this_genome
        
      }
      
    }
    
    # finalize this sample time
    these_rows <- 1:(lineage_index_ctr - 1)
    these_lineages <- lineages[these_rows,]
    these_genomes  <- genome_list[these_rows]
    
    # round off the end times
    these_lineages$end_time[these_lineages$status == "active"] <- this_observation_time
    these_lineages$status[these_lineages$status == "active"]   <- "extant"
      
    # store the values
    observed_lineages[[i]] <- these_lineages
    observed_genomes[[i]]  <- these_genomes
    
  }
 
  # return
  return(list(lineages = observed_lineages, genomes = observed_genomes))
  
}










