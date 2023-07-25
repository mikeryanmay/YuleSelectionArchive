library(ape)
library(phangorn)
library(parallel)

simulatePolymorphicNeutralSites <- function(tree, gamma, num_sites) {
  
  # upscale the number of sites to simulate
  empirical_estimated_num_sites <- num_sites
  
  # simulate the sites
  sites <- as.DNAbin(simSeq(tree, l = empirical_estimated_num_sites, rate = gamma / 3))
  
  repeat {
    
    # determine which sites are variable
    sites <- sites[,apply(sites, 2, function(x) length(unique(x))) > 1]
    
    # check if there are enough
    if ( ncol(sites) < num_sites ) {
      
      # simulate some new sites
      new_sites <- as.DNAbin(simSeq(tree, l = empirical_estimated_num_sites, rate = gamma / 3))
      
      # add to existing sites
      sites <- cbind(sites, new_sites)
      
    } else {
      sites <- sites[,1:num_sites]
      break
    }
    
  }
  
  return(sites)
  
}

simulateYuleSelection <- function(S, lambda0, delta, phi, gamma, n, min_num_gains = 1, model = "additive", verbose = FALSE) {
  
  num_sims  <- 0
  num_gains <- 0
  while (any(num_gains < min_num_gains)) {
    
    sim <- NULL
    while ( is.null(sim) ) {
      sim <- simulateYuleSelectionTree(S, lambda0, delta, phi, gamma, n, model)  
    }
    
    num_gains <- sim$num_gains
    num_sims  <- num_sims + 1
    
    if ( verbose == TRUE ) {
      cat(num_sims, "\n")
    }
    
  }
  
  # append the number of simulations
  sim$num_sims <- num_sims
  
  # done
  return(sim)
  
}

simulateYuleSelection2 <- function(S, lambda0, delta, phi, gamma, n, model = "additive", verbose = FALSE) {
  
  num_sims <- 0
  num_extinct <- 0
  num_no_gains <- 0
  num_reversions <- 0
  num_multiple_sites <- 0
  num_no_extant_mutants <- 0
  not_enough_samples <- 0
  repeat {
    
    # increment the number of simulations
    num_sims  <- num_sims + 1
    
    # simulate
    sim <- simulateYuleSelectionTree(S, lambda0, delta, phi, gamma, n, model)  
    if ( is.null(sim) ) {
      num_extinct <- num_extinct + 1
      next
    }
    
    # get the number of gains
    num_gains <- sim$num_gains
    
    # check for no hits
    if ( sum(num_gains) == 0 ) {
      num_no_gains <- num_no_gains + 1 
    }
    
    # check for multiple hits at one site
    if ( any(num_gains == 2) ) {
      num_reversions <- num_reversions + 1
    }
    
    # check for multiple hits at different sites
    if ( sum(num_gains > 0) > 1 ) {
      num_multiple_sites <- num_multiple_sites + 1
    }
    
    # # only retain if there is extant variation
    # has_variation <- FALSE
    # if ( sum(num_gains) == 1 ) {
    #   if ( S[num_gains == 1] %in% sim$seq[sim$extant, num_gains == 1] ) {
    #     has_variation <- TRUE
    #   } else {
    #     num_no_extant_mutants <- num_no_extant_mutants + 1
    #   }
    # }
    
    # only retain 
    enough_samples <- FALSE
    if ( sum(num_gains) == 1 ) {
      if ( mean(sim$seq[,num_gains == 1] == S[num_gains == 1]) > 0.1 & mean(sim$seq[,num_gains == 1] == S[num_gains == 1]) < 0.9 ) {
        enough_samples <- TRUE
      } else {
        not_enough_samples <- not_enough_samples + 1
      }
    }
    
    if ( verbose == TRUE ) {
      cat(num_sims, "\n")
    }
    
    # if ( sum(num_gains) == 1 & has_variation ) {
    if ( sum(num_gains) == 1 & enough_samples ) {
      break
    }
    
  }
  
  # append the info
  sim$num_sims              <- num_sims
  sim$num_extinct           <- num_extinct
  sim$num_no_gains          <- num_no_gains
  sim$num_reversions        <- num_reversions
  sim$num_multiple_sites    <- num_multiple_sites
  sim$num_no_extant_mutants <- num_no_extant_mutants
  sim$not_enough_samples    <- not_enough_samples
  sim$site                  <- which(sim$num_gains == 1)
  
  # done
  return(sim)
  
}

simulateYuleSelectionTime <- function(S, lambda0, delta, phi, gamma, t, max_size = 1000, model = "additive", verbose = FALSE) {
  
  num_sims <- 0
  num_extinct <- 0
  num_no_gains <- 0
  num_reversions <- 0
  num_multiple_sites <- 0
  num_too_large <- 0
  num_no_extant_mutants <- 0
  repeat {
  
    # increment the number of simulations
    num_sims  <- num_sims + 1
    
    # simulate
    sim <- simulateYuleSelectionTreeTime(S, lambda0, delta, phi, gamma, t, max_size, model)  
    if ( is.null(sim) ) {
      num_extinct <- num_extinct + 1
      next
    }
    
    # get the number of gains
    num_gains <- sim$num_gains
    
    # check for no hits
    if ( sum(num_gains) == 0 ) {
      num_no_gains <- num_no_gains + 1 
    }
    
    # check for multiple hits at one site
    if ( any(num_gains == 2) ) {
      num_reversions <- num_reversions + 1
    }
    
    # check for multiple hits at different sites
    if ( sum(num_gains > 0) > 1 ) {
      num_multiple_sites <- num_multiple_sites + 1
    }
    
    # reject if too big
    if ( length(sim$phy$tip.label) > max_size ) {
      num_too_large <- num_too_large + 1
    }
    
    # only retain if there is extant variation
    has_variation <- FALSE
    if ( sum(num_gains) == 1 ) {
      if ( S[num_gains == 1] %in% sim$seq[sim$extant, num_gains == 1] ) {
        has_variation <- TRUE
      } else {
        num_no_extant_mutants <- num_no_extant_mutants + 1
      }
    }
    
    if ( verbose == TRUE ) {
      cat(num_sims, "\n")
    }

    if ( sum(num_gains) == 1 & length(sim$phy$tip.label) <= max_size & has_variation ) {
      break
    }
    
  }
  
  # append the info
  sim$num_sims              <- num_sims
  sim$num_extinct           <- num_extinct
  sim$num_no_gains          <- num_no_gains
  sim$num_reversions        <- num_reversions
  sim$num_multiple_sites    <- num_multiple_sites
  sim$num_too_large         <- num_too_large
  sim$num_no_extant_mutants <- num_no_extant_mutants
  sim$site                  <- which(sim$num_gains == 1)
  
  # done
  return(sim)
  
}

getFitnessForSequence <- function(x, S, lambda, delta, model) {
  
  # get the number of matches
  num_matches <- nchar(x) - adist(x, S, costs = c(Inf, Inf, 1))[1,1]
  
  if ( model == "additive" ) {
    fitness <- lambda + delta * num_matches
  } else if ( model == "multiplicative" ) {
    fitness <- lambda0 * delta^num_matches
  } else {
    stop("Invalid fitness model")
  }
  
  return(fitness)
  
}

additiveFitnessFunction <- function(S, lambda0, delta) {
  
  # get the number of sites
  num_sites <- length(S)
  
  # enumerate all state combinations
  nucleotides <- c("A","C","G","T")
  state_combos <- expand.grid(lapply(1:num_sites, function(x) nucleotides), stringsAsFactors = FALSE)
  states <- apply(state_combos, 1, paste0, collapse = "")
  
  # compute the fitness for each state
  fitness <- apply(state_combos, 1, function(x) {
    lambda0 + delta * sum(x == S)
  })
  names(fitness) <- states
  
  return(fitness)
  
}

multiplicativeFitnessFunction <- function(S, lambda0, delta) {

  # get the number of sites
  num_sites <- length(S)
  
  # enumerate all state combinations
  nucleotides <- c("A","C","G","T")
  state_combos <- expand.grid(lapply(1:num_sites, function(x) nucleotides), stringsAsFactors = FALSE)
  states <- apply(state_combos, 1, paste0, collapse = "")
  
  # compute the fitness for each state
  fitness <- apply(state_combos, 1, function(x) {
    lambda0 * prod(ifelse(x == S, delta, 1))
  })
  names(fitness) <- states
  
  return(fitness)
  
}

simulateYuleSelectionTree <- function(S, lambda0, delta, phi, gamma, n, model = "additive") {
  
  # get the number of sites
  num_sites <- length(S)
  nucleotides <- c("A","C","G","T")
  
  # scale the mutation rate by the number of sites
  gamma <- gamma * num_sites
  
  # choose a random starting state
  init_state <- initial_state(nucleotides, num_sites, S)
  init_state <- paste0(init_state, collapse = "")
  fit_state  <- paste0(S, collapse = "")
  
  # compute the (maximum) number of lineages (including the stem)
  num_lineages <- 2 * n - 1
  
  # pre-allocate the lineage container
  lineages   <- data.frame(anc = rep(NA, num_lineages), desc = NA, start = NA, end = NA, state = NA, fitness = NA, active = FALSE)
  
  # initialize the stem
  lineages[1,]$anc     <- 0
  lineages[1,]$start   <- 0
  lineages[1,]$state   <- init_state
  lineages[1,]$fitness <- getFitnessForSequence(init_state, fit_state, lambda0, delta, model)
  lineages[1,]$active  <- TRUE
  
  # incrementers
  num_lineages <- 1
  node_index   <- n + 1
  num_samples  <- 0
  num_gains    <- numeric(num_sites)
  gain_times   <- numeric(num_sites)
  current_time <- 0
  
  # simulate
  repeat {
    
    # determine the active lineages
    active_lineage_index <- which(lineages$active)
    
    # compute the number of active lineages
    num_active <- length(active_lineage_index)
    
    # cat(num_active, "\t", num_samples, "\n", sep = "")
    
    # stop if the number of samples is n (complete) or n is 0 (extinct)
    if ( num_active == 0 || num_active + num_samples == n ) {
      break
    } else if ( num_active + num_samples > n ) {
      stop("Too many samples... what happened?")
    }
    
    # otherwise, choose a waiting time for each lineage
    active_fitnesses <- lineages$fitness[active_lineage_index]
    active_rates     <- active_fitnesses + phi + gamma
    
    # simulate waiting times
    next_event_times <- rexp(num_active, rate = active_rates)
    current_time     <- current_time + min(next_event_times)
    
    # choose the lineage on which an event occurs
    event_lineage <- active_lineage_index[which.min(next_event_times)]
    
    # choose the type of event
    this_lineage_fitness <- lineages$fitness[event_lineage]
    relative_event_probs <- c(this_lineage_fitness, phi, gamma) / (this_lineage_fitness + phi + gamma)
    event_type           <- sample.int(3, size = 1, prob = relative_event_probs)
    
    # perform the event
    if (event_type == 1) {
      # speciation event
      
      # terminate the current lineage
      lineages[event_lineage,]$desc   <- node_index
      lineages[event_lineage,]$end    <- current_time
      lineages[event_lineage,]$active <- FALSE
      
      # add the two new lineages
      lineages[num_lineages + 1,]$anc     <- node_index
      lineages[num_lineages + 1,]$start   <- current_time
      lineages[num_lineages + 1,]$state   <- lineages[event_lineage,]$state
      lineages[num_lineages + 1,]$fitness <- lineages[event_lineage,]$fitness
      lineages[num_lineages + 1,]$active  <- TRUE
      
      lineages[num_lineages + 2,]$anc     <- node_index
      lineages[num_lineages + 2,]$start   <- current_time
      lineages[num_lineages + 2,]$state   <- lineages[event_lineage,]$state
      lineages[num_lineages + 2,]$fitness <- lineages[event_lineage,]$fitness
      lineages[num_lineages + 2,]$active  <- TRUE
      
      # increment the counters
      node_index   <- node_index + 1
      num_lineages <- num_lineages + 2
      
    } else if (event_type == 2) {
      # sampling event
      num_samples <- num_samples + 1
      lineages[event_lineage,]$desc   <- num_samples
      lineages[event_lineage,]$end    <- current_time
      lineages[event_lineage,]$active <- FALSE
    } else if (event_type == 3) {
      # mutation event

      # get the current state
      current_state <- lineages[event_lineage,]$state
      
      # choose which site mutates
      site_index <- sample.int(num_sites, size = 1)
      
      # get the state at this site
      current_site_state <- substr(current_state, site_index, site_index)
      
      # get the new state for this site
      new_site_state <- sample(nucleotides[nucleotides != current_site_state], size = 1)
      
      # create the new state
      new_state <- current_state
      substr(new_state, site_index, site_index) <- new_site_state
      
      # keep track if this is a new gain
      if ( new_site_state == S[site_index] ) {
        num_gains[site_index]  <- num_gains[site_index] + 1
        gain_times[site_index] <- current_time
      }
      
      # update the lineage
      lineages[event_lineage,]$state   <- new_state
      lineages[event_lineage,]$fitness <- getFitnessForSequence(new_state, fit_state, lambda0, delta, model)
      
    } else {
      stop("Invalid event type... what happened?")
    }
    
  }
  
  # get to the end
  num_active <- sum(lineages$active)
  if ( (num_samples + num_active) != n ) {
    # dead simulation
    return(NULL)
  } else {
    # we survived!    
  }
  
  # add some extra time
  active_fitnesses <- lineages$fitness[active_lineage_index]
  active_rates     <- active_fitnesses + phi + gamma
  next_event_times <- rexp(num_active, rate = active_rates)
  current_time     <- current_time + min(next_event_times)
  
  # finalize the lineages
  lineages[lineages$active,]$desc <- num_samples + 1:num_active
  lineages[lineages$active,]$end  <- current_time
  
  # create the phylo object
  edge        <- as.matrix(lineages[-1,1:2])
  edge.length <- (lineages$end - lineages$start)[-1]
  tip.label   <- paste0("t_", 1:n)
  Nnode       <- length(tip.label) - 1
  root.edge   <- lineages[1,]$end
  phy         <- list(edge        = edge, 
                      edge.length = edge.length,
                      tip.label   = tip.label,
                      Nnode       = Nnode,
                      root.edge   = root.edge)
  class(phy) <- "phylo"
  
  # create the sequence alignment
  seq <- matrix(lineages$state[lineages$desc <= n], ncol = 1)
  rownames(seq) <- paste0("t_", lineages$desc[lineages$desc <= n])
  seq <- t(t(seq[tip.label,]))
  
  # get sequences for nodes
  nodes <- sort(lineages$desc[lineages$desc > n])
  nodeseq <- matrix(lineages$state[match(nodes, lineages$desc)], ncol = 1)
  rownames(nodeseq) <- nodes
  
  # if we have more than one site we need to split the sequence string into a matrix
  if ( num_sites > 1 ) {
    
    # tip data
    seq_tax <- rownames(seq)
    seq <- do.call(rbind, strsplit(seq, ""))
    rownames(seq) <- seq_tax
    
    # internal data
    nodeseq_tax <- rownames(nodeseq)
    nodeseq <- do.call(rbind, strsplit(nodeseq, ""))
    rownames(nodeseq) <- nodeseq_tax
    
  }
  
  # get the extant tips
  extant_tips <- tip.label[edge[which(lineages$active) - 1,2]]
  
  # return tree and alignment
  res <- list(phy = phy, seq = seq, nodeseq = nodeseq, num_gains = num_gains, extant = extant_tips,
              gain_times = gain_times, time = current_time)
  return(res)
  
}

simulateYuleSelectionTreeTime <- function(S, lambda0, delta, phi, gamma, t, max_size = 1000, model = "additive") {
  
  # get the number of sites
  num_sites <- length(S)
  nucleotides <- c("A","C","G","T")
  
  # scale the mutation rate by the number of sites
  gamma <- gamma * num_sites
  
  # choose a random starting state
  init_state <- initial_state(nucleotides, num_sites, S)
  init_state <- paste0(init_state, collapse = "")
  fit_state  <- paste0(S, collapse = "")
  
  # compute the (maximum) number of lineages (including the stem)
  lineage_chunk <- 2 * max_size
  
  # pre-allocate the lineage container
  lineages <- data.frame(anc = rep(NA, lineage_chunk), desc = NA, start = NA, end = NA, state = NA, fitness = NA, active = FALSE)
  
  # initialize the stem
  lineages[1,]$anc     <- 0
  lineages[1,]$start   <- 0
  lineages[1,]$state   <- init_state
  lineages[1,]$fitness <- getFitnessForSequence(init_state, fit_state, lambda0, delta, model)
  lineages[1,]$active  <- TRUE
  
  # incrementers
  num_lineages <- 1
  node_index   <- 1
  num_samples  <- 0
  num_gains    <- numeric(num_sites)
  current_time <- 0
  
  # simulate
  repeat {
    
    # determine the active lineages
    active_lineage_index <- which(lineages$active)
    
    # compute the number of active lineages
    num_active <- length(active_lineage_index)
    
    # cat(num_active, "\t", num_samples, "\n", sep = "")
    
    # stop if the number of samples is 0 (extinct)
    if ( num_active == 0 ) {
      break
    }
    
    # stop if the number of samples is too big
    if ( num_active + num_samples > max_size ) {
      current_time <- t
      break
    }
    
    # otherwise, choose a waiting time for each lineage
    active_fitnesses <- lineages$fitness[active_lineage_index]
    active_rates     <- active_fitnesses + phi + gamma
    
    # simulate waiting times
    next_event_times <- rexp(num_active, rate = active_rates)
    current_time     <- current_time + min(next_event_times)
    
    # terminate if too long
    if ( current_time > t ) {
      current_time <- t
      break
    }
    
    # choose the lineage on which an event occurs
    event_lineage <- active_lineage_index[which.min(next_event_times)]
    
    # choose the type of event
    this_lineage_fitness <- lineages$fitness[event_lineage]
    relative_event_probs <- c(this_lineage_fitness, phi, gamma) / (this_lineage_fitness + phi + gamma)
    event_type           <- sample.int(3, size = 1, prob = relative_event_probs)
    
    # perform the event
    if (event_type == 1) {
      # speciation event
      
      # terminate the current lineage
      lineages[event_lineage,]$desc   <- node_index
      lineages[event_lineage,]$end    <- current_time
      lineages[event_lineage,]$active <- FALSE
      
      # add the two new lineages
      lineages[num_lineages + 1,]$anc     <- node_index
      lineages[num_lineages + 1,]$start   <- current_time
      lineages[num_lineages + 1,]$state   <- lineages[event_lineage,]$state
      lineages[num_lineages + 1,]$fitness <- lineages[event_lineage,]$fitness
      lineages[num_lineages + 1,]$active  <- TRUE
      
      lineages[num_lineages + 2,]$anc     <- node_index
      lineages[num_lineages + 2,]$start   <- current_time
      lineages[num_lineages + 2,]$state   <- lineages[event_lineage,]$state
      lineages[num_lineages + 2,]$fitness <- lineages[event_lineage,]$fitness
      lineages[num_lineages + 2,]$active  <- TRUE
      
      # increment the counters
      node_index   <- node_index + 1
      num_lineages <- num_lineages + 2
      
    } else if (event_type == 2) {
      
      # sampling event
      lineages[event_lineage,]$desc   <- node_index
      lineages[event_lineage,]$end    <- current_time
      lineages[event_lineage,]$active <- FALSE

      # increment counters
      num_samples <- num_samples + 1
      node_index  <- node_index + 1
      
    } else if (event_type == 3) {
      
      # mutation event
      
      # get the current state
      current_state <- lineages[event_lineage,]$state
      
      # choose which site mutates
      site_index <- sample.int(num_sites, size = 1)
      
      # get the state at this site
      current_site_state <- substr(current_state, site_index, site_index)
      
      # get the new state for this site
      new_site_state <- sample(nucleotides[nucleotides != current_site_state], size = 1)
      
      # create the new state
      new_state <- current_state
      substr(new_state, site_index, site_index) <- new_site_state
      
      # keep track if this is a new gain
      if ( new_site_state == S[site_index] ) {
        num_gains[site_index] <- num_gains[site_index] + 1
      }
      
      # update the lineage
      lineages[event_lineage,]$state   <- new_state
      lineages[event_lineage,]$fitness <- getFitnessForSequence(new_state, fit_state, lambda0, delta, model)
      
    } else {
      stop("Invalid event type... what happened?")
    }
    
    # check if we need to make the matrix larger
    if ( num_lineages %% lineage_chunk == 0 ) {
      new_lineages <- data.frame(anc = rep(NA, lineage_chunk), desc = NA, start = NA, end = NA, state = NA, fitness = NA, active = FALSE)
      lineages <- rbind(lineages, new_lineages)
    }
    
  }
  
  # drop NAs
  lineages <- lineages[is.na(lineages$anc) == FALSE,]
  
  # check if we made it to the end
  num_active <- sum(lineages$active)
  if ( num_active == 0 ) {
    # dead simulation
    return(NULL)
  }
  
  # check that we have a tree
  if ( num_lineages < 2 ) {
    # no tree
    return(NULL)
  }
  
  # finalize the lineages
  lineages[lineages$active,]$desc <- node_index + 1:num_active - 1
  lineages[lineages$active,]$end  <- current_time
  
  # sort by age
  lineages <- lineages[order(lineages$start, decreasing = FALSE),]
  
  # recode the indices
  edge         <- as.matrix(lineages[-1,1:2])
  rownames(edge) <- colnames(edge) <- NULL
  tip_indices  <- edge[,2] %in% edge[,1] == FALSE
  old_tip_num  <- edge[tip_indices,2]
  num_tips     <- length(old_tip_num)
  new_tip_num  <- 1:num_tips
  
  # replace tip numbers
  tmp <- edge
  tmp[match(old_tip_num, edge[,2]),2] <- new_tip_num
  
  # replace node numbers
  old_node_numbers <- unique(edge[,1])
  new_node_numbers <- num_tips + 1:length(old_node_numbers)
  node_indices     <- tip_indices == FALSE
  tmp[,1] <- new_node_numbers[match(edge[,1], old_node_numbers)]
  tmp[node_indices,2] <- new_node_numbers[match(edge[node_indices,2], old_node_numbers)]
  
  # replace edge matrix
  edge <- tmp
  
  # create the phylo object
  edge.length <- (lineages$end - lineages$start)[-1]
  tip.label   <- paste0("t_", 1:num_tips)
  Nnode       <- length(tip.label) - 1
  root.edge   <- lineages[1,]$end
  phy         <- list(edge        = edge, 
                      edge.length = edge.length,
                      tip.label   = tip.label,
                      Nnode       = Nnode,
                      root.edge   = root.edge)
  class(phy) <- "phylo"
  
  # create the sequence alignment
  seq <- matrix(lineages$state[which(tip_indices) + 1], ncol = 1)
  rownames(seq) <- tip.label
  
  # if we have more than one site we need to split the sequence string into a matrix
  if ( num_sites > 1 ) {
    
    # tip data
    seq_tax <- rownames(seq)
    seq <- do.call(rbind, strsplit(seq, ""))
    rownames(seq) <- seq_tax
    
  }

  # get the extant tips
  extant_tips <- tip.label[which(lineages[which(tip_indices),]$active)]
  
  # return tree and alignment
  res <- list(phy = phy, seq = seq, num_gains = num_gains, extant = extant_tips)
  return(res)
  
}



simulateYuleSelectionTreeSingleton <- function(S, lambda0, delta, phi, gamma, n, model = "additive") {
  
  # get the number of sites
  num_sites <- length(S)
  nucleotides <- c("A","C","G","T")
  
  # scale the mutation rate by the number of sites
  gamma_per_site <- gamma
  gamma <- gamma * num_sites
  
  if ( model == "additive" ) {
    fitness <- additiveFitnessFunction(S, lambda0, delta)
  } else if ( model == "multiplicative" ) {
    fitness <- multiplicativeFitnessFunction(S, lambda0, delta)
  } else {
    stop("Please specify either additive or multiplicative model.")
  }
  
  # choose a random starting state
  init_state <- sample(names(fitness[fitness == lambda0]), size = 1)
  
  # compute the (maximum) number of lineages (including the stem)
  num_lineages <- 2 * n - 1
  
  # pre-allocate the lineage container
  lineages   <- data.frame(anc = rep(NA, num_lineages), desc = NA, start = NA, end = NA, state = NA, fitness = NA, active = FALSE)
  
  # initialize the stem
  lineages[1,]$anc     <- 0
  lineages[1,]$start   <- 0
  lineages[1,]$state   <- init_state
  lineages[1,]$fitness <- fitness[init_state]
  lineages[1,]$active  <- TRUE
  
  # incrementers
  num_lineages <- 1
  node_index   <- n + 1
  num_samples  <- 0
  num_gains    <- numeric(num_sites)
  current_time <- 0
  
  # simulate
  repeat {
    
    # determine the active lineages
    active_lineage_index <- which(lineages$active)
    
    # compute the number of active lineages
    num_active <- length(active_lineage_index)
    
    # cat(num_active, "\t", num_samples, "\n", sep = "")
    
    # stop if the number of samples is n (complete) or n is 0 (extinct)
    if ( num_active == 0 || num_active + num_samples == n ) {
      break
    } else if ( num_active + num_samples > n ) {
      stop("Too many samples... what happened?")
    }
    
    # otherwise, choose a waiting time for each lineage
    active_fitnesses <- lineages$fitness[active_lineage_index]
    active_rates     <- active_fitnesses + phi + gamma
    
    # simulate waiting times
    next_event_times <- rexp(num_active, rate = active_rates)
    current_time     <- current_time + min(next_event_times)
    
    # choose the lineage on which an event occurs
    event_lineage <- active_lineage_index[which.min(next_event_times)]
    
    # choose the type of event
    this_lineage_fitness <- lineages$fitness[event_lineage]
    relative_event_probs <- c(this_lineage_fitness, phi, gamma) / (this_lineage_fitness + phi + gamma)
    event_type           <- sample.int(3, size = 1, prob = relative_event_probs)
    
    # perform the event
    if (event_type == 1) {
      # speciation event
      
      # terminate the current lineage
      lineages[event_lineage,]$desc   <- node_index
      lineages[event_lineage,]$end    <- current_time
      lineages[event_lineage,]$active <- FALSE
      
      # add the two new lineages
      lineages[num_lineages + 1,]$anc     <- node_index
      lineages[num_lineages + 1,]$start   <- current_time
      lineages[num_lineages + 1,]$state   <- lineages[event_lineage,]$state
      lineages[num_lineages + 1,]$fitness <- lineages[event_lineage,]$fitness
      lineages[num_lineages + 1,]$active  <- TRUE
      
      lineages[num_lineages + 2,]$anc     <- node_index
      lineages[num_lineages + 2,]$start   <- current_time
      lineages[num_lineages + 2,]$state   <- lineages[event_lineage,]$state
      lineages[num_lineages + 2,]$fitness <- lineages[event_lineage,]$fitness
      lineages[num_lineages + 2,]$active  <- TRUE
      
      # increment the counters
      node_index   <- node_index + 1
      num_lineages <- num_lineages + 2
      
    } else if (event_type == 2) {
      # sampling event
      num_samples <- num_samples + 1
      lineages[event_lineage,]$desc   <- num_samples
      lineages[event_lineage,]$end    <- current_time
      lineages[event_lineage,]$active <- FALSE
    } else if (event_type == 3) {
      
      # mutation event
      
      # get the current state
      current_state <- lineages[event_lineage,]$state
      
      # choose which site mutates
      site_index <- sample.int(num_sites, size = 1, prob = as.numeric(num_gains == 0))
      
      # get the state at this site
      current_site_state <- substr(current_state, site_index, site_index)
      
      # get the new state for this site
      new_site_state <- S[site_index]
      
      # create the new state
      new_state <- current_state
      substr(new_state, site_index, site_index) <- new_site_state
      
      # keep track if this is a new gain
      num_gains[site_index] <- 1
      
      # decrement total mutation rate
      gamma <- gamma - gamma_per_site
      
      # update the lineage
      lineages[event_lineage,]$state   <- new_state
      lineages[event_lineage,]$fitness <- fitness[new_state]
      
    } else {
      stop("Invalid event type... what happened?")
    }
    
  }
  
  # get to the end
  num_active <- sum(lineages$active)
  if ( (num_samples + num_active) != n ) {
    # dead simulation
    return(NULL)
  } else {
    # we survived!    
  }
  
  # add some extra time
  active_fitnesses <- lineages$fitness[active_lineage_index]
  active_rates     <- active_fitnesses + phi + gamma
  next_event_times <- rexp(num_active, rate = active_rates)
  current_time     <- current_time + min(next_event_times)
  
  # finalize the lineages
  lineages[lineages$active,]$desc <- num_samples + 1:num_active
  lineages[lineages$active,]$end  <- current_time
  
  # create the phylo object
  edge        <- as.matrix(lineages[-1,1:2])
  edge.length <- (lineages$end - lineages$start)[-1]
  tip.label   <- paste0("t_", 1:n)
  Nnode       <- length(tip.label) - 1
  root.edge   <- lineages[1,]$end
  phy         <- list(edge        = edge, 
                      edge.length = edge.length,
                      tip.label   = tip.label,
                      Nnode       = Nnode,
                      root.edge   = root.edge)
  class(phy) <- "phylo"
  
  # create the sequence alignment
  seq <- matrix(lineages$state[lineages$desc <= n], ncol = 1)
  rownames(seq) <- paste0("t_", lineages$desc[lineages$desc <= n])
  seq <- t(t(seq[tip.label,]))
  
  # get sequences for nodes
  nodes <- sort(lineages$desc[lineages$desc > n])
  nodeseq <- matrix(lineages$state[match(nodes, lineages$desc)], ncol = 1)
  rownames(nodeseq) <- nodes
  
  # if we have more than one site we need to split the sequence string into a matrix
  if ( num_sites > 1 ) {
    
    # tip data
    seq_tax <- rownames(seq)
    seq <- do.call(rbind, strsplit(seq, ""))
    rownames(seq) <- seq_tax
    
    # internal data
    nodeseq_tax <- rownames(nodeseq)
    nodeseq <- do.call(rbind, strsplit(nodeseq, ""))
    rownames(nodeseq) <- nodeseq_tax
    
  }
  
  # return tree and alignment
  res <- list(phy = phy, seq = seq, nodeseq = nodeseq, num_gains = num_gains)
  return(res)
  
}


# make the initial sequence
initial_state <- function(states, L, M) {
  
  # sample an initial state with minimal fitness
  sim <- numeric(L)
  for(i in 1:L) {
    sim[i] <- sample(states[states != M[i]], size = 1)
  }
  return(sim)
  
}

computeNumGains <- function(tree, optimal_sequence) {
  
  # store the data
  data <- tree$dat
  
  # find all the mutation events
  mutation_edges <- which(data$status == "mutation")
  num_mutation_edges <- length(mutation_edges)
  
  # determine if each mutation event is to the selected state
  n <- 0
  for(i in seq_len(num_mutation_edges)) {
    
    # get this edge
    this_edge <- mutation_edges[i]
    
    # get this node index
    this_idx <- data[this_edge,]$desc
    
    # find the state of the descendant edge
    this_state <- data[data$anc == this_idx,]$seq[[1]]
    
    # increment the counter if it matches the optimal state
    if (all(this_state == optimal_sequence)) {
      n <- n + 1
    }
    
  }
  
  return(n)
  
}

