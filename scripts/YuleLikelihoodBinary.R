library(ape)
library(expm)

YuleLikelihoodBinary <- setRefClass(
  
  "YuleLikelihoodBinary",
  
  fields = c(
    
    "tree",
    "raw_seq",
    "sel_seq",
    "sel_states",
    "initial_state",
    "singleton",
    "max_time",
    
    "lambda0",
    "lambda1",
    "gamma01",
    "gamma10",
    "phi",
    "fitnesses",
    "fitness_function",
    
    "model",
    "S",
    "edge_matrix",
    "extant_tips",
    
    "selected_sites",
    "sel_lik_container",
    "sel_lik_dirty",
    "sel_lik_container_dirty",
    "sel_rescale_log",
    "selected_site_log_likelihood",
    
    "P_sel"
    
  ),
  
  methods = list(
    
    initialize = function(tree_,
                          data_,
                          model_,
                          lambda0_,
                          lambda1_,
                          gamma01_,
                          gamma10_,
                          phi_,
                          initial_state_ = "0",
                          singleton_     = FALSE) {
      
      # store the data
      tree          <<- tree_
      raw_seq       <<- toupper(data_)
      initial_state <<- initial_state_
      singleton     <<- singleton_
      
      # store model and parameters
      model            <<- model_
      lambda0          <<- lambda0_
      lambda1          <<- lambda1_
      gamma01          <<- gamma01_
      gamma10          <<- gamma10_
      phi              <<- phi_
      
      # set everything as dirty
      sel_lik_dirty <<- TRUE
      sel_lik_container_dirty <<- TRUE
      
      # do something special if the tree is a singleton
      if ( singleton ) {
        
        # find the extant tips
        extant_tips <<- tree$tip.label

        # create the tree structure
        edge_matrix <<- data.frame(index = 2, anc = 2, desc = 1, bl = 0)

        # append the root edge
        edge_matrix <<- rbind(data.frame(index = 1, anc = 0, desc = 2, bl = tree$edge.length), edge_matrix)
        
        # create the descendant indexes
        desc_list <- vector("list", nrow(edge_matrix))
        for(i in 1:nrow(edge_matrix)) {
          # get the descendant indexes
          desc_list[[i]] <- list(edge_matrix$index[edge_matrix$anc == edge_matrix$desc[i]])
        }
        edge_matrix <<- cbind(edge_matrix, data.frame(desc_index = I(desc_list)))
        
        # store the max time
        max_time <<- 0
        
      } else {

        # find the extant tips
        extant_tips <<- tree$tip.label[tree$tip.label %in% is.extinct(tree) == FALSE]
        
        # create the tree structure
        edge_matrix <<- data.frame(index = 1 + 1:nrow(tree$edge), anc = tree$edge[,1], desc = tree$edge[,2], bl = tree$edge.length)
        
        # append the root edge
        edge_matrix <<- rbind(data.frame(index = 1, anc = 0, desc = tree$Nnode + 2, bl = tree$root.edge), edge_matrix)
        
        # create the descendant indexes
        desc_list <- vector("list", nrow(edge_matrix))
        for(i in 1:nrow(edge_matrix)) {
          # get the descendant indexes
          desc_list[[i]] <- list(edge_matrix$index[edge_matrix$anc == edge_matrix$desc[i]])
        }
        edge_matrix <<- cbind(edge_matrix, data.frame(desc_index = I(desc_list)))
        
        # store the max time
        max_time <<- max(node.depth.edgelength(tree))
        
      }
      
      # initialize the containers
      initSelectedSitesContainers()
      
      # DONE INITIALIZING
      
    },
    
    setTree = function(tree_) {
      
      # set the tree
      tree <<- tree_
      
      # update the edge matrix
      edge_matrix$bl <<- c(tree$root.edge, tree$edge.length)
      
      # set likelihood as dirty
      sel_lik_dirty <<- TRUE
      
    },
    
    initSelectedSitesContainers = function() {
      
      # make the matrix
      S <<- matrix(0, 3, 3, dimnames = list(c("0","1","A"), c("0","1","A")))
      S[1,2] <<- gamma01
      S[1,3] <<- lambda0 + phi
      S[2,1] <<- gamma10
      S[2,3] <<- lambda1 + phi
      diag(S) <<- -rowSums(S)
      
      # make fitnesses
      fitnesses <<- c(lambda0, lambda1, 0)
      
      if ( sel_lik_container_dirty ) {
        
        # update sel seq
        sel_seq <<- ifelse(raw_seq == model, "1", "0")
        
        # update the containers
        sel_lik_container <<- vector("list", nrow(edge_matrix))
        likelihood_container <- matrix(0, nrow = 3, ncol = 1)
        rownames(likelihood_container) <- c("0", "1", "A")
        
        # populate the containers
        if ( singleton  ) {
          
          # initialize the container
          this_container <- likelihood_container
          
          # get tip data
          this_sel_data <- sel_seq[tree$tip.label,]
          if ( tree$tip.label %in% extant_tips  ) {
            this_container[this_sel_data,1] <- 1.0
          } else {
            this_container[this_sel_data,1] <- phi
          }
          
          # store the likelihood
          sel_lik_container[[2]] <<- this_container
          
        } else {

          tips <- tree$tip.label
          for (i in 1:length(tips)) {
            
            # get the tip
            this_tip <- tips[i]
            this_row <- edge_matrix$index[edge_matrix$desc == i]
            
            # initialize the container
            this_container <- likelihood_container
            
            # get the data for this tip
            this_sel_data <- sel_seq[this_tip,]
            if (this_tip %in% extant_tips) {
              this_container[this_sel_data,1] <- 1.0
            } else {
              this_container[this_sel_data,1] <- phi
            }
            
            # store the likelihood
            sel_lik_container[[this_row]] <<- this_container
            
          }
          
        }
        
        # set clean
        sel_lik_container_dirty <<- FALSE
        
      }
      
    },
    
    computeSelectedTransitionProbabilities = function() {
      
      # eigen decomposition
      eigen <- eigen(S)
      vals  <- eigen$values
      vec   <- eigen$vectors
      vecI  <- solve(vec)
      tpFunc <- function(t) computeTransitionProbability(vec, vecI, exp(t * vals))
      
      # populate transition probabilities per branch
      P_sel <<- vector("list", nrow(edge_matrix))
      for(i in 1:nrow(edge_matrix)) {
        
        # get the branch length
        this_bl <- edge_matrix$bl[i]
        
        # store P
        P_sel[[i]] <<- tpFunc(this_bl)
        
      }
      
    },
    
    setModel = function(model_) {
      model <<- model_
      sel_lik_dirty <<- TRUE
      sel_lik_container_dirty <<- TRUE
    },
    
    setLambda0 = function(lambda0_) {
      lambda0 <<- lambda0_
      sel_lik_dirty <<- TRUE
    },

    setLambda1 = function(lambda1_) {
      lambda1 <<- lambda1_
      sel_lik_dirty <<- TRUE
    },
    
    setGamma01 = function(gamma01_) {
      gamma01 <<- gamma01_
      sel_lik_dirty <<- TRUE
    },

    setGamma10 = function(gamma10_) {
      gamma10 <<- gamma10_
      sel_lik_dirty <<- TRUE
    },
    
    setPhi = function(phi_) {
      phi <<- phi_
      sel_lik_dirty <<- TRUE
      sel_lik_container_dirty <<- TRUE
    },
    
    computeLikelihood = function() {
      
      # first, update the likelihoods
      computeSelectedLikelihood()
      
      # now, collect the appropriate values at the root
      ll <- getSelectedLikelihood()
      return(ll)
      
    },
    
    getSelectedLikelihood = function() {
      # only use the likelihoods for sites under selection
      return(selected_site_log_likelihood)
    },
    
    computeLikelihoodIntegrated = function(stem_ages, stem_probs) {

      # update containers
      initSelectedSitesContainers()
      
      # update transition probabilities
      computeSelectedTransitionProbabilities()

      # eigen decomposition
      eigen <- eigen(S)
      vals  <- eigen$values
      vec   <- eigen$vectors
      vecI  <- solve(vec)
      tpFunc <- function(t) computeTransitionProbability(vec, vecI, exp(t * vals))
            
      if ( singleton ) {
      
        # reset the scalar
        sel_rescale_log <<- 0
        
        # get the descendant likelihood
        desc_CL <- sel_lik_container[[2]]
        
        # for each stem age, compute the likelihood
        num_stem_ages <- length(stem_ages)
        likelihoods   <- numeric(num_stem_ages)
        for(i in 1:num_stem_ages) {
          
          # get the stem age
          this_stem_branch_length <- stem_ages[i] - max_time
          
          # compute the transition probability
          this_P <- tpFunc(this_stem_branch_length)
          
          # compute the conditional likelihood
          this_CL <- this_P %*% desc_CL
          
          # accumulate probability at the origin
          likelihoods[i] <- log(this_CL[as.numeric(initial_state) + 1]) + sel_rescale_log
          
        }
        
      } else {
        
        # compute the likelihood for the selected site
        sel_rescale_log <<- 0
        recursiveComputeConditionalLikelihoodSelected(edge_matrix$index[1])

        # get the conditional likelihoods
        desc_CL <- sel_lik_container[[1]]

        # for each stem age, compute the likelihood
        num_stem_ages <- length(stem_ages)
        likelihoods   <- numeric(num_stem_ages)
        for(i in 1:num_stem_ages) {
          
          # get the stem age
          this_stem_branch_length <- stem_ages[i] - max_time
          
          # compute the transition probability
          this_P <- tpFunc(this_stem_branch_length)
          
          # compute the conditional likelihood
          this_CL <- this_P %*% desc_CL
          
          # accumulate probability at the origin
          likelihoods[i] <- log(this_CL[as.numeric(initial_state) + 1]) + sel_rescale_log
          
        }
        
      }
      
      # include the prior
      log_posterior <- likelihoods + log(stem_probs)
      
      # rescale
      max_log_post  <- max(log_posterior)
      log_posterior <- log_posterior - max_log_post
      
      # exponentiate
      posterior <- exp(log_posterior)
      
      # compute the marginal likelihood
      if ( num_stem_ages > 1 ) {
        
        # add the sample point
        posterior <- c(0, posterior)
        stem_ages <- c(max_time, stem_ages)
        
        # integrate the likelihood
        selected_site_log_likelihood <<- log(integrateLikelihood(stem_ages, posterior)) + max_log_post  
        
      } else {
        
        # just use the one likelihood
        selected_site_log_likelihood <<- log(posterior) + max_log_post
        
      }
      
      # update the dirty flag
      sel_lik_dirty <<- TRUE
        
      return(selected_site_log_likelihood)
      
    },
    
    integrateLikelihood = function(points, values) {
      
      # get the number of nodes
      num_nodes <- length(points)
      
      # use parabola
      simp_nodes <- seq(2, num_nodes - 1, 1)
      simp_estimates <- numeric(num_nodes - 1)
      for(i in 1:length(simp_nodes)) {
        
        if (i == 1) { # left boundary
          
          # only need one parabola
          pt_idx <- 1:3 + simp_nodes[i] - 2
          xs <- points[pt_idx]
          ys <- matrix(values[pt_idx])
          
          # solve for the coefficients
          X <- matrix(c(xs^2, xs, xs^0), ncol = 3)
          coeffs <- (solve(X) %*% ys)[,1]
          
          # integrate from xs[1] to xs[2]
          simp_estimates[1] <- integrateParabola(coeffs, xs[1:2])
          
          # integrate from x[2] to xs[3]
          simp_estimates[2] <- integrateParabola(coeffs, xs[2:3])
          
        } else if (i == length(simp_nodes)) { # right boundary
          
          # get other parabola
          pt_idx <- 1:3 + simp_nodes[i] - 2
          xs <- points[pt_idx]
          ys <- matrix(values[pt_idx])
          
          # solve for the coefficients
          X <- matrix(c(xs^2, xs, xs^0), ncol = 3)
          coeffs <- (solve(X) %*% ys)[,1]
          
          # integrate from xs[1] to xs[2]
          # averaging with first estimate
          simp_estimates[i] <- (integrateParabola(coeffs, xs[1:2]) + simp_estimates[i]) / 2
          simp_estimates[i + 1] <- integrateParabola(coeffs, xs[2:3])
          
        } else { # interior, take average
          
          # get final parabola
          pt_idx <- 1:3 + simp_nodes[i] - 2
          xs <- points[pt_idx]
          ys <- matrix(values[pt_idx])
          
          # solve for the coefficients
          X <- matrix(c(xs^2, xs, xs^0), ncol = 3)
          coeffs <- (solve(X) %*% ys)[,1]
          
          # integrate from xs[1] to xs[2]
          simp_estimates[i] <- (integrateParabola(coeffs, xs[1:2]) + simp_estimates[i]) / 2
          
          # integrate from x[2] to xs[3]
          simp_estimates[i + 1] <- integrateParabola(coeffs, xs[2:3])
          
        }
        
      }
      
      return(sum(simp_estimates))
      
    },
    
    integrateParabola = function(coeffs, bounds) {
      diff(((coeffs[1] * bounds^3) / 3) + ((coeffs[2] * bounds^2) / 2) + coeffs[3] * bounds)
    },
    
    sampleAlleleAges = function(nreps, verbose = TRUE) {
      
      # first, make sure the likelihood is computed
      computeLikelihood()

      # we always assume the root begins in state 0
      root_state <- ifelse(initial_state == "0", 1, 2)

      # progress bar
      if ( verbose ) {
        bar <- txtProgressBar(style = 3, width = 40)
      }
      
      # simulate
      ages <- numeric(nreps)
      for(i in 1:nreps) {
        
        # simulate the age
        ages[i] <- recursiveSampleAlleleState(1, root_state, 0)
        
        # increment progress bar
        if ( verbose ) {
          setTxtProgressBar(bar, i / nreps)
        }
        
      }
      
      # remember we are measuring backward
      ages <- max(node.depth.edgelength(tree)) - ages
      
      # close the progress bar
      if ( verbose ) {
        close(bar)
      }
      
      # return the ages
      return(ages)
            
    },
    
    recursiveSampleAlleleState = function(index, initial_state, t0) {
      
      # get the transition probability up to this branch
      this_P <- P_sel[[index]][initial_state,]
      
      # get the conditional likelihoods at this branch
      this_CL <- sel_lik_container[[index]]
      
      # get the probability of each state
      these_state_probs <- this_P * this_CL[,1]
      
      # sample the new state
      end_state <- sample.int(3, size = 1, prob = these_state_probs)
      
      # simulate along the edge
      sim <- simulateBranchConditional(initial_state, end_state, edge_matrix$bl[index])

      # if there is a gain, return early
      if ( sim["has_gained"] == TRUE ) {
        
        # get the age of the gain
        gain_time <- t0 + sim["first_gain_time"]
        # cat("gained at time ", gain_time, "\n", sep = "")
        
        # stop
        return(gain_time)
        
      }
      
      # otherwise, simulate the descendants
      descendants <- edge_matrix$desc_index[index][[1]][[1]]
      if ( length(descendants) == 0 ) {
        # nothing to simulate
        return(NULL)
      } else {
        
        # increment time
        t0 <- t0 + edge_matrix$bl[index]
        
        # do left descendant
        left_sim  <- recursiveSampleAlleleState( descendants[1], end_state, t0)
        right_sim <- recursiveSampleAlleleState( descendants[2], end_state, t0)
        sim_age   <- c(left_sim, right_sim)
        
        # return the smaller of the two
        if ( !is.null(sim_age) ) {
          min_age <- min(sim_age)
        } else {
          min_age <- NULL
        }
        return(min_age)
        
      }
            
    },
    
    simulateBranchConditional = function(start_state, end_state, bl) {
      
      # if the start state and the end state match, 
      # first compute the probability of no events
      if ( start_state == end_state ) {
        
        # get the rates
        if ( start_state == 1 ) {
          total_rate <- lambda0 + gamma01 + phi
        } else {
          total_rate <- lambda1 + gamma10 + phi
        }
        
        # probability of no events
        no_event_prob <- dpois(0, total_rate * bl)
        
        # if there is no event, just return
        if ( runif(1) < no_event_prob ) {
          return(c(end_state = end_state, has_gained = FALSE))
        }
        
      }
        
      repeat {
        
        # simulate the branch
        branch_sim <- simulateBranch(start_state, bl)
        
        # successfully ended without a speciation or sampling event
        if ( !is.null(branch_sim) ) { 
          
          # ended in correct state
          if ( branch_sim["end_state"] == end_state ) {
            break
          }
          
        }
        
      }
      
      return(branch_sim)
      
    },
    
    simulateBranch = function(start_state, bl) {
      
      # containers
      has_gained            <- FALSE
      first_gain_time       <- numeric()
      num_sample_events     <- 0
      num_speciation_events <- 0
      
      # states
      current_time  <- 0
      current_state <- start_state
      
      # relative rates
      state_0_rates <- c(lambda0, gamma01, phi)
      state_0_rate  <- sum(state_0_rates)
      state_1_rates <- c(lambda1, gamma10, phi)
      state_1_rate  <- sum(state_1_rates)
      
      # simulate along the branch
      repeat {
        
        # get the rate
        if ( current_state == 1 ) {
          rates      <- state_0_rates
          total_rate <- state_0_rate
        } else {
          rates      <- state_1_rates
          total_rate <- state_1_rate
        }
        
        # simulate a waiting time
        current_time <- current_time + rexp(1, total_rate)
        
        # check if we finish
        if ( current_time > bl ) {
          break 
        }
        
        # otherwise, do an event
        event_probs <- rates / total_rate
        event_type  <- sample.int(3, size = 1, prob = event_probs)
        if ( event_type == 1 ) {
          # speciation event -> reject prematurely
          return(NULL)
        } else if ( event_type == 2 ) {
          # change type
          if ( current_state == 1 ) {
            current_state <- 2
            if ( has_gained == FALSE ) {
              has_gained <- TRUE
              first_gain_time <- current_time
            }
          } else {
            current_state <- 1
          }
        } else if ( event_type == 3 ) {
          # sampling event -> reject prematurely
          return(NULL)
        }
        
      }
      
      # return
      return(c(end_state = current_state, has_gained = has_gained, first_gain_time = first_gain_time))
      
    },
    
    computeSelectedLikelihood = function() {
      
      if (sel_lik_dirty) {
        
        # update containers
        initSelectedSitesContainers()
        
        # update transition probabilities
        computeSelectedTransitionProbabilities()
        if ( singleton ) {
          
          # get the transition probability
          this_P <- P_sel[[1]]
          
          # get the descendant likelihood
          desc_CL <- sel_lik_container[[2]]
          
          # compute the conditional likelihood
          this_CL <- this_P %*% desc_CL
          
          # accumulate probability at the origin
          selected_site_log_likelihood <<- log(this_CL[as.numeric(initial_state) + 1])
          
        } else {

          # compute the likelihood for the selected site
          sel_rescale_log <<- 0
          recursiveComputeConditionalLikelihoodSelected(edge_matrix$index[1])
          
          # get the transition probability
          this_P <- P_sel[[1]]
          
          # get the conditional likelihoods
          desc_CL <- sel_lik_container[[1]]
          
          # compute the conditional likelihood
          this_CL <- this_P %*% desc_CL
          
          # accumulate probability at the origin
          selected_site_log_likelihood <<- log(this_CL[as.numeric(initial_state) + 1]) + sel_rescale_log
          
        }
        
        # update the dirty flag
        sel_lik_dirty <<- FALSE
        
      }
      
    },
    
    recursiveComputeConditionalLikelihoodSelected = function(index) {
      
      # get all the descendants
      descendants <- edge_matrix$desc_index[index][[1]][[1]]
      
      if ( length(descendants) == 0 ) {
        # this is a tip
      } else {
        
        # get the branch indices
        this_index <- index
        desc_index <- descendants
        
        # this is an internal node
        for(i in 1:length(descendants)) {
          recursiveComputeConditionalLikelihoodSelected(descendants[i])
        }
        
        # compute the likelihood for the left descendants
        left_P       <- P_sel[[desc_index[1]]]
        left_CL      <- sel_lik_container[[desc_index[1]]]
        left_partial <- left_P %*% left_CL
        
        # compute the likelihood for the right descendants
        right_P       <- P_sel[[desc_index[2]]]
        right_CL      <- sel_lik_container[[desc_index[2]]]
        right_partial <- right_P %*% right_CL
        
        # compute the likelihood at the focal node
        this_CL <- left_partial * right_partial
        
        # add the speciation rate
        this_CL <- 2 * fitnesses * this_CL
        
        # rescale the likelihood
        scalar <- max(this_CL)
        this_CL <- this_CL / scalar
        sel_rescale_log <<- sel_rescale_log + log(scalar)
        
        # store the conditional likelihood
        sel_lik_container[[this_index]] <<- this_CL
        
      }
      
    }
    
  )
  
)