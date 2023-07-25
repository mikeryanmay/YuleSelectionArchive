library(ape)
library(dispRity)
library(Matrix)
library(expm)
library(deSolve)
library(pracma)

# source the dependent stuff
source("utils.R")
source("countMatrixClass.R")
source("countMatrixClassBackward.R")
source("sampleMatrixClass.R")
source("sampleMatrixClassBackward.R")

# define the class
countModelBackward <- setRefClass(
  
  "countModelBackward",
  
  fields = c(
    
    # data
    "tree",           # the tree
    "seq",            # the alignment
    "nmax",           # the maximum number of states
    "initial_state",  # the state at the root of the tree
    "singleton",      # whether there is a single sample
    "stem",           # whether the process begins at the stem

    # models
    "N",              # the nucleotide state of the derived allele
    "Q",              # the rate matrix
    "R0",             # the sample matrix for ancestral allele
    "R1",             # the sample matrix for derived allele
    
    # counts over time
    "temporal_data",     # all the samples and times
    "sample_times",      # the times of samples
    "sample_states",     # the states at each sample time
    "first_sample_time", # the time at the end of the process
    "final_state",       # the state at the end of the process as a duple
    
    # probabilities
    "p_init",         # the initial probability vector
    "p_final",        # the final probability vector
    "scalar",         # probability rescalar for underflow
    "likelihood",     # the probability of the data
    
    # flags
    "likelihood_dirty"
    
  ),
  
  methods = list(
    
    initialize = function(tree_,
                          seq_,
                          N_,
                          lambda0_, 
                          lambda1_, 
                          phi_, 
                          gamma01_, 
                          gamma10_,
                          initial_state_ = "0",
                          n_extra_times = 0,
                          singleton_ = FALSE,
                          stem_      = FALSE) {
      
      # store data
      tree      <<- tree_
      seq       <<- toupper(seq_)
      nmax      <<- length(tree$tip.label)
      singleton <<- singleton_
      stem      <<- stem_
      
      # make the initial state
      if ( initial_state_ == "0" ) {
        initial_state <<- "(1,0)"
      } else if ( initial_state_ == "1" ) {
        initial_state <<- "(0,1)"
      }

      # create the model
      Q  <<- countMatrixModelBackward(nmax, lambda0_, lambda1_, phi_, gamma01_, gamma10_)
      R0 <<- sampleMatrixModelBackward(nmax, phi_, "0")
      R1 <<- sampleMatrixModelBackward(nmax, phi_, "1")
      
      # create the count data over time
      if ( singleton ) {
        temporal_data     <<- data.frame(ages = 0, elements = tree$tip.label, type = "extant", state = seq[tree$tip.label,])
        first_sample_time <<- 0
      } else {
        temporal_data     <<- getSampleDataBackward(tree, seq)
        first_sample_time <<- max(temporal_data$ages)
      }
      
      # create the data
      setDerivedNucleotide(N_, n_extra_times)
      
      # create probabilities
      likelihood_dirty <<- TRUE
      
    },
    
    sampleAlleleAges = function(nreps, root_ages, root_priors, method = "ode45", verbose = interactive(), ...) {
      
      #########################################################################
      # step 1: compute the likelihood backward, storing values at each point #
      #########################################################################
      
      # get the maximum age of the process      
      max_root <- max(root_ages)
      
      # reset scalar in integrated p
      scalar <<- 0
      p <- p_init

      # make the container for the likelihoods
      num_sample_times     <- length(sample_times)
      num_root_times       <- length(root_ages)
      unscaled_likelihoods <- vector("list", length = num_sample_times + num_root_times + 1)
      likelihood_times     <- numeric(num_sample_times + num_root_times + 1)
      
      # store the first likelihood
      k <- 1
      unscaled_likelihoods[[k]] <- p
      likelihood_times[k] <- 0
      
      # make a progress bar
      if ( verbose ) {
        cat("Calculating likelihoods.\n")
        bar <- txtProgressBar(style = 3, width = 40)
      }
      
      # integrate over sample times
      current_time     <- 0
      for(i in seq_len(num_sample_times)) {
        
        # get next sample time
        next_time <- sample_times[i]
        
        # integrate
        p <- Q$solve(p, next_time - current_time, method = method, ...)
        
        # apply sample event
        if ( sample_states[i] == 0 ) {
          p <- R0$doSampleEvent(p)
        } else if ( sample_states[i] == 1 ) {
          p <- R1$doSampleEvent(p)
        } else if ( sample_states[i] == 3 ) {
          # do nothing
        } else {
          stop("Oops, couldn't find sample event.")
        }
        
        # rescale
        max_p <- max(p)
        p <- p / max_p
        scalar <<- scalar + log(max_p)
        
        # store the likelihood
        k <- k + 1
        unscaled_likelihoods[[k]] <- p
        likelihood_times[k] <- next_time
        
        # increment time
        current_time <- next_time
        if ( verbose ) {
          setTxtProgressBar(bar, current_time / max_root)
        }
        
      }
      
      # compute likelihood at each root time
      num_root_ages <- length(root_ages)
      likelihoods   <- numeric(num_root_ages)
      for(i in seq_len(num_root_ages)) {
      
        # get the next potential root age
        next_time <- root_ages[i]
        
        # integrate
        p <- Q$solve(p, next_time - current_time, method = method, ...)
        
        # compute the likelihood for this root age
        likelihoods[i] <- log(p[initial_state] * Q$lambda0) + scalar
        
        # rescale
        max_p <- max(p)
        p <- p / max_p
        scalar <<- scalar + log(max_p)
        
        # store the likelihood
        k <- k + 1
        unscaled_likelihoods[[k]] <- p
        likelihood_times[k] <- next_time
        
        # increment time
        current_time <- next_time
        if ( verbose ) {
          setTxtProgressBar(bar, current_time / max_root)
        }
        
      }
      
      # close the progress bar
      if ( verbose ) {
        setTxtProgressBar(bar, 1)
        close(bar)
      }
      
      # include the prior
      log_posterior <- likelihoods + log(root_priors)
      
      # rescale
      max_log_post  <- max(log_posterior)
      log_posterior <- log_posterior - max_log_post
      
      # exponentiate
      posterior <- exp(log_posterior)
      
      # compute the marginal likelihood
      if ( num_root_ages > 1 ) {
        
        # add the sample point
        posterior <- c(0, posterior)
        root_ages <- c(first_sample_time, root_ages)
        
        # integrate the likelihood
        marginal_likelihood <- log(integrateLikelihood(root_ages, posterior)) + max_log_post  
        
      } else {
        
        # just use the one likelihood
        marginal_likelihood <- log(posterior) + max_log_post
        
      }
      
      # compute the posterior distribution of ages
      posterior <- posterior / sum(posterior)
      cumulative_posterior <- cumsum(posterior)
      
      # approximate the cumulative posterior
      posterior_quantile_function <- approxfun(cumulative_posterior, root_ages)

      ######################################################
      # step 2: step forward, drawing states at each point #
      ######################################################
      
      # make some forward ODE objects
      Q_forward  <- countMatrixModel(nmax, Q$lambda0, Q$lambda1, Q$phi, Q$gamma01, Q$gamma10)
      R0_forward <- sampleMatrixModel(nmax, Q$phi, "0")
      R1_forward <- sampleMatrixModel(nmax, Q$phi, "1")
      
      # make an initial p vector
      p_init_forward <- numeric()
      p_init_forward <- numeric(Q$num_states)
      names(p_init_forward) <- Q$labels
      if ( initial_state == "(1,0)" ) {
        root_state < "(2,0)"
      } else if ( initial_state == "(0,1)" ) {
        root_state < "(0,2)"
      }
      p_init_forward[root_state] <- 1.0 
      
      # keep track of the gain time
      first_gain_age <- max(sample_times[sample_states == 1])
      
      # progress bar
      if ( verbose ) {
        cat("Simulating stochastic ages.\n")
        bar <- txtProgressBar(style = 3, width = 40)
      }
      
      # repeat for each replicate
      root_times <- numeric(nreps)
      gain_times <- numeric(nreps)
      for(r in 1:nreps) {
        
        # initialize the probability
        this_p <- p_init_forward
        
        # sample a start time
        this_root_age <- posterior_quantile_function(runif(1))
        root_times[r] <- this_root_age
        
        # figure out which time points to sample
        these_times <- likelihood_times[likelihood_times < this_root_age]
        nt <- length(these_times)
        
        # make a container for the state
        current_state <- root_state
        current_state_index <- which(rownames(Q$states) == current_state)
          
        # set the current time
        current_time <- this_root_age
        
        # loop over times
        for(i in nt:1) {
          
          # get the next time
          next_time <- these_times[i]
          
          # integrate forward a bit
          this_p <- Q_forward$solve(this_p, current_time - next_time, method = method, ...)
          
          # get the likelihood for this time point
          this_l <- unscaled_likelihoods[[i]]
          
          # compute the relative posterior probability of states at this time
          this_post <- this_p * this_l
          this_post <- this_post / sum(this_post)
          this_post <- pmax(this_post, 0)
          
          # sample a state
          next_state_index <- sample.int(Q$num_states, size = 1, prob = this_post)
          next_state <- names(this_post[next_state_index])
          
          # update p
          this_p[] <- 0
          this_p[next_state] <- 1.0
          
          # check if we evolved the derived allele
          if ( Q$states[next_state_index,2] > 0 ) {
            # use this age
            gain_time <- runif(1, next_time, current_time)
            break
          }
          
          # apply a sampling event, if necessary
          if (next_time %in% sample_times) {
            # get the type of event
            this_sample_type <- sample_states[sample_times == next_time]
            if ( this_sample_type == 0 ) {
              this_p <- R0$doSampleEvent(this_p)
            } else if ( this_sample_type == 1 ) {
              this_p <- R1$doSampleEvent(this_p)
            }
          }
          
          # increment time
          current_state_index <- next_state_index
          current_state       <- next_state
          current_time        <- next_time
          
        }
        
        if ( verbose ) {
          setTxtProgressBar(bar, r / nreps)
        }
        
        gain_times[r] <- gain_time
        
      }
      
      # close the progress bar
      if ( verbose ) {
        setTxtProgressBar(bar, 1)
        close(bar)
      }

      return(list(root_times = root_times, gain_times = gain_times, first_derived_sample = first_gain_age))
            
    },
    
    simulateForwardRejection = function(start_state, end_state, start_time, end_time) {
      
      repeat {
      
        cat("*")
        
        # simulate
        sim <- simulateForward(start_state, start_time, end_time)

        # make sure we have the right number of lineages, and there were no samples
        pass <- all(sim$current_state == end_state) & sim$num_samples == 0
        if ( pass ) {
          break
        }        
          
      }
      
    },
    
    simulateForward = function(start_state, start_time, end_time) {
      
      # get the parameters
      lambda0 <- Q$lambda0
      lambda1 <- Q$lambda1
      gamma01 <- Q$gamma01
      gamma10 <- Q$gamma10
      phi     <- Q$phi
      
      # set the current time
      current_time <- start_time
      
      # set the current state
      current_state <- start_state

      # counters
      num_gains   <- 0
      gain_time   <- numeric()
      num_samples <- 0
      
      # simulate forward
      repeat {
        
        # get the rate
        event_rates <- c(current_state[1] * lambda0,
                         current_state[2] * lambda1,
                         current_state[1] * gamma01,
                         current_state[2] * gamma10,
                         sum(current_state) * phi)
        
        total_rate <- sum(event_rates)
          
        # simulate time to next event
        current_time <- current_time + rexp(1, total_rate)
        
        # terminate
        if ( current_time > end_time ) {
          break
        }
        
        # otherwise, simulate an event
        event_probs <- event_rates / total_rate
        event_type  <- sample.int(5, size = 1, prob = event_probs)
        
        # do the event
        if ( event_type == 1 ) {
          current_state[1] <- current_state[1] + 1
        } else if ( event_type == 2 ) {
          current_state[2] <- current_state[2] + 1
        } else if ( event_type == 3 ) {
          current_state[1] <- current_state[1] - 1
          current_state[2] <- current_state[2] + 1
          # if this is the first gain, record the time
          if (num_gains == 0) {
            num_gains <- 1
            gain_time <- current_time
          }
        } else if ( event_type == 4 ) {
          current_state[1] <- current_state[1] + 1
          current_state[2] <- current_state[2] - 1
        } else if ( event_type == 5 ) {
          num_samples <- num_samples + 1
        }
        
      }
      
      # wrap up
      res <- list(current_state = current_state,
                  num_gains     = num_gains,
                  gain_time     = gain_time,
                  num_samples  = num_samples)
      
      return(res)
      
    },
    
    computeLikelihood = function(root_ages, root_priors, method = "ode45", verbose = interactive(), ...) {
      
      # only update if dirty
      if (likelihood_dirty == TRUE) {
        updateLikelihood(root_ages, root_priors, method, verbose, ...)
      }
      
      return(likelihood)
      
    },
    
    updateLikelihood = function(root_ages, root_priors, method = "ode45", verbose = interactive(), ...) {

      # get the maximum age of the process      
      max_root <- max(root_ages)
      
      # reset scalar in integrated p
      scalar <<- 0
      p <- p_init
      
      # make a progress bar
      if ( verbose ) {
        bar <- txtProgressBar(style = 3, width = 40)
      }
      
      # integrate over sample times
      current_time     <- 0
      num_sample_times <- length(sample_times)
      for(i in seq_len(num_sample_times)) {
        
        # get next sample time
        next_time <- sample_times[i]
        
        # integrate
        p <- Q$solve(p, next_time - current_time, method = method, ...)
        
        # apply sample event
        if ( sample_states[i] == 0 ) {
          p <- R0$doSampleEvent(p)
        } else if ( sample_states[i] == 1 ) {
          p <- R1$doSampleEvent(p)
        } else if ( sample_states[i] == 3 ) {
          # do nothing
        } else {
          stop("Oops, couldn't find sample event.")
        }
        
        # rescale
        max_p <- max(p)
        p <- p / max_p
        scalar <<- scalar + log(max_p)
        
        # cat(scalar, "\n")
        # cat(i, "\t", next_time, "\t", sum(p), "\n", sep = "")
        
        # increment time
        current_time <- next_time
        if ( verbose ) {
          setTxtProgressBar(bar, current_time / max_root)
        }
        
      }
      
      # compute likelihood at each root time
      num_root_ages <- length(root_ages)
      likelihoods   <- numeric(num_root_ages)
      for(i in seq_len(num_root_ages)) {
        
        # get the next potential root age
        next_time <- root_ages[i]
        
        # integrate
        p <- Q$solve(p, next_time - current_time, method = method, ...)
        
        # compute the likelihood for this age
        if ( singleton ) {
          likelihoods[i] <- log(p[initial_state]) + scalar
        } else {
          if ( stem ) {
            if ( initial_state == "(0,1)" ) {
              likelihoods[i] <- log(p[initial_state]) + scalar  
            } else if ( initial_state == "(1,0)" ) {
              likelihoods[i] <- log(p[initial_state]) + scalar    
            } else {
              stop("invalid initial state")
            }
          } else {
            if ( initial_state == "(0,1)" ) {
              likelihoods[i] <- log(p[initial_state] * Q$lambda1) + scalar  
            } else if ( initial_state == "(1,0)" ) {
              likelihoods[i] <- log(p[initial_state] * Q$lambda0) + scalar    
            } else {
              stop("invalid initial state")
            }
          }
        }
        
        # rescale
        max_p <- max(p)
        p <- p / max_p
        scalar <<- scalar + log(max_p)
        
        # increment time
        current_time <- next_time
        if ( verbose ) {
          setTxtProgressBar(bar, current_time / max_root)
        }
        
      }
      
      # close the progress bar
      if ( verbose ) {
        setTxtProgressBar(bar, 1)
        close(bar)
      }
      
      # include the prior
      log_posterior <- likelihoods + log(root_priors)
      
      # rescale
      max_log_post  <- max(log_posterior)
      log_posterior <- log_posterior - max_log_post
      
      # exponentiate
      posterior <- exp(log_posterior)
      
      # compute the marginal likelihood
      if ( num_root_ages > 1 ) {
        
        # add the sample point
        posterior <- c(0, posterior)
        root_ages <- c(first_sample_time, root_ages)
        
        # integrate the likelihood
        likelihood <<- log(integrateLikelihood(root_ages, posterior)) + max_log_post  
        
      } else {
        
        # just use the one likelihood
        likelihood <<- log(posterior) + max_log_post
        
      }
      
      # compute the integrated likelihood
      likelihood_dirty <<- FALSE
      
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
    
    setDerivedNucleotide = function(N_, n_extra_times) {
      
      # set the sampled nucleotide
      N <<- toupper(N_)
      
      # update temporal data
      temporal_data$allele <<- ifelse(temporal_data$state == N, 1, 0)
      
      # update sample times and states
      sample_times  <<- temporal_data$age[temporal_data$type == "sample"]
      sample_states <<- temporal_data$allele[temporal_data$type == "sample"]
      
      # add some extra times
      if ( n_extra_times > 0 ) {
        extra_sample_times  <- seq(0, first_sample_time, length.out = n_extra_times + 2)[2:(n_extra_times + 1)]
        extra_sample_states <- rep(3, length(extra_sample_times))
        sample_times        <<- c(sample_times, extra_sample_times)
        sample_states       <<- c(sample_states, extra_sample_states)
      }
      
      # reverse time
      sample_states <<- sample_states[order(sample_times, decreasing = FALSE)]
      sample_times  <<- sample_times[order(sample_times, decreasing = FALSE)]
      
      # update final state
      final_states <- temporal_data$allele[temporal_data$type == "extant"]
      final_state <<- paste0("(", sum(final_states == 0), ",", sum(final_states == 1), ")")
      
      # update p_init
      p_init <<- numeric(Q$num_states)
      names(p_init) <<- Q$labels
      p_init[final_state] <<- 1.0
      
      # make the likelihood dirty
      likelihood_dirty <<- TRUE
      
    },

    setLambda0 = function(lambda0_) {
      
      # update the integrator
      Q$setLambda0(lambda0_)
      
      # make likelihood dirty
      likelihood_dirty <<- TRUE
      
    },
    
    setLambda1 = function(lambda1_) {
      
      # update the integrator
      Q$setLambda1(lambda1_)
      
      # make likelihood dirty
      likelihood_dirty <<- TRUE
      
    },
    
    setPhi = function(phi_) {
      
      # update the integrator
      Q$setPhi(phi_)
      
      # update the sample matrices
      R0$setPhi(phi_)
      R1$setPhi(phi_)
      
      # make likelihood dirty
      likelihood_dirty <<- TRUE
      
    },
    
    setGamma01 = function(gamma01_) {
      
      # update the integrator
      Q$setGamma01(gamma01_)
      
      # make likelihood dirty
      likelihood_dirty <<- TRUE
      
    },
    
    setGamma10 = function(gamma10_) {
      
      # update the integrator
      Q$setGamma01(gamma10_)
      
      # make likelihood dirty
      likelihood_dirty <<- TRUE
      
    }
    
  )
  
)

















