library(ape)
library(dispRity)
library(Matrix)
library(expm)
library(deSolve)
library(pracma)

# source the dependent stuff
source("utils.R")
source("countMatrixClass.R")
source("countMatrixClass2.R")
source("sampleMatrixClass.R")

# define the class
countModel <- setRefClass(
  
  "countModel",
  
  fields = c(
    
    # data
    "tree",           # the tree
    "seq",            # the alignment
    "nmax",           # the maximum number of states
    "has_stem",       # whether the process begins with 1 lineage (true) or 2 (false)

    # models
    "N",              # the nucleotide state of the derived allele
    "Q",              # the rate matrix
    "R0",             # the sample matrix for ancestral allele
    "R1",             # the sample matrix for derived allele
    
    # counts over time
    "temporal_data",  # all the samples and times
    "sample_times",   # the times of samples
    "sample_states",  # the states at each sample time
    "final_time",     # the time at the end of the process
    "final_state",    # the state at the end of the process as a duple
    
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
                          gamma10_) {
      
      # store data
      tree     <<- tree_
      seq      <<- toupper(seq_)
      nmax     <<- length(tree$tip.label)
      has_stem <<- is.null(tree$root.edge) == FALSE

      # create the model
      # Q  <<- countMatrixModel(nmax, lambda0_, lambda1_, phi_, gamma01_, gamma10_)
      Q  <<- countMatrixModel2(nmax, lambda0_, lambda1_, phi_, gamma01_, gamma10_)
      R0 <<- sampleMatrixModel(nmax, phi_, "0")
      R1 <<- sampleMatrixModel(nmax, phi_, "1")
      
      # create the count data over time
      temporal_data <<- getSampleData(tree, seq)
      
      final_time <<- max(temporal_data$ages)
      setDerivedNucleotide(N_)
      
      # create probabilities
      likelihood_dirty <<- TRUE
      
    },
    
    computeLikelihood = function(method = "ode45", verbose = interactive(), ...) {
      
      # only update if dirty
      if (likelihood_dirty == TRUE) {
        updateLikelihood(method, verbose, ...)
      }
      
      return(likelihood)
      
    },
    
    updateLikelihood = function(method = "ode45", verbose = interactive(), ...) {
      
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
        # cat(next_time, "\t", sum(p), "\n", sep = "")
        
        # increment time
        current_time <- next_time

        if ( verbose ) {
          setTxtProgressBar(bar, current_time / final_time)
        }
        
      }
      
      # integrate to the present
      p <- Q$solve(p, final_time - current_time, method = method, ...)
      
      if ( verbose ) {
        setTxtProgressBar(bar, 1)
        close(bar)
      }
      
      # get log-likelihood, re-add scalar
      likelihood <<- as.numeric(log(p[final_state])) + scalar
      
      # include root event if necessary
      if ( has_stem == FALSE ) {
        likelihood <<- likelihood + log(lambda0)
      }
      
      likelihood_dirty <<- FALSE
      
    },
    
    setDerivedNucleotide = function(N_) {
      
      # set the sampled nucleotide
      N <<- toupper(N_)
      
      # update temporal data
      temporal_data$allele <<- ifelse(temporal_data$state == N, 1, 0)
      
      # update sample times and states
      sample_times  <<- temporal_data$age[temporal_data$type == "sample"]
      sample_states <<- temporal_data$allele[temporal_data$type == "sample"]
      
      # if there are no samples, let's throw in some dummy time points
      # to stop the integrator, anyway
      if ( length(sample_times) == 0 ) {
        sample_times  <<- seq(0, final_time, length.out = 101)[2:100]
        sample_states <<- rep(3, 99) # 3 is a flag for a nonsense state
      }
      
      # update final state
      final_states <- temporal_data$allele[temporal_data$type == "extant"]
      final_state <<- paste0("(", sum(final_states == 0), ",", sum(final_states == 1), ")")
      
      # update p_init
      p_init <<- numeric(Q$num_states)
      names(p_init) <<- Q$labels
      if ( has_stem == TRUE ) {
        p_init["(1,0)"] <<- 1.0
      } else {
        p_init["(2,0)"] <<- 1.0
      }
      
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

















