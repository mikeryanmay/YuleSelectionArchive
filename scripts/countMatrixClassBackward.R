library(Rcpp)
library(RcppEigen)

sourceCpp("computeDerivative.cpp")

countMatrixModelBackward <- setRefClass(
  
  "countMatrixModelBackward",
  
  fields = c(
    
    # state space
    "nmax",
    "states",
    "labels",
    "num_states",
    
    # parameters
    "lambda0", # birth of ancestral individual
    "lambda1", # birth of derived individual
    "phi",     # sample of any individual
    "gamma01", # mutation from ancestral to derived
    "gamma10", # mutation from derived to ancestral
    
    # the index table
    "index_table",
    "derivative_holder"
    
  ),
  
  methods = list(
    
    initialize = function(nmax_, 
                          lambda0_, 
                          lambda1_, 
                          phi_, 
                          gamma01_, 
                          gamma10_) {
      
      ##############
      # parameters #
      ##############
      
      setLambda0(lambda0_)
      setLambda1(lambda1_)
      setPhi(phi_)
      setGamma01(gamma01_)
      setGamma10(gamma10_)
      
      ###############
      # make states #
      ###############
      
      nmax       <<- nmax_
      states     <<- enumerateStates(nmax)
      labels     <<- rownames(states)
      num_states <<- nrow(states)
      
      ###########################
      # create the factor table #
      ###########################
      
      index_table <<- data.frame(states,
                                 lambda0_index = 0, lambda0_factor = 0,
                                 lambda1_index = 0, lambda1_factor = 0,
                                 gamma01_index = 0, gamma01_factor = 0,
                                 gamma10_index = 0, gamma10_factor = 0)

      ##############################
      # create indices for lambda0 #
      ##############################

      # determine connected states
      lambda0_states     <- states
      lambda0_states[,1] <- lambda0_states[,1] + 1
      lambda0_index      <- stateToReducedIndex(nmax, lambda0_states)
      lambda0_factor     <- lambda0_states[,1] - 1
      
      # truncate
      to_remove      <- rowSums(lambda0_states) > nmax | lambda0_factor <= 0
      lambda0_states <- lambda0_states[to_remove == FALSE,]
      lambda0_index  <- lambda0_index[to_remove == FALSE]
      lambda0_factor <- lambda0_factor[to_remove == FALSE]
      
      # append
      index_table[rownames(lambda0_states), "lambda0_index"]  <<- lambda0_index - 1 # zero indexing on c-side
      index_table[rownames(lambda0_states), "lambda0_factor"] <<- lambda0_factor
      
      ##############################
      # create indices for lambda1 #
      ##############################

      # determine connected states
      lambda1_states     <- states
      lambda1_states[,2] <- lambda1_states[,2] + 1
      lambda1_index      <- stateToReducedIndex(nmax, lambda1_states)
      lambda1_factor     <- lambda1_states[,2] - 1
      
      # truncate
      to_remove      <- rowSums(lambda1_states) > nmax | lambda1_factor <= 0
      lambda1_states <- lambda1_states[to_remove == FALSE,]
      lambda1_index  <- lambda1_index[to_remove == FALSE]
      lambda1_factor <- lambda1_factor[to_remove == FALSE]
      
      # append
      index_table[rownames(lambda1_states), "lambda1_index"]  <<- lambda1_index - 1 # zero indexing on c-side
      index_table[rownames(lambda1_states), "lambda1_factor"] <<- lambda1_factor
      
      ###########################################
      # create indices for mutation from 0 to 1 #
      ###########################################
      
      # determine connected states
      gamma01_states     <- states
      gamma01_states[,1] <- gamma01_states[,1] - 1
      gamma01_states[,2] <- gamma01_states[,2] + 1
      
      # drop badness
      gamma01_states <- gamma01_states[gamma01_states[,1] >= 0,]
      
      # get indexes and factors
      gamma01_index      <- stateToReducedIndex(nmax, gamma01_states)
      gamma01_factor     <- gamma01_states[,1] + 1
      
      # truncate
      to_remove      <- rowSums(gamma01_states) > nmax | gamma01_factor <= 0 | gamma01_factor > nmax
      gamma01_states <- gamma01_states[to_remove == FALSE,]
      gamma01_index  <- gamma01_index[to_remove == FALSE]
      gamma01_factor <- gamma01_factor[to_remove == FALSE]
      
      # append
      index_table[rownames(gamma01_states), "gamma01_index"]  <<- gamma01_index - 1 # zero indexing on c-side
      index_table[rownames(gamma01_states), "gamma01_factor"] <<- gamma01_factor
      
      ###########################################
      # create indices for mutation from 1 to 0 #
      ###########################################

      # determine connected states
      gamma10_states     <- states
      gamma10_states[,1] <- gamma10_states[,1] + 1
      gamma10_states[,2] <- gamma10_states[,2] - 1
      
      # drop badness
      gamma10_states <- gamma10_states[gamma10_states[,2] >= 0,]
      
      # get indexes and factors
      gamma10_index      <- stateToReducedIndex(nmax, gamma10_states)
      gamma10_factor     <- gamma10_states[,2] + 1
      
      # truncate
      to_remove      <- rowSums(gamma10_states) > nmax | gamma10_factor <= 0 | gamma10_factor > nmax | gamma10_index < 1
      gamma10_states <- gamma10_states[to_remove == FALSE,]
      gamma10_index  <- gamma10_index[to_remove == FALSE]
      gamma10_factor <- gamma10_factor[to_remove == FALSE]
      
      # append
      index_table[rownames(gamma10_states), "gamma10_index"]  <<- gamma10_index - 1 # zero indexing on c-side
      index_table[rownames(gamma10_states), "gamma10_factor"] <<- gamma10_factor

      #############
      # finish up #
      #############
      
      # reformat table
      index_table <<- matrix(as.integer(as.matrix(index_table)), nrow = num_states)
      derivative_holder <<- numeric(num_states)
      
    },
    
    computeDerivative = function(p) {
      computeDerativeCpp(p, index_table, derivative_holder, lambda0, lambda1, gamma01, gamma10, phi)
      # derivative_holder <<- -1.0 * derivative_holder
    },
    
    solve = function(p, t, method = "ode45", atol = 1e-16, rtol = 1e-16, ...) {
      
      if ( method == "Higham08" ) {

        stop("Higham08 is not available for this class.")        

      } else if ( method == "pracma" ) {
        
        # solve with pracma ode45
        new_p <- ode45(function(t, y, ...) {
          as.matrix(computeDerivative(p))
        }, t0 = 0, tfinal = t, y0 = p, rtol = rtol, atol = atol, ...)
        new_p <- new_p$y[nrow(new_p$y),]
        names(new_p) <- labels
        
      } else {

        new_p <- deSolve::ode(y = p, times = c(0, t), func = function(t, y, ...) {
          computeDerivative(y)
          list(derivative_holder)
        }, rtol = rtol, atol = atol, method = method, parms = list(), ...)
        new_p <- new_p[nrow(new_p),-1]

      }
      
      # truncate at 0
      new_p <- pmax(new_p, 0)
      
      return(new_p)
      
    },
    
    setLambda0 = function(lambda0_) {
      lambda0        <<- lambda0_
    },
    
    setLambda1 = function(lambda1_) {
      lambda1        <<- lambda1_
    },

    setPhi = function(phi_) {
      phi            <<- phi_
    },
    
    setGamma01 = function(gamma01_) {
      gamma01        <<- gamma01_
    },

    setGamma10 = function(gamma10_) {
      gamma10        <<- gamma10_
    }
    
  )
  
)