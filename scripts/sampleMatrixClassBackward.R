sampleMatrixModelBackward <- setRefClass(
  
  "sampleMatrixModelBackward",
  
  fields = c(
    
    # state space
    "nmax",
    "states",
    "labels",
    "num_states",
    "sampled_state",
    
    # parameters
    "phi", # sample of any individual
    
    # index helpers
    "phi_indexes",
    "phi_factors",
    
    # the matrix
    "R",
    
    # flags
    "phi_dirty"
    
  ),
  
  methods = list(
    
    initialize = function(nmax_, 
                          phi_,
                          sampled_state_) {
      
      ##############
      # parameters #
      ##############
      
      setPhi(phi_)
      
      ###############
      # make states #
      ###############
      
      nmax          <<- nmax_
      states        <<- enumerateStates(nmax)
      labels        <<- rownames(states)
      num_states    <<- nrow(states)
      sampled_state <<- sampled_state_
      
      ##########################
      # create indices for phi #
      ##########################
      
      # make dummy states
      phi_states <- states
      if ( sampled_state == "0" ) {
        phi_states[,1] <- phi_states[,1] - 1  
      } else if ( sampled_state == "1" ) {
        phi_states[,2] <- phi_states[,2] - 1
      } else {
        stop("Oops, invalid sample state.")
      }
      phi_states <- phi_states[rowSums(phi_states) <= nmax & rowSums(phi_states) > 0,]
      
      # get from and to indexes
      phi_to_index   <- stateToReducedIndex(nmax, phi_states)
      phi_from_index <- match(rownames(phi_states), rownames(states))
      
      # get factors
      if ( sampled_state == "0" ) {
        phi_factors <<- states[phi_from_index,1]        
      } else if ( sampled_state == "1" ) {
        phi_factors <<- states[phi_from_index,2]
      } else {
        stop("Oops, invalid sample state.")
      }
      
      # drop anything with rate 0
      phi_from_index <-  phi_from_index[phi_factors != 0]
      phi_to_index   <-  phi_to_index[phi_factors != 0]
      phi_factors    <<- phi_factors[phi_factors != 0]
      phi_indexes    <<- cbind(from = phi_to_index, to = phi_from_index)
      
      #####################
      # create the matrix #
      #####################
      
      # make an empty matrix
      R <<- Matrix(data = 0, nrow = num_states, ncol = num_states, dimnames = list(labels, labels), sparse = TRUE)
      
      # update
      updatePhi()
      
    },
    
    show = function() {
      print(R)
    },
    
    doSampleEvent = function(p) {
      
      # just multiply
      (p %*% R)[1,]
      
    },
    
    setPhi = function(phi_) {
      phi        <<- phi_
      phi_dirty  <<- TRUE
    },
    
    updatePhi = function() {
      if ( phi_dirty == TRUE ) {
        populatePhi()
      } 
      phi_dirty <<- FALSE
    },
    
    populatePhi = function() {
      R[phi_indexes] <<- phi * phi_factors
    }
    
  )
  
)