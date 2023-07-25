countMatrixModel <- setRefClass(
  
  "countMatrixModel",
  
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
    
    # index helpers
    # k is the flattened-state index
    # i is the number of ancestral alleles
    # j is the number of derived alleles
    "lambda0_indexes", # + 1
    "lambda0_factors", # i
    "lambda1_indexes", # + k + 1
    "lambda1_factors", # j
    "gamma01_indexes", # + k
    "gamma01_factors", # i
    "gamma10_indexes", # - k
    "gamma10_factors", # j
    
    # the matrix
    "Q",
    "J", # the jacobian of the matrix
    
    # flags
    "lambda0_dirty",
    "lambda1_dirty",
    # phi doesn't have its own flag because it only affects the diagonal
    "gamma01_dirty",
    "gamma10_dirty",
    "diagonal_dirty"
    
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
      
      ##############################
      # create indices for lambda0 #
      ##############################
      
      # make dummy states
      lambda0_states     <- states
      lambda0_states[,1] <- lambda0_states[,1] + 1
      lambda0_states     <- lambda0_states[rowSums(lambda0_states) <= nmax,]
      
      # get from and to indexes
      lambda0_to_index   <- stateToReducedIndex(nmax, lambda0_states)
      lambda0_from_index <- match(rownames(lambda0_states), rownames(states))
      
      # get factors
      lambda0_factors <<- states[lambda0_from_index,1]
      
      # drop anything with rate 0
      lambda0_from_index <-  lambda0_from_index[lambda0_factors != 0]
      lambda0_to_index   <-  lambda0_to_index[lambda0_factors != 0]
      lambda0_factors    <<- lambda0_factors[lambda0_factors != 0]
      lambda0_indexes    <<- cbind(from = lambda0_from_index, to = lambda0_to_index)
      
      ##############################
      # create indices for lambda1 #
      ##############################
      
      # make dummy states
      lambda1_states     <- states
      lambda1_states[,2] <- lambda1_states[,2] + 1
      lambda1_states     <- lambda1_states[rowSums(lambda1_states) <= nmax,]
      
      # get from and to indexes
      lambda1_to_index   <- stateToReducedIndex(nmax, lambda1_states)
      lambda1_from_index <- match(rownames(lambda1_states), rownames(states))
      
      # get factors
      lambda1_factors <<- states[lambda1_from_index,2]
      
      # drop anything with rate 0
      lambda1_from_index <-  lambda1_from_index[lambda1_factors != 0]
      lambda1_to_index   <-  lambda1_to_index[lambda1_factors != 0]
      lambda1_factors    <<- lambda1_factors[lambda1_factors != 0]
      lambda1_indexes    <<- cbind(from = lambda1_from_index, to = lambda1_to_index)
      
      ###########################################
      # create indices for mutation from 0 to 1 #
      ###########################################
      
      # make dummy states
      gamma01_states     <- states
      gamma01_states[,1] <- gamma01_states[,1] - 1
      gamma01_states[,2] <- gamma01_states[,2] + 1
      gamma01_states     <- gamma01_states[rowSums(gamma01_states) <= nmax & gamma01_states[,1] >= 0,]
      
      # get from and to indexes
      gamma01_to_index   <- stateToReducedIndex(nmax, gamma01_states)
      gamma01_from_index <- match(rownames(gamma01_states), rownames(states))
      
      # get factors
      gamma01_factors <<- states[gamma01_from_index,1]
      
      # drop anything with rate 0
      gamma01_from_index <-  gamma01_from_index[gamma01_factors != 0]
      gamma01_to_index   <-  gamma01_to_index[gamma01_factors != 0]
      gamma01_factors    <<- gamma01_factors[gamma01_factors != 0]
      gamma01_indexes    <<- cbind(from = gamma01_from_index, to = gamma01_to_index)
      
      ###########################################
      # create indices for mutation from 1 to 0 #
      ###########################################
      
      # make dummy states
      gamma10_states     <- states
      gamma10_states[,1] <- gamma10_states[,1] + 1
      gamma10_states[,2] <- gamma10_states[,2] - 1
      gamma10_states     <- gamma10_states[rowSums(gamma10_states) <= nmax & gamma10_states[,2] >= 0,]
      
      # get from and to indexes
      gamma10_to_index   <- stateToReducedIndex(nmax, gamma10_states)
      gamma10_from_index <- match(rownames(gamma10_states), rownames(states))
      
      # get factors
      gamma10_factors <<- states[gamma10_from_index,2]
      
      # drop anything with rate 0
      gamma10_from_index <-  gamma10_from_index[gamma10_factors != 0]
      gamma10_to_index   <-  gamma10_to_index[gamma10_factors != 0]
      gamma10_factors    <<- gamma10_factors[gamma10_factors != 0]
      gamma10_indexes    <<- cbind(from = gamma10_from_index, to = gamma10_to_index)
      
      #####################
      # create the matrix #
      #####################
      
      # make an empty matrix
      Q <<- Matrix(data = 0, nrow = num_states, ncol = num_states, dimnames = list(labels, labels), sparse = TRUE)
      J <<- Matrix(data = 0, nrow = num_states, ncol = num_states, dimnames = list(labels, labels), sparse = TRUE)
      
      # update the matrix
      updateAll()
      
    },
    
    solve = function(p, t, method = "ode45", atol = 1e-16, rtol = 1e-16, ...) {
      
      # make sure matrix is up to date
      updateAll()
      
      if ( method == "Higham08" ) {
        
        # matrix exponentiate
        new_p <- (p %*% expm.Higham08(Q * t, balancing = FALSE))[1,]
        # new_p <- (expm.Higham08(Q * t, balancing = FALSE) %*% t(t(p)))[,1]
        
      } else if ( method == "pracma" ) {
        
        # solve with pracma ode45
        new_p <- ode45(function(t, y, ...) {
          as.matrix(y[,1] %*% Q)
        }, t0 = 0, tfinal = t, y0 = p, rtol = rtol, atol = atol, ...)
        new_p <- new_p$y[nrow(new_p$y),]
        names(new_p) <- labels
        
      } else if ( method == "lsode" ) {
        
        new_p <- deSolve::ode(y = p, times = c(0, t), func = function(t, y, ...) {
          list((y %*% Q)[1,])
        }, rtol = rtol, atol = atol, method = "lsode", mf = 21, jactype = "fullusr", jacfunc = function(t, y, parms) { return(as.matrix(J)) }, parms = list(), ...)
        new_p <- new_p[nrow(new_p),-1]
        
      } else if ( method == "lsodes" ) {

        new_p <- deSolve::ode(y = p, times = c(0, t), func = function(t, y, ...) {
          list((y %*% Q)[1,])
        }, rtol = rtol, atol = atol, method = "lsodes", jacvec = function(t, y, j, parms) { return(J[,j]) }, parms = list(), ...)
        new_p <- new_p[nrow(new_p),-1]
                
      } else if ( method == "radau" ) {
        
        new_p <- deSolve::ode(y = p, times = c(0, t), func = function(t, y, ...) {
          list((y %*% Q)[1,])
        }, rtol = rtol, atol = atol, method = "radau", jactype = "fullusr", jacfunc = function(t, y, parms) { return(as.matrix(J)) }, parms = list(), ...)
        new_p <- new_p[nrow(new_p),-1]
        
      } else {

        new_p <- suppressWarnings(deSolve::ode(y = p, times = c(0, t), func = function(t, y, ...) {
          list((y %*% Q)[1,])
        }, rtol = rtol, atol = atol, method = method, parms = list(), ...))
        new_p <- new_p[nrow(new_p),-1]
        
      }
      
      # else if ( method == "lsoda" ) {
      #   
      #   new_p <- deSolve::lsoda(y = p, times = c(0, t), func = function(t, y, ...) {
      #     list((y %*% Q)[1,])
      #   }, parms = list(), jacfunc = function(t, y, parms) { return(J) },...)
      #   new_p <- new_p[nrow(new_p),-1]
      #   
      # } else if ( method == "lsode" ) {
      #   
      #   # solve with deSolve lsodes
      #   new_p <- deSolve::lsode(y = p, times = c(0, t), func = function(t, y, ...) {
      #     list((y %*% Q)[1,])
      #   }, parms = list(), jacfunc = function(t, y, parms) { return(J) }, ...)
      #   new_p <- new_p[nrow(new_p),-1]
      #   
      # } else if ( method == "lsodes" ) {
      #   
      #   # solve with deSolve lsodes
      #   new_p <- deSolve::lsodes(y = p, times = c(0, t), func = function(t, y, ...) {
      #     list((y %*% Q)[1,])
      #   }, parms = list(), jacvec = function(t, y, j, parms) { return(J[,j]) }, ...)
      #   new_p <- new_p[nrow(new_p),-1]
      #   
      # } else if ( method == "radau" ) {
      #   
      #   # solve with deSolve lsodes
      #   new_p <- deSolve::radau(y = p, times = c(0, t), func = function(t, y, ...) {
      #     list((y %*% Q)[1,])
      #   }, parms = list(), jacfunc = function(t, y, parms) { return(J) }, ...)
      #   new_p <- new_p[nrow(new_p),-1]
      #   
      # } else if ( method == "ode45" ) {
      #   
      #   # solve with deSolve ode45
      #   new_p <- deSolve::ode(y = p, times = c(0, t), func = function(t, y, ...) {
      #     list((y %*% Q)[1,])
      #   }, method = "ode45", parms = list(), ...)
      #   new_p <- new_p[nrow(new_p),-1]
      #   
      # } else if ( method == "ode23" ) {
      #   
      #   # solve with deSolve ode45
      #   new_p <- deSolve::ode(y = p, times = c(0, t), func = function(t, y, ...) {
      #     list((y %*% Q)[1,])
      #   }, method = "ode23", parms = list(), ...)
      #   new_p <- new_p[nrow(new_p),-1]
      #   
      # } else if ( method == "impAdams" ) {
      #   
      #   # solve with deSolve lsodes
      #   new_p <- deSolve::lsode(y = p, times = c(0, t), func = function(t, y, ...) {
      #     list((y %*% Q)[1,])
      #   }, parms = list(), jacfunc = function(t, y, parms) { return(J) }, ...)
      #   new_p <- new_p[nrow(new_p),-1]
      #   
      # } else if ( method == "rk4" ) {
      #   
      #   # solve with deSolve rk4
      #   new_p <- deSolve::ode(y = p, times = c(0, t), func = function(t, y, ...) {
      #     list((y %*% Q)[1,])
      #   }, method = "rk4", parms = list(), ...)
      #   new_p <- new_p[nrow(new_p),-1]
      #   
      # } else {
      #   stop("Choose a different method")
      # }
      
      # truncate at 0
      # new_p <- pmax(new_p, 0)
      
      return(new_p)
      
    },
    
    show = function() {
      print(Q)
    },
    
    updateAll = function() {
      updateLambda0()
      updateLambda1()
      updateGamma01()
      updateGamma10()
      updateDiagonal()
    },
    
    setLambda0 = function(lambda0_) {
      lambda0        <<- lambda0_
      lambda0_dirty  <<- TRUE
      diagonal_dirty <<- TRUE
    },
    
    updateLambda0 = function() {
      if ( lambda0_dirty == TRUE ) {
        populateLambda0()
      } 
      lambda0_dirty <<- FALSE
    },
    
    populateLambda0 = function() {
      Q[lambda0_indexes] <<- lambda0 * lambda0_factors
      J[lambda0_indexes[,2:1]] <<- lambda0 * lambda0_factors
    },
    
    setLambda1 = function(lambda1_) {
      lambda1        <<- lambda1_
      lambda1_dirty  <<- TRUE
      diagonal_dirty <<- TRUE
    },
    
    updateLambda1 = function() {
      if ( lambda1_dirty == TRUE ) {
        populateLambda1()
      } 
      lambda1_dirty <<- FALSE
    },
    
    populateLambda1 = function() {
      Q[lambda1_indexes] <<- lambda1 * lambda1_factors
      J[lambda1_indexes[,2:1]] <<- lambda1 * lambda1_factors
    },
    
    setPhi = function(phi_) {
      phi            <<- phi_
      diagonal_dirty <<- TRUE
    },
    
    setGamma01 = function(gamma01_) {
      gamma01        <<- gamma01_
      gamma01_dirty  <<- TRUE
      diagonal_dirty <<- TRUE
    },
    
    updateGamma01 = function() {
      if ( gamma01_dirty == TRUE ) {
        populateGamma01()
      } 
      gamma01_dirty <<- FALSE
    },
    
    populateGamma01 = function() {
      Q[gamma01_indexes] <<- gamma01 * gamma01_factors
      J[gamma01_indexes[,2:1]] <<- gamma01 * gamma01_factors
    },
    
    setGamma10 = function(gamma10_) {
      gamma10        <<- gamma10_
      gamma10_dirty  <<- TRUE
      diagonal_dirty <<- TRUE
    },
    
    updateGamma10 = function() {
      if ( gamma10_dirty == TRUE ) {
        populateGamma10()
      } 
      gamma10_dirty <<- FALSE
    },
    
    populateGamma10 = function() {
      Q[gamma10_indexes] <<- gamma10 * gamma10_factors
      J[gamma10_indexes[,2:1]] <<- gamma10 * gamma10_factors
    },
    
    updateDiagonal = function() {
      if ( diagonal_dirty == TRUE ) {
        populateDiagonal()
      }
      diagonal_dirty <<- FALSE
    },
    
    populateDiagonal = function() {
      leaving_rate <- phi * rowSums(states) + # sampled
        (lambda0 + gamma01) * states[,1]    + # speciation or mutation in state 0 
        (lambda1 + gamma10) * states[,2]      # speciation or mutation in state 1 
      diag(Q) <<- -leaving_rate
      diag(J) <<- -leaving_rate
    }
    
  )
  
)