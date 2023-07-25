library(cubature)
library(parallel)
source("../scripts/countClass.R", chdir = TRUE)
source("../scripts/YuleLikelihoodBinary.R", chdir = TRUE)

# settings
r0      <- 2.5
phi     <- 0.1
lambda0 <- 0.5
lambda1 <- 1.5 * lambda0
gamma   <- 0.2
nmax    <- 4

# the data
t0 <- 0   # start time
t1 <- 0.5 # time of first sample
t2 <- 1.0 # time of second sample
t3 <- 1.5 # time of third sample
t4 <- 2.0 # time of fourth sample

x0 <- "0" # state of the original lineage
x1 <- "1" # state of sample at time t1
x2 <- "0" # state of sample at time t2
x3 <- "0" # state of sample at time t3
x4 <- "1" # state of sample at time t4

################################
# likelihood under count model #
################################

count_likelihood <- function(nmax, lambda0, lambda1, phi, gamma) {
  
  # make the matrices
  Q  <- countMatrixModel(nmax, lambda0, lambda1, phi, gamma / 3, gamma)
  R1 <- sampleMatrixModel(nmax, phi, x1)
  R2 <- sampleMatrixModel(nmax, phi, x2)
  R3 <- sampleMatrixModel(nmax, phi, x3)
  
  # make init p
  start_state <- c("0" = 0, "1" = 0)
  start_state[x0] <- 1
  start_state <- paste0("(", start_state[1], ",", start_state[2], ")")
  
  p <- numeric(Q$num_states)
  names(p) <- Q$labels
  p[start_state] <- 1
  
  # integrate to time t1
  p <- Q$solve(p, t1, method = "Higham08")
  
  # do the sample event
  p <- R1$doSampleEvent(p)
  
  # integrate to time t2
  p <- Q$solve(p, t2 - t1, method = "Higham08")
  
  # do the sample event
  p <- R2$doSampleEvent(p)
  
  # integrate to time t3
  p <- Q$solve(p, t3 - t2, method = "Higham08")

  # do the sample event
  p <- R3$doSampleEvent(p)

  # integrate to time t4
  p <- Q$solve(p, t4 - t3, method = "Higham08")
  
  # get probability of present sample
  end_state <- c("0" = 0, "1" = 0)
  end_state[x4] <- 1
  end_state <- paste0("(", end_state[1], ",", end_state[2], ")")
  lik_count <- as.numeric(p[end_state])
  
  return(lik_count)
  
}

lik_count <- count_likelihood(nmax, lambda0, lambda1, phi, gamma)

###############################
# likelihood under tree model #
###############################

# data matrix
data <- matrix(c(x1, x2, x3, x4), ncol = 1)
rownames(data) <- c("x1", "x2", "x3", "x4")
data <- ifelse(data == 1, "A", "C") # convert to arbitrary nucleotides

## tree 1 ##

# ages to phylo
tree_1_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x3:T3,x4:T4):T5,x2:T2):T6,x1:T1):T7;"
  
  # compute branch lengths
  T1 <- t1 - a1
  T2 <- t2 - a2
  T3 <- t3 - a3
  T4 <- t4 - a3
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
 
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
   
}

# make initial trees
tree_1 <- tree_1_ages_to_phylo(0.2, 0.3, 0.4)

# create the calculator
tree_1_likelihood_calculator <- YuleLikelihoodBinary(tree_1, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_1_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_1_likelihood_calculator$setTree(tree_1_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_1_likelihood_calculator$computeLikelihood())
  
}

tree_1_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_1_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_1_likelihood_calculator$setLambda0(lambda0)
  tree_1_likelihood_calculator$setLambda1(lambda1)
  tree_1_likelihood_calculator$setGamma01(gamma / 3)
  tree_1_likelihood_calculator$setGamma10(gamma)
  tree_1_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_1_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t2, zmin = function(x, y) y, zmax = t3, reltol = tol)
  
  return(int)
  
}

tree_1_likelihood(lambda0, lambda1, phi, gamma)

# plot(tree_1, root.edge = TRUE, direction = "up")

## tree 2 ##

# ages to phylo
tree_2_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x3:T3,x4:T4):T5,x1:T1):T6,x2:T2):T7;"
  
  # compute branch lengths
  T1 <- t1 - a2
  T2 <- t2 - a1
  T3 <- t3 - a3
  T4 <- t4 - a3
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_2 <- tree_2_ages_to_phylo(0.2, 0.3, 0.4)

# create the calculator
tree_2_likelihood_calculator <- YuleLikelihoodBinary(tree_2, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_2_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_2_likelihood_calculator$setTree(tree_2_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_2_likelihood_calculator$computeLikelihood())
  
}

tree_2_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_2_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_2_likelihood_calculator$setLambda0(lambda0)
  tree_2_likelihood_calculator$setLambda1(lambda1)
  tree_2_likelihood_calculator$setGamma01(gamma / 3)
  tree_2_likelihood_calculator$setGamma10(gamma)
  tree_2_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_2_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1,
  zmin = function(x, y) {
    rep(y, length(y))
  }, zmax = t3, reltol = tol)
  
  return(int)
  
}

tree_2_likelihood(lambda0, lambda1, phi, gamma)

# plot(tree_2, root.edge = TRUE, direction = "up")

## tree 3 ##

# ages to phylo
tree_3_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "((x1:T1,x2:T2):T6,(x3:T3,x4:T4):T5):T7;"

  # compute branch lengths
  T1 <- t1 - a3
  T2 <- t2 - a3
  T3 <- t3 - a2
  T4 <- t4 - a2
  T5 <- a2 - a1
  T6 <- a3 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_3 <- tree_3_ages_to_phylo(0.2, 0.3, 0.4)

# create the calculator
tree_3_likelihood_calculator <- YuleLikelihoodBinary(tree_3, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_3_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_3_likelihood_calculator$setTree(tree_3_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_3_likelihood_calculator$computeLikelihood())
  
}

tree_3_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_3_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_3_likelihood_calculator$setLambda0(lambda0)
  tree_3_likelihood_calculator$setLambda1(lambda1)
  tree_3_likelihood_calculator$setGamma01(gamma / 3)
  tree_3_likelihood_calculator$setGamma10(gamma)
  tree_3_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_3_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, 
  xmin = 0,
  xmax = t1, 
  ymin = function(x) {
    x
  },
  ymax = t3,
  zmin = function(x, y) {
    rep(x, length(y))
  },
  zmax = t1,
  reltol = tol)
  
  return(int)
  
}

# tree_3_likelihood(lambda0, lambda1, phi, gamma)

# plot(tree_3, root.edge = TRUE, direction = "up")

## tree 4 ##

# ages to phylo
tree_4_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x1:T1,x2:T2):T5,x4:T4):T6,x3:T3):T7;"
  
  # compute branch lengths
  T1 <- t1 - a3
  T2 <- t2 - a3
  T3 <- t3 - a1
  T4 <- t4 - a2
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_4 <- tree_4_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_4, root.edge = TRUE, direction = "up")


# create the calculator
tree_4_likelihood_calculator <- YuleLikelihoodBinary(tree_4, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_4_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_4_likelihood_calculator$setTree(tree_4_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_4_likelihood_calculator$computeLikelihood())
  
}

tree_4_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_4_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_4_likelihood_calculator$setLambda0(lambda0)
  tree_4_likelihood_calculator$setLambda1(lambda1)
  tree_4_likelihood_calculator$setGamma01(gamma / 3)
  tree_4_likelihood_calculator$setGamma10(gamma)
  tree_4_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_4_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, zmin = function(x, y) rep(y, length(y)), zmax = t1, reltol = tol)
  
  return(int)
  
}

# tree_4_likelihood(lambda0, lambda1, phi, gamma)

## tree 5 ##

# ages to phylo
tree_5_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x1:T1,x2:T2):T5,x3:T3):T6,x4:T4):T7;"
  
  # compute branch lengths
  T1 <- t1 - a3
  T2 <- t2 - a3
  T3 <- t3 - a2
  T4 <- t4 - a1
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_5 <- tree_5_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_5, root.edge = TRUE, direction = "up")


# create the calculator
tree_5_likelihood_calculator <- YuleLikelihoodBinary(tree_5, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_5_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_5_likelihood_calculator$setTree(tree_5_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_5_likelihood_calculator$computeLikelihood())
  
}

tree_5_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_5_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_5_likelihood_calculator$setLambda0(lambda0)
  tree_5_likelihood_calculator$setLambda1(lambda1)
  tree_5_likelihood_calculator$setGamma01(gamma / 3)
  tree_5_likelihood_calculator$setGamma10(gamma)
  tree_5_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_5_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, zmin = function(x, y) rep(y, length(y)), zmax = t1, reltol = tol)
  
  return(int)
  
}

# tree_5_likelihood(lambda0, lambda1, phi, gamma)

## tree 6 ##

# ages to phylo
tree_6_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x2:T2,x3:T3):T5,x4:T4):T6,x1:T1):T7;"
  
  # compute branch lengths
  T1 <- t1 - a1
  T2 <- t2 - a3
  T3 <- t3 - a3
  T4 <- t4 - a2
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_6 <- tree_6_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_6, root.edge = TRUE, direction = "up")

# create the calculator
tree_6_likelihood_calculator <- YuleLikelihoodBinary(tree_6, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_6_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_6_likelihood_calculator$setTree(tree_6_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_6_likelihood_calculator$computeLikelihood())
  
}

tree_6_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_6_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_6_likelihood_calculator$setLambda0(lambda0)
  tree_6_likelihood_calculator$setLambda1(lambda1)
  tree_6_likelihood_calculator$setGamma01(gamma / 3)
  tree_6_likelihood_calculator$setGamma10(gamma)
  tree_6_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_6_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t2, zmin = function(x, y) rep(y, length(y)), zmax = t2, reltol = tol)
  
  return(int)
  
}

# tree_6_likelihood(lambda0, lambda1, phi, gamma)

## tree 7 ##

# ages to phylo
tree_7_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x2:T2,x3:T3):T5,x1:T1):T6,x4:T4):T7;"
  
  # compute branch lengths
  T1 <- t1 - a2
  T2 <- t2 - a3
  T3 <- t3 - a3
  T4 <- t4 - a1
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_7 <- tree_7_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_7, root.edge = TRUE, direction = "up")

# create the calculator
tree_7_likelihood_calculator <- YuleLikelihoodBinary(tree_7, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_7_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_7_likelihood_calculator$setTree(tree_7_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_7_likelihood_calculator$computeLikelihood())
  
}

tree_7_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_7_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_7_likelihood_calculator$setLambda0(lambda0)
  tree_7_likelihood_calculator$setLambda1(lambda1)
  tree_7_likelihood_calculator$setGamma01(gamma / 3)
  tree_7_likelihood_calculator$setGamma10(gamma)
  tree_7_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_7_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, zmin = function(x, y) y, zmax = t2, reltol = tol)
  
  return(int)
  
}

# tree_7_likelihood(lambda0, lambda1, phi, gamma)

## tree 8 ##

# ages to phylo
tree_8_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "((x1:T1,x4:T4):T6,(x2:T2,x3:T3):T5):T7;"
  
  # compute branch lengths
  T1 <- t1 - a3
  T2 <- t2 - a2
  T3 <- t3 - a2
  T4 <- t4 - a3
  T5 <- a2 - a1
  T6 <- a3 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_8 <- tree_8_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_8, root.edge = TRUE, direction = "up")

# create the calculator
tree_8_likelihood_calculator <- YuleLikelihoodBinary(tree_8, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_8_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_8_likelihood_calculator$setTree(tree_8_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_8_likelihood_calculator$computeLikelihood())
  
}

tree_8_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_8_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_8_likelihood_calculator$setLambda0(lambda0)
  tree_8_likelihood_calculator$setLambda1(lambda1)
  tree_8_likelihood_calculator$setGamma01(gamma / 3)
  tree_8_likelihood_calculator$setGamma10(gamma)
  tree_8_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_8_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, 
  xmin = 0,
  xmax = t1, 
  ymin = function(x) {
    x
  },
  ymax = t2,
  zmin = function(x, y) {
    rep(x, length(y))
  },
  zmax = t1,
  reltol = tol)
  
  return(int)
  
}

# tree_8_likelihood(lambda0, lambda1, phi, gamma)

## tree 9 ##

# ages to phylo
tree_9_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x1:T1,x4:T4):T5,x3:T3):T6,x2:T2):T7;"
  
  # compute branch lengths
  T1 <- t1 - a3
  T2 <- t2 - a1
  T3 <- t3 - a2
  T4 <- t4 - a3
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_9 <- tree_9_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_9, root.edge = TRUE, direction = "up")

# create the calculator
tree_9_likelihood_calculator <- YuleLikelihoodBinary(tree_9, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_9_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_9_likelihood_calculator$setTree(tree_9_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_9_likelihood_calculator$computeLikelihood())
  
}

tree_9_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_9_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_9_likelihood_calculator$setLambda0(lambda0)
  tree_9_likelihood_calculator$setLambda1(lambda1)
  tree_9_likelihood_calculator$setGamma01(gamma / 3)
  tree_9_likelihood_calculator$setGamma10(gamma)
  tree_9_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_9_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, zmin = function(x, y) y, zmax = t1, reltol = tol)
  
  return(int)
  
}

# tree_9_likelihood(lambda0, lambda1, phi, gamma)


## tree 10 ##

# ages to phylo
tree_10_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x1:T1,x4:T4):T5,x2:T2):T6,x3:T3):T7;"
  
  # compute branch lengths
  T1 <- t1 - a3
  T2 <- t2 - a2
  T3 <- t3 - a1
  T4 <- t4 - a3
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_10 <- tree_10_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_10, root.edge = TRUE, direction = "up")

# create the calculator
tree_10_likelihood_calculator <- YuleLikelihoodBinary(tree_10, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_10_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_10_likelihood_calculator$setTree(tree_10_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_10_likelihood_calculator$computeLikelihood())
  
}

tree_10_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_10_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_10_likelihood_calculator$setLambda0(lambda0)
  tree_10_likelihood_calculator$setLambda1(lambda1)
  tree_10_likelihood_calculator$setGamma01(gamma / 3)
  tree_10_likelihood_calculator$setGamma10(gamma)
  tree_10_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_10_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, zmin = function(x, y) y, zmax = t1, reltol = tol)
  
  return(int)
  
}

# tree_10_likelihood(lambda0, lambda1, phi, gamma)


## tree 11 ##

# ages to phylo
tree_11_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x2:T2,x4:T4):T5,x3:T3):T6,x1:T1):T7;"
  
  # compute branch lengths
  T1 <- t1 - a1
  T2 <- t2 - a3
  T3 <- t3 - a2
  T4 <- t4 - a3
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_11 <- tree_11_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_11, root.edge = TRUE, direction = "up")

# create the calculator
tree_11_likelihood_calculator <- YuleLikelihoodBinary(tree_11, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_11_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_11_likelihood_calculator$setTree(tree_11_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_11_likelihood_calculator$computeLikelihood())
  
}

tree_11_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_11_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_11_likelihood_calculator$setLambda0(lambda0)
  tree_11_likelihood_calculator$setLambda1(lambda1)
  tree_11_likelihood_calculator$setGamma01(gamma / 3)
  tree_11_likelihood_calculator$setGamma10(gamma)
  tree_11_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_11_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t2, zmin = function(x, y) y, zmax = t2, reltol = tol)
  
  return(int)
  
}

# tree_11_likelihood(lambda0, lambda1, phi, gamma)



## tree 12 ##

# ages to phylo
tree_12_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x2:T2,x4:T4):T5,x1:T1):T6,x3:T3):T7;"
  
  # compute branch lengths
  T1 <- t1 - a2
  T2 <- t2 - a3
  T3 <- t3 - a1
  T4 <- t4 - a3
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_12 <- tree_12_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_12, root.edge = TRUE, direction = "up")

# create the calculator
tree_12_likelihood_calculator <- YuleLikelihoodBinary(tree_12, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_12_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_12_likelihood_calculator$setTree(tree_12_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_12_likelihood_calculator$computeLikelihood())
  
}

tree_12_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_12_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_12_likelihood_calculator$setLambda0(lambda0)
  tree_12_likelihood_calculator$setLambda1(lambda1)
  tree_12_likelihood_calculator$setGamma01(gamma / 3)
  tree_12_likelihood_calculator$setGamma10(gamma)
  tree_12_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_12_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, zmin = function(x, y) y, zmax = t2, reltol = tol)
  
  return(int)
  
}

# tree_12_likelihood(lambda0, lambda1, phi, gamma)

## tree 13 ##
tree_13_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "((x1:T1,x3:T3):T6,(x2:T2,x4:T4):T5):T7;"
  
  # compute branch lengths
  T1 <- t1 - a3
  T2 <- t2 - a2
  T3 <- t3 - a3
  T4 <- t4 - a2
  T5 <- a2 - a1
  T6 <- a3 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_13 <- tree_13_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_13, root.edge = TRUE, direction = "up")

# create the calculator
tree_13_likelihood_calculator <- YuleLikelihoodBinary(tree_13, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_13_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_13_likelihood_calculator$setTree(tree_13_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_13_likelihood_calculator$computeLikelihood())
  
}

# tree_13_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_13_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_13_likelihood_calculator$setLambda0(lambda0)
  tree_13_likelihood_calculator$setLambda1(lambda1)
  tree_13_likelihood_calculator$setGamma01(gamma / 3)
  tree_13_likelihood_calculator$setGamma10(gamma)
  tree_13_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_13_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, 
  xmin = 0,
  xmax = t1, 
  ymin = function(x) {
    x
  },
  ymax = t2,
  zmin = function(x, y) {
    rep(x, length(y))
  },
  zmax = t1,
  reltol = tol)
  
  return(int)
  
}

# tree_13_likelihood(lambda0, lambda1, phi, gamma)

## tree 14 ##

# ages to phylo
tree_14_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x1:T1,x3:T3):T5,x4:T4):T6,x2:T2):T7;"
  
  # compute branch lengths
  T1 <- t1 - a3
  T2 <- t2 - a1
  T3 <- t3 - a3
  T4 <- t4 - a2
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_14 <- tree_14_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_14, root.edge = TRUE, direction = "up")

# create the calculator
tree_14_likelihood_calculator <- YuleLikelihoodBinary(tree_14, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_14_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_14_likelihood_calculator$setTree(tree_14_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_14_likelihood_calculator$computeLikelihood())
  
}

# tree_14_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_14_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_14_likelihood_calculator$setLambda0(lambda0)
  tree_14_likelihood_calculator$setLambda1(lambda1)
  tree_14_likelihood_calculator$setGamma01(gamma / 3)
  tree_14_likelihood_calculator$setGamma10(gamma)
  tree_14_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_14_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, zmin = function(x, y) y, zmax = t1, reltol = tol)
  
  return(int)
  
}

# tree_14_likelihood(lambda0, lambda1, phi, gamma)

## tree 15 ##

# ages to phylo
tree_15_ages_to_phylo <- function(a1, a2, a3) {
  
  # base newick string
  newick <- "(((x1:T1,x3:T3):T5,x2:T2):T6,x4:T4):T7;"
  
  # compute branch lengths
  T1 <- t1 - a3
  T2 <- t2 - a2
  T3 <- t3 - a3
  T4 <- t4 - a1
  T5 <- a3 - a2
  T6 <- a2 - a1
  T7 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  newick <- gsub("T6", T6, newick)
  newick <- gsub("T7", T7, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_15 <- tree_15_ages_to_phylo(0.2, 0.3, 0.4)
# plot(tree_15, root.edge = TRUE, direction = "up")

# create the calculator
tree_15_likelihood_calculator <- YuleLikelihoodBinary(tree_15, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# make the likelihood function as a function of time
tree_15_likelihood_time <- function(a1, a2, a3) {
  
  # set the tree
  tree_15_likelihood_calculator$setTree(tree_15_ages_to_phylo(a1, a2, a3))
  
  # return the likelihood
  exp(tree_15_likelihood_calculator$computeLikelihood())
  
}

tree_15_likelihood_time(0.2, 0.3, 0.4)

# create the integrator
tree_15_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-8) {
  
  # set parameters
  tree_15_likelihood_calculator$setLambda0(lambda0)
  tree_15_likelihood_calculator$setLambda1(lambda1)
  tree_15_likelihood_calculator$setGamma01(gamma / 3)
  tree_15_likelihood_calculator$setGamma10(gamma)
  tree_15_likelihood_calculator$setPhi(phi)
  
  int <- integral3(function(x, y, z) {
    
    # cat("*")
    r <- y
    for(i in 1:ncol(r)) {
      for(j in 1:nrow(r)) {
        z[i,j] <- tree_15_likelihood_time(x, y[i,j], z[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, zmin = function(x, y) y, zmax = t1, reltol = tol)
  
  return(int)
  
}

# tree_15_likelihood(lambda0, lambda1, phi, gamma)

# combine the functions together
tree_likelihood <- function(lambda0, lambda1, phi, gamma) {
  tree_1_likelihood(lambda0, lambda1, phi, gamma) +
  tree_2_likelihood(lambda0, lambda1, phi, gamma) +
  tree_3_likelihood(lambda0, lambda1, phi, gamma) +
  tree_4_likelihood(lambda0, lambda1, phi, gamma) +
  tree_5_likelihood(lambda0, lambda1, phi, gamma) +
  tree_6_likelihood(lambda0, lambda1, phi, gamma) +
  tree_7_likelihood(lambda0, lambda1, phi, gamma) +
  tree_8_likelihood(lambda0, lambda1, phi, gamma) +
  tree_9_likelihood(lambda0, lambda1, phi, gamma) +
  tree_10_likelihood(lambda0, lambda1, phi, gamma) +
  tree_11_likelihood(lambda0, lambda1, phi, gamma) +
  tree_12_likelihood(lambda0, lambda1, phi, gamma) +
  tree_13_likelihood(lambda0, lambda1, phi, gamma) +
  tree_14_likelihood(lambda0, lambda1, phi, gamma) +
  tree_15_likelihood(lambda0, lambda1, phi, gamma)
}


lik_tree <- tree_likelihood(lambda0, lambda1, phi, gamma)
lik_count


print(log(lik_count))
print(log(lik_tree))


print(lik_count / lik_tree)



# iterate over values of lambda1
lambda1   <- seq(0, 2, length.out = 101)
lik_count <- numeric(length(lambda1))
lik_tree  <- numeric(length(lambda1))
# bar <- txtProgressBar(style = 3, width = 40)
# for(i in 1:length(lambda1)) {
#   this_lambda1 <- lambda1[i]
#   lik_count[i] <- count_likelihood(nmax, lambda0, this_lambda1, phi, gamma)
#   lik_tree[i]  <- tree_likelihood(lambda0, this_lambda1, phi, gamma)
#   setTxtProgressBar(bar, i / length(lambda1))
# }

lik_count <- as.numeric(mclapply(lambda1, function(x) {
  count_likelihood(nmax, lambda0, x, phi, gamma)
}, mc.cores = 8))

lik_tree <- as.numeric(mclapply(lambda1, function(x) {
  tree_likelihood(lambda0, x, phi, gamma)
}, mc.cores = 8))


# plot
dir.create("figures", showWarnings = FALSE)
pdf("figures/validate_likelihood_4_samples.pdf", height = 4, width = 8)
par(mar=c(4,4,0,0)+0.1)
plot(lambda1, log(lik_count), pch = 3, xaxt = "n", yaxt = "n", xlab = bquote(lambda[1]), ylab = "log likelihood")
points(lambda1, log(lik_tree), pch = 4)
axis(1, lwd = 0, lwd.tick = 1)
axis(2, lwd = 0, lwd.tick = 1, las = 1)
legend("center", legend = c("count", "tree (integrated)"), title = "method", pch = c(3,4), bty = "n")
dev.off()

# save results
dir.create("results", showWarnings = FALSE)
res <- data.frame(lambda1 = lambda1, count = log(lik_count), tree = log(lik_tree))
write.table(res, file = "results/compare_likelihoods_4_samples.tsv", col.names = TRUE, quote = FALSE, row.names = FALSE)






