library(Rcpp)
library(RcppEigen)

cppFunction("
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision

// [[Rcpp::export]]
VectorXd colMax(Map<MatrixXd> M) {
  return M.colwise().maxCoeff();
}")

cppFunction("
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision

// [[Rcpp::export]]
MatrixXd scaleMatrix(Map<MatrixXd> M, Map<VectorXd> D) {
  MatrixXd R = M;
  for(size_t i = 0; i < R.cols(); ++i) {
    R.col(i) /= D(i);
  }
  return R;
}")

cppFunction("
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision

// [[Rcpp::export]]
MatrixXd computeTransitionProbability(Map<MatrixXd> U, Map<MatrixXd> Ui, Map<VectorXd> D) {
  MatrixXd diag = D.asDiagonal();
  return U * diag * Ui;
}")

cppFunction("
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision

// [[Rcpp::export]]
MatrixXd computeTransitionProbability2(Map<MatrixXd> U, Map<MatrixXd> Ui, Map<VectorXd> D, double t) {
  // VectorXd exponential = (D * t).exp();
  VectorXd exponential = D * t;
  for(size_t i = 0; i < exponential.size(); ++i) {
    exponential(i) = std::exp(exponential(i));
  }
  MatrixXd diag = exponential.asDiagonal();
  return U * diag * Ui;
}")

drop.extinct <- function (phy, tol = NULL) {
  if (!"phylo" %in% class(phy)) {
    stop("\"phy\" is not of class \"phylo\".")
  }
  if (is.null(phy$edge.length)) {
    stop("\"phy\" does not have branch lengths.")
  }
  if (is.null(tol)) {
    tol <- min(phy$edge.length)/100
  }
  aa <- is.extinct(phy = phy, tol = tol)
  if (length(aa) > 0) {
    phy <- .drop.tip(phy, aa)
  }
  return(phy)
}

is.extinct <- function (phy, tol = NULL) {
  if (!"phylo" %in% class(phy)) {
    stop("\"phy\" is not of class \"phylo\".")
  }
  if (is.null(phy$edge.length)) {
    stop("\"phy\" does not have branch lengths.")
  }
  if (is.null(tol)) {
    tol <- min(phy$edge.length)/100
  }
  phy <- reorder(phy)
  xx <- numeric(Ntip(phy) + phy$Nnode)
  for (i in 1:length(phy$edge[, 1])) {
    xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
  }
  aa <- max(xx[1:Ntip(phy)]) - xx[1:Ntip(phy)] > tol
  if (any(aa)) {
    return(phy$tip.label[which(aa)])
  }
  else {
    return(NULL)
  }
}

.drop.tip <- function (phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0,
                       rooted = is.rooted(phy))
{
  if (missing(tip))
    return(phy)
  if (is.character(tip))
    tip <- which(phy$tip.label %in% tip)
  if (!length(tip))
    return(phy)
  phy = as.phylo(phy)
  Ntip <- length(phy$tip.label)
  tip = tip[tip %in% c(1:Ntip)]
  if (!length(tip))
    return(phy)
  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- nrow(phy$edge)
  wbl <- !is.null(phy$edge.length)
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !(edge2 %in% tip)
  ints <- edge2 > Ntip
  repeat {
    sel <- !(edge2 %in% edge1[keep]) & ints & keep
    if (!sum(sel))
      break
    keep[sel] <- FALSE
  }
  phy2 <- phy
  phy2$edge <- phy2$edge[keep, ]
  if (wbl)
    phy2$edge.length <- phy2$edge.length[keep]
  TERMS <- !(phy2$edge[, 2] %in% phy2$edge[, 1])
  oldNo.ofNewTips <- phy2$edge[TERMS, 2]
  n <- length(oldNo.ofNewTips)
  idx.old <- phy2$edge[TERMS, 2]
  phy2$edge[TERMS, 2] <- rank(phy2$edge[TERMS, 2])
  phy2$tip.label <- phy2$tip.label[-tip]
  if (!is.null(phy2$node.label))
    phy2$node.label <- phy2$node.label[sort(unique(phy2$edge[,
                                                             1])) - Ntip]
  phy2$Nnode <- nrow(phy2$edge) - n + 1L
  i <- phy2$edge > n
  phy2$edge[i] <- match(phy2$edge[i], sort(unique(phy2$edge[i]))) +
    n
  storage.mode(phy2$edge) <- "integer"
  collapse.singles(phy2)
}

makeYuleRateMatrixBinary <- function(S, fitness, phi, gamma) {
  
  # get the number of sites
  num_sites <- length(S)
  dim <- 2^num_sites + 1
  
  # enumerate all state combinations
  state_combos <- expand.grid(lapply(1:num_sites, function(x) c("0","1")), stringsAsFactors = FALSE)
  states <- c(apply(state_combos, 1, paste0, collapse = ""), "A")
  
  # create the rate matrix
  M <- matrix(0, dim, dim)
  rownames(M) <- colnames(M) <- states
  for (i in 1:(dim-1)) {
    
    # get the fitness for this state
    this_fitness <- fitness[i]
    
    # get the state
    this_state <- states[i]
    
    # get the connected states
    connected_states <- adist(this_state, states, costs = list("insertions" = Inf, "deletions" = Inf))[1,] == 1
    
    # fill in the mutation events
    M[i, connected_states] <- gamma
    
    # speciation events and sampling events are transitions
    # to an absorbing state
    M[i, dim] <- this_fitness + phi
    
  }
  
  # fill in diagonals
  diag(M) <- -rowSums(M)
  
  return(M)
  
}

additiveFitnessFunctionBinary <- function(lambda0, delta) {
  
  fitness <- c("0" = lambda, "1" = lambda + delta)
  return(fitness)
  
}

multiplicativeFitnessFunctionBinary <- function(lambda0, delta) {
  
  fitness <- c("0" = lambda, "1" = lambda * delta)
  return(fitness)
  
}

getSampleData <- function(tree, seq) {

  # get ages from tree
  init_age <- max(tree.age(tree, order = "present", digits = 6)$age)
  if ( is.null(tree$root.edge) == FALSE ) {
    age  <- init_age + tree$root.edge
    ages <- tree.age(tree, order = "present", age = age, digits = 6)
  } else {
    ages <- tree.age(tree, order = "present", age = init_age, digits = 6)
  }

  # get the time of each event
  sample_ages <- ages[grepl("t", ages$elements),]
  
  # sort the ages
  sample_ages <- sample_ages[order(sample_ages$ages),]
  
  # get the time and state for each sample
  max_age <- max(sample_ages$age)
  sample_ages$type  <- ifelse(sample_ages$ages == max_age, "extant", "sample")
  sample_ages$state <- toupper(seq[sample_ages$elements,])
  
  return(sample_ages)
  
}

getSampleDataBackward <- function(tree, seq) {
  
  # get ages from tree
  init_age <- max(tree.age(tree, order = "present", digits = 6)$age)
  if ( is.null(tree$root.edge) == FALSE ) {
    age  <- init_age + tree$root.edge
    ages <- tree.age(tree, order = "present", age = age, digits = 6)
  } else {
    ages <- tree.age(tree, order = "present", age = init_age, digits = 6)
  }
  
  # get the time of each event
  sample_ages <- ages[grepl("t", ages$elements),]
  
  # sort the ages
  sample_ages <- sample_ages[order(sample_ages$ages),]
  
  # get the time and state for each sample
  max_age <- max(sample_ages$age)
  sample_ages$type  <- ifelse(sample_ages$ages == max_age, "extant", "sample")
  sample_ages$state <- toupper(seq[sample_ages$elements,])
  
  # reverse the ages
  sample_ages$ages <- max_age - sample_ages$ages
  
  return(sample_ages)
  
}

enumerateStates <- function(kmax) {
  
  # the states for one allele
  states <- 0:kmax
  
  # combined states
  combined_states <- expand.grid("0" = states, "1" = states)
  
  # remove the first state
  combined_states <- combined_states[-1,]
  
  # remove any state with too large a sum
  combined_states <- combined_states[rowSums(combined_states) <= kmax,]
  
  return(combined_states)
  
}

makeStateLabelsFromStates <- function(states) {
  
  paste0("(",states[,1], ",", states[,2],")")
  
}

makeStateLabels <- function(kmax) {
  
  # make the states
  combined_states <- enumerateStates(kmax)
  
  # make the labels
  state_labels <- makeStateLabelsFromStates(combined_states)
  
  return(state_labels)
  
}

forwardSimulate <- function(t, init, lambda0, lambda1, gamma) {
  
  # start at time 0
  current_time <- 0
  
  # store number of individuals of each type
  nums <- init
  
  # simulate until end
  repeat {
    
    # get current rate
    current_rate <- nums[1] * (lambda0 + gamma) + nums[2] * (lambda1 + gamma)
    
    # increment time
    current_time <- current_time + rexp(1, current_rate)
    
    # if we exceed time, break
    if ( current_time > t ) {
      break
    }
    
    # determine type of event
    probs <- c(nums[1] * lambda0, nums[1] * gamma, nums[2] * lambda1, nums[2] * gamma) / current_rate
    event <- sample.int(4, size = 1, prob = probs)
    
    # execute event
    if (event == 1) {
      nums[1] <- nums[1] + 1
    } else if (event == 2) {
      nums[1] <- nums[1] - 1
      nums[2] <- nums[2] + 1
    } else if (event == 3) {
      nums[2] <- nums[2] + 1
    } else if (event == 4) {
      nums[1] <- nums[1] + 1
      nums[2] <- nums[2] - 1
    } else {
      stop("WTF.")
    }
    
  }
  
  return(nums)
  
}

enumerateStates <- function(nmax) {
  
  # the states for one allele
  states <- 0:nmax
  
  # combined states
  combined_states <- expand.grid("1" = states, "0" = states)
  combined_states <- combined_states[,2:1]
  
  # remove the first state
  combined_states <- combined_states[-1,]
  
  # remove any state with too large a sum
  combined_states <- combined_states[rowSums(combined_states) <= nmax,]
  
  # make names
  labels <- paste0("(",combined_states[,1], ",", combined_states[,2],")")
  rownames(combined_states) <- labels
  
  return(combined_states)
  
}

stateToReducedIndex <- function(nmax, state) {
  
  # the reduced index is in the state space constrained so that 0 < i + j <= n
  # given a state, what is the reduced index?
  # it's the full index of the state minus the number of states that have been removed below it,
  # which is the running sum of (i - 1)
  
  # get states i and j
  state_i <- state$`0`
  state_j <- state$`1`
  
  # get full index
  full_index <- state_j + state_i * (nmax + 1)
  
  # get reduced index
  offset        <- state_i - 1 + choose(state_i - 1, 2)
  reduced_index <- full_index - offset
  
  return(reduced_index)
  
}
