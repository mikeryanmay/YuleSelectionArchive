library(ape)
library(expm)
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

makeYuleRateMatrix <- function(S, fitness, phi, gamma) {

  # get the number of sites
  num_sites <- length(S)
  dim <- 4^num_sites + 1

  # enumerate all state combinations
  state_combos <- expand.grid(lapply(1:num_sites, function(x) c("A","C","G","T")), stringsAsFactors = FALSE)
  states <- c(apply(state_combos, 1, paste0, collapse = ""), "0")

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
    M[i, connected_states] <- gamma / 3

    # speciation events and sampling events are transitions
    # to an absorbing state
    M[i, dim] <- this_fitness + phi

  }

  # fill in diagonals
  diag(M) <- -rowSums(M)

  return(M)

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

YuleLikelihood <- setRefClass(

  "YuleLikelihood",

  fields = c(

    "tree",
    "raw_seq",
    "neu_seq",
    "num_sites",
    "sel_seq",
    "sel_states",

    "lambda0",
    "gamma",
    "delta",
    "mu",
    "phi",
    "fitnesses",
    "fitness_function",

    "model",
    "S",
    "edge_matrix",
    "extant_tips",

    "neutral_sites",
    "selected_sites",
    "neu_lik_container",
    "neu_lik_dirty",
    "sel_lik_container",
    "sel_lik_dirty",
    "sel_lik_container_dirty",
    "neu_rescale_log",
    "sel_rescale_log",

    "neutral_site_log_likelihoods",
    "selected_site_log_likelihood",

    "P_neu",
    "P_sel"

  ),

  methods = list(

    initialize = function(tree_,
                          data_,
                          model_,
                          lambda0_,
                          gamma_,
                          delta_,
                          mu_,
                          phi_,
                          fitness_function_ = "additive") {

      # store the data
      tree    <<- tree_
      raw_seq <<- toupper(data_)

      # find the extant tips
      extant_tips <<- tree$tip.label[tree$tip.label %in% is.extinct(tree) == FALSE]

      # store model and parameters
      model            <<- model_
      lambda0          <<- lambda0_
      gamma            <<- gamma_
      delta            <<- delta_
      mu               <<- mu_
      phi              <<- phi_
      fitness_function <<- fitness_function_

      # set everything as dirty
      neu_lik_dirty <<- TRUE
      sel_lik_dirty <<- TRUE
      sel_lik_container_dirty <<- TRUE

      # process the data
      neutral_sites  <<- which(model == "-")
      selected_sites <<- which(model != "-")
      neu_seq        <<- raw_seq
      num_sites      <<- ncol(neu_seq)

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

      # initialize the containers
      initNeutralSitesContainers()
      initSelectedSitesContainers()

      # DONE INITIALIZING

    },

    initNeutralSitesContainers = function() {

      # create the likelihood containers
      likelihood_container <- matrix(0, nrow = 4, ncol = num_sites)
      rownames(likelihood_container) <- c("A","C","G","T")
      neu_lik_container <<- vector("list", nrow(edge_matrix))

      # populate the containers
      tips <- tree$tip.label
      for (i in 1:length(tips)) {

        # get the tip
        this_tip <- tips[i]
        this_row <- edge_matrix$index[edge_matrix$desc == i]

        # initialize the container
        this_container <- likelihood_container

        # get the neutral data for this tip
        this_neu_data <- neu_seq[this_tip,]

        # fill in conditional likelihood
        this_container[cbind(match(this_neu_data, rownames(this_container)), 1:num_sites)] <- 1.0

        # store the likelihood
        neu_lik_container[[this_row]] <<- this_container

      }

    },

    initSelectedSitesContainers = function() {

      if ( all(model == "-")  ) {

        no_selected_sites <- TRUE
        
        # make a temporary model
        tmp_model <- "A"
        
        # get the fitnesses
        if ( fitness_function == "additive" ) {
          fitnesses <<- additiveFitnessFunction(tmp_model, lambda0, delta)
        } else if (fitness_function == "multiplicative" ) {
          fitnesses <<- multiplicativeFitnessFunction(tmp_model, lambda0, delta)
        } else {
          stop("Please choose either additive or multiplicative fitness function.")
        }
        fitnesses <<- c(fitnesses, 0)
        
        # create the rate matrix
        S <<- makeYuleRateMatrix(tmp_model, fitnesses, phi, gamma)

        state_combos <- expand.grid(lapply(1:1, function(x) c("A","C","G","T")), stringsAsFactors = FALSE)
        sel_states  <<- c(apply(state_combos, 1, paste0, collapse = ""), "0")
        # fitnesses <<- apply(state_combos, 1, function(x) {
        #   lambda0 + delta * sum(x == "A")
        # })
        # fitnesses <<- c(fitnesses, 0)

      } else {

        selected_sites <<- which(model != "-")
        no_selected_sites <- FALSE
        sel_seq <<- neu_seq[,selected_sites]
        if ( inherits(sel_seq, "character") ) {
          sel_seq <<- t(t(sel_seq))
        }

        state_combos <- expand.grid(lapply(1:ncol(sel_seq), function(x) c("A","C","G","T")), stringsAsFactors = FALSE)
        sel_states  <<- c(apply(state_combos, 1, paste0, collapse = ""), "0")

        # compute the fitness for each state
        if ( fitness_function == "additive" ) {
          fitnesses <<- additiveFitnessFunction(model[selected_sites], lambda0, delta)
        } else if (fitness_function == "multiplicative" ) {
          fitnesses <<- multiplicativeFitnessFunction(model[selected_sites], lambda0, delta)
        } else {
          stop("Please choose either additive or multiplicative fitness function.")
        }
        fitnesses <<- c(fitnesses, 0)
        
        # recode the selected sites
        sel_seq <<- t(t(apply(sel_seq, 1, paste0, collapse = "")))

        # create the rate matrix
        S <<- makeYuleRateMatrix(model[selected_sites], fitnesses, phi, gamma)

      }

      if ( sel_lik_container_dirty == TRUE ) {

        # create the likelihood containers
        if (no_selected_sites) {
          likelihood_container <- matrix(0, nrow = 5, ncol = 1)
          rownames(likelihood_container) <- c("A","C","G","T","0")
        } else {
          likelihood_container <- matrix(0, nrow = length(sel_states), ncol = 1)
          rownames(likelihood_container) <- sel_states
        }
        sel_lik_container <<- vector("list", nrow(edge_matrix))

        # populate the containers
        tips <- tree$tip.label
        for (i in 1:length(tips)) {

          # get the tip
          this_tip <- tips[i]
          this_row <- edge_matrix$index[edge_matrix$desc == i]

          # initialize the container
          this_container <- likelihood_container

          # fill in conditional likelihood
          if (no_selected_sites) {
            if (this_tip %in% extant_tips) {
              this_container[,1] <- c(1,1,1,1,0)
            } else {
              this_container[,1] <- c(phi,phi,phi,phi,0)
            }
          } else {
            # get the data for this tip
            this_sel_data <- sel_seq[this_tip,]
            if (this_tip %in% extant_tips) {
              this_container[this_sel_data,1] <- 1.0
            } else {
              this_container[this_sel_data,1] <- phi
            }
          }

          # store the likelihood
          sel_lik_container[[this_row]] <<- this_container

        }

        # set clean
        sel_lik_container_dirty <<- FALSE

      }

    },

    computeNeutralTransitionProbabilities = function() {

      P_neu <<- vector("list", nrow(edge_matrix))
      Q <- matrix(gamma / 3, 4, 4)
      diag(Q) <- -gamma

      # eigen decomposition
      eigen <- eigen(Q)
      vals  <- eigen$values
      vec   <- eigen$vectors
      vecI  <- solve(vec)
      tpFunc <- function(t) computeTransitionProbability(vec, vecI, exp(t * vals))

      # populate transition probabilities per branch
      for(i in 1:nrow(edge_matrix)) {

        # get the branch length
        this_bl <- edge_matrix$bl[i]

        # store P
        P_neu[[i]] <<- tpFunc(this_bl)

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

    setGamma = function(gamma_) {
      gamma <<- gamma_
      neu_lik_dirty <<- TRUE
      sel_lik_dirty <<- TRUE
    },

    setDelta = function(delta_) {
      delta <<- delta_
      sel_lik_dirty <<- TRUE
    },

    setPhi = function(phi_) {
      phi <<- phi_
      sel_lik_dirty <<- TRUE
      sel_lik_container_dirty <<- TRUE
    },

    computeLikelihood = function() {

      # first, update the likelihoods
      computeNeutralLikelihood()
      computeSelectedLikelihood()

      # now, collect the appropriate values at the root
      ll <- getNeutralLikelihood() + getSelectedLikelihood()
      return(ll)

    },

    getNeutralLikelihood = function() {
      # only use the likelihoods for neutral sites
      return(sum(neutral_site_log_likelihoods[model == "-"]))
    },

    getSelectedLikelihood = function() {
      # only use the likelihoods for sites under selection
      return(selected_site_log_likelihood)
    },

    computeNeutralLikelihood = function() {

      if (neu_lik_dirty) {

        # update transition probabilities
        computeNeutralTransitionProbabilities()

        # compute the likelihood for each site
        neu_rescale_log <<- rep(0, num_sites)
        recursiveComputeConditionalLikelihoodNeutral(edge_matrix$index[1])

        # get likelihoods at the root
        neu_lik_root <- neu_lik_container[[1]]

        # compute per-site log likelihoods
        neutral_site_log_likelihoods <<- log(colMeans(neu_lik_root)) + neu_rescale_log

        # update the dirty flag
        neu_lik_dirty <<- FALSE

      }

    },

    computeSelectedLikelihood = function() {

      if (sel_lik_dirty) {

        # update containers
        initSelectedSitesContainers()

        # update transition probabilities
        computeSelectedTransitionProbabilities()

        # compute the likelihood for the selected site
        sel_rescale_log <<- 0
        recursiveComputeConditionalLikelihoodSelected(edge_matrix$index[1])

        # get likelihoods at the root
        sel_lik_root <- sel_lik_container[[1]]

        # get the transition probability
        this_P <- P_sel[[1]]

        # get the conditional likelihoods
        desc_CL <- sel_lik_container[[1]]

        # compute the conditional likelihood
        this_CL <- this_P %*% desc_CL

        # accumulate probability at the origin
        # excluding the
        # this_CL <- this_CL[-nrow(this_CL),]
        # selected_site_log_likelihood <<- log(mean(this_CL)) + sel_rescale_log
        selected_site_log_likelihood <<- log(mean(this_CL[fitnesses == lambda0])) + sel_rescale_log

        # update the dirty flag
        sel_lik_dirty <<- FALSE

      }

    },

    recursiveComputeConditionalLikelihoodNeutral = function(index) {

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
          recursiveComputeConditionalLikelihoodNeutral(descendants[i])
        }

        # compute the likelihood for the left descendants
        left_P       <- P_neu[[desc_index[1]]]
        left_CL      <- neu_lik_container[[desc_index[1]]]
        left_partial <- left_P %*% left_CL

        # compute the likelihood for the right descendants
        right_P       <- P_neu[[desc_index[2]]]
        right_CL      <- neu_lik_container[[desc_index[2]]]
        right_partial <- right_P %*% right_CL

        # compute the likelihood at the focal node
        this_CL <- left_partial * right_partial

        # rescale the likelihood
        scalar  <- colMax(this_CL)
        this_CL <- scaleMatrix(this_CL, scalar)

        # scalar <- numeric(num_sites)
        # for(i in 1:num_sites) {
        #   scalar[i] <- max(this_CL[,i])
        #   this_CL[,i] <- this_CL[,i] / scalar[i]
        # }

        neu_rescale_log <<- neu_rescale_log + log(scalar)

        # store the conditional likelihood
        neu_lik_container[[this_index]] <<- this_CL

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
        this_CL <- fitnesses * this_CL

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

fitGamma = function(phylo, aln, lambda0, gamma, delta, mu, phi) {

  # create the calculator
  nsites <- ncol(aln)
  model <- rep("-", nsites)
  calculator <- YuleLikelihood(phylo, aln, model, lambda0, gamma, delta, 0, phi)

  # create the objective function for gamma
  objective_gamma <- function(gamma) {

    # set gamma
    calculator$setGamma(gamma)

    # compute likelihood
    calculator$computeLikelihood()

    # get just the neutral likelihood
    neutral_likelihood <- calculator$getNeutralLikelihood()

    return(neutral_likelihood)

  }

  # estimate gamma
  fit_gamma <- optimize(objective_gamma, lower = 1e-16, upper = 1, maximum = TRUE)$maximum
  calculator$setGamma(fit_gamma)
  gamma_lik <- calculator$computeLikelihood()

  return(c(gamma = fit_gamma, lik = gamma_lik))

}

fitDelta = function(calculator, model) {

  # set the model
  calculator$setModel(model)

  # create the objective function for delta
  objective_delta <- function(delta) {

    # set gamma
    calculator$setDelta(delta)

    # compute likelihood
    likelihood <- calculator$computeLikelihood()

    return(likelihood)

  }

  fit_delta <- optimize(objective_delta, lower = 1e-16, upper = 1, maximum = TRUE)$maximum
  calculator$setDelta(fit_delta)
  delta_lik <- calculator$computeLikelihood()

  res <- data.frame(lik = delta_lik, gamma = calculator$gamma, delta = fit_delta)
  return(res)

}

fitModel = function(phylo, aln, model, lambda0, gamma, delta, mu, phi) {

  # create the calculator
  calculator <- YuleLikelihood(phylo, aln, model, lambda0, gamma, delta, 0, phi)

  # create the objective function for gamma
  objective_gamma <- function(gamma) {

    # set gamma
    calculator$setGamma(gamma)

    # compute likelihood
    calculator$computeLikelihood()

    # get just the neutral likelihood
    neutral_likelihood <- calculator$getNeutralLikelihood()

    return(neutral_likelihood)

  }

  # estimate gamma
  fit_gamma <- optimize(objective_gamma, lower = 1e-16, upper = 1, maximum = TRUE)$maximum
  calculator$setGamma(fit_gamma)
  gamma_lik <- calculator$computeLikelihood()

  if ( any(model != "-") == FALSE ) {

    # just fit the constant-rate model
    calculator$setDelta(0)
    fit_delta <- 0
    delta_lik <- calculator$computeLikelihood()

  } else {

    # create the objective function for delta
    objective_delta <- function(delta) {

      # set gamma
      calculator$setDelta(delta)

      # compute likelihood
      likelihood <- calculator$computeLikelihood()

      return(likelihood)

    }

    fit_delta <- optimize(objective_delta, lower = 1e-16, upper = 1, maximum = TRUE)$maximum
    calculator$setDelta(fit_delta)
    delta_lik <- calculator$computeLikelihood()

  }

  res <- data.frame(lik = delta_lik, gamma = fit_gamma, delta = fit_delta)
  return(res)

}
