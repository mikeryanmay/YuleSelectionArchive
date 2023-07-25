library(cubature)
library(parallel)
source("../scripts/countClass.R", chdir = TRUE)
source("../scripts/YuleLikelihoodBinary.R", chdir = TRUE)

# settings
lambda0 <- 1
lambda1 <- 2.5 * lambda0
phi     <- 0.1
gamma   <- 0.5
nmax    <- 2

# the data
t0 <- 0   # start time
t1 <- 0.5 # time of sample
t2 <- 1.0 # end time
x0 <- "0" # state of the original lineage
x1 <- "1" # state of sample at time t1
x2 <- "0" # state of sample at time t2

################################
# likelihood under count model #
################################

count_likelihood <- function(nmax, lambda0, lambda1, phi, gamma) {
  
  # make the matrices
  Q <- countMatrixModel(nmax, lambda0, lambda1, phi, gamma / 3, gamma)
  R <- sampleMatrixModel(nmax, phi, x1)
  
  # choose the start state
  start_state <- c("0" = 0, "1" = 0)
  start_state[x0] <- 1
  start_state <- paste0("(", start_state[1], ",", start_state[2], ")")
  
  # make init p
  p <- numeric(R$num_states)
  names(p) <- R$labels
  p[start_state] <- 1
  
  # integrate to time t1
  p <- Q$solve(p, t1)
  
  # do the sample event
  p <- R$doSampleEvent(p)
  
  # integrate to time t2
  p <- Q$solve(p, t2 - t1)
  
  # determine the end state
  end_state <- c("0" = 0, "1" = 0)
  end_state[x2] <- 1
  end_state <- paste0("(", end_state[1], ",", end_state[2], ")")
  
  # get probability of present sample
  lik_count <- as.numeric(p[end_state])
 
  return(lik_count)
   
}

lik_count <- count_likelihood(nmax, lambda0, lambda1, phi, gamma)

###############################
# likelihood under tree model #
###############################

# data matrix
data <- matrix(c(x1, x2), ncol = 1)
rownames(data) <- c("x1", "x2")
data <- ifelse(data == 1, "A", "C") # convert to arbitrary nucleotides

tree_1_ages_to_phylo <- function(a1) {
 
  # base newick string
  newick <- "(x1:T1, x2:T2):T3;"
  
  # compute branch lengths
  T1 <- t1 - a1
  T2 <- t2 - a1
  T3 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

tree_1 <- tree_1_ages_to_phylo(0.2)

# create the calculator
tree_likelihood_calculator <- YuleLikelihoodBinary(tree_1, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)
tree_likelihood_calculator$computeLikelihood()

# create the likelihood function as a function of t
tree_likelihood_time <- function(t) {
  
  # set the tree
  tree_likelihood_calculator$setTree(tree_1_ages_to_phylo(t))
  
  # return the likelihood
  tree_likelihood_calculator$computeLikelihood()
  
}

# create the integrated likelihood function
tree_likelihood <- function(lambda0, lambda1, phi, gamma) {
 
  # set parameters
  tree_likelihood_calculator$setLambda0(lambda0)
  tree_likelihood_calculator$setLambda1(lambda1)
  tree_likelihood_calculator$setGamma01(gamma / 3)
  tree_likelihood_calculator$setGamma10(gamma)
  tree_likelihood_calculator$setPhi(phi)
 
  # integrate
  hcubature(function(t) {
    exp(tree_likelihood_time(t))
  }, lower = 0, upper = t1)$integral
    
}

# compute the likelihood once
lik_tree <- tree_likelihood(lambda0, lambda1, phi, gamma)

# compare
lik_tree; lik_count
lik_tree / lik_count

# iterate over values of lambda1
lambda1   <- seq(0, 2, length.out = 101)
lik_count <- numeric(length(lambda1))
lik_tree  <- numeric(length(lambda1))
bar <- txtProgressBar(style = 3, width = 40)
for(i in 1:length(lambda1)) {
  this_lambda1 <- lambda1[i]
  lik_count[i] <- count_likelihood(nmax, lambda0, this_lambda1, phi, gamma)
  lik_tree[i]  <- tree_likelihood(lambda0, this_lambda1, phi, gamma)
  setTxtProgressBar(bar, i / length(lambda1))
}

# plot
dir.create("figures", showWarnings = FALSE)
pdf("figures/validate_likelihood_2_samples.pdf", height = 4, width = 8)
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
write.table(res, file = "results/compare_likelihoods_2_samples.tsv", col.names = TRUE, quote = FALSE, row.names = FALSE)








