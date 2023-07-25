library(cubature)
library(parallel)
source("../scripts/countClass.R", chdir = TRUE)
source("../scripts/YuleLikelihoodBinary.R", chdir = TRUE)

# settings
r0      <- 2.5
phi     <- 0.2
lambda0 <- 1.0
lambda1 <- 1.5 * lambda0
gamma   <- 0.5
nmax    <- 3

# the data
t0 <- 0   # start time
t1 <- 0.5 # time of first sample
t2 <- 1.0 # time of second sample
t3 <- 1.5 # time of third sample

x0 <- "0" # state of the original lineage
x1 <- "1" # state of sample at time t1
x2 <- "0" # state of sample at time t2
x3 <- "0" # state of sample at time t3

################################
# likelihood under count model #
################################

count_likelihood <- function(nmax, lambda0, lambda1, phi, gamma) {
  
  # make the matrices
  Q  <- countMatrixModel(nmax, lambda0, lambda1, phi, gamma / 3, gamma)
  R1 <- sampleMatrixModel(nmax, phi, x1)
  R2 <- sampleMatrixModel(nmax, phi, x2)
  
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
  
  # get probability of present sample
  end_state <- c("0" = 0, "1" = 0)
  end_state[x3] <- 1
  end_state <- paste0("(", end_state[1], ",", end_state[2], ")")
  lik_count <- as.numeric(p[end_state])
  
  return(lik_count)
  
}

lik_count <- count_likelihood(nmax, lambda0, lambda1, phi, gamma)

###############################
# likelihood under tree model #
###############################

# data matrix
data <- matrix(c(x1, x2, x3), ncol = 1)
rownames(data) <- c("x1", "x2", "x3")
data <- ifelse(data == 1, "A", "C") # convert to arbitrary nucleotides

# ages to phylo
tree_1_ages_to_phylo <- function(a1, a2) {
  
  # base newick string
  newick <- "((x3:T3, x2:T2):T4, x1:T1):T5;"
  
  # compute branch lengths
  T1 <- t1 - a1
  T2 <- t2 - a2
  T3 <- t3 - a2
  T4 <- a2 - a1
  T5 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
 
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
   
}

tree_2_ages_to_phylo <- function(a1, a2) {
  
  # base newick string
  newick <- "((x2:T2, x1:T1):T4, x3:T3):T5;"
  
  # compute branch lengths
  T1 <- t1 - a2
  T2 <- t2 - a2
  T3 <- t3 - a1
  T4 <- a2 - a1
  T5 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

tree_3_ages_to_phylo <- function(a1, a2) {
  
  # base newick string
  newick <- "((x3:T3, x1:T1):T4, x2:T2):T5;"
  
  # compute branch lengths
  T1 <- t1 - a2
  T2 <- t2 - a1
  T3 <- t3 - a2
  T4 <- a2 - a1
  T5 <- a1
  
  # substitute branch lengths
  newick <- gsub("T1", T1, newick)
  newick <- gsub("T2", T2, newick)
  newick <- gsub("T3", T3, newick)
  newick <- gsub("T4", T4, newick)
  newick <- gsub("T5", T5, newick)
  
  # make phylo
  phylo <- read.tree(text = newick)
  return(phylo)
  
}

# make initial trees
tree_1 <- tree_1_ages_to_phylo(0.2, 0.3)
tree_2 <- tree_2_ages_to_phylo(0.2, 0.3)
tree_3 <- tree_3_ages_to_phylo(0.2, 0.3)

# par(mar = c(0,0,0,0), mfrow = c(1,3))
# plot(tree_1, root.edge = TRUE, direction = "up")
# plot(tree_2, root.edge = TRUE, direction = "up")
# plot(tree_3, root.edge = TRUE, direction = "up")

# create the calculator
tree_1_likelihood_calculator <- YuleLikelihoodBinary(tree_1, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)
tree_2_likelihood_calculator <- YuleLikelihoodBinary(tree_2, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)
tree_3_likelihood_calculator <- YuleLikelihoodBinary(tree_3, data, "A", lambda0, lambda1, gamma / 3, gamma, phi)

# create the likelihood functions as a function of t
tree_1_likelihood_time <- function(a1, a2) {
  
  # set the tree
  tree_1_likelihood_calculator$setTree(tree_1_ages_to_phylo(a1, a2))
  
  # return the likelihood
  exp(tree_1_likelihood_calculator$computeLikelihood())
  
}

tree_2_likelihood_time <- function(a1, a2) {
  
  # set the tree
  tree_2_likelihood_calculator$setTree(tree_2_ages_to_phylo(a1, a2))
  
  # return the likelihood
  exp(tree_2_likelihood_calculator$computeLikelihood())
  
}

tree_3_likelihood_time <- function(a1, a2) {
  
  # set the tree
  tree_3_likelihood_calculator$setTree(tree_3_ages_to_phylo(a1, a2))
  
  # return the likelihood
  exp(tree_3_likelihood_calculator$computeLikelihood())
  
}

# tree_1_likelihood_time(0.2, 0.3)
# tree_2_likelihood_time(0.2, 0.3)
# tree_3_likelihood_time(0.2, 0.3)

# create the integrators
tree_1_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-10) {
  
  # set parameters
  tree_1_likelihood_calculator$setLambda0(lambda0)
  tree_1_likelihood_calculator$setLambda1(lambda1)
  tree_1_likelihood_calculator$setGamma01(gamma / 3)
  tree_1_likelihood_calculator$setGamma10(gamma)
  tree_1_likelihood_calculator$setPhi(phi)
  
  int <- integral2(function(x, y) {
    
    # cat("*")
    z <- x
    for(i in 1:ncol(x)) {
      for(j in 1:nrow(x)) {
        z[i,j] <- tree_1_likelihood_time(x[i,j], y[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t2, reltol = tol, abstol = tol)$Q

  return(int)
    
}

tree_2_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-10) {
  
  # set parameters
  tree_2_likelihood_calculator$setLambda0(lambda0)
  tree_2_likelihood_calculator$setLambda1(lambda1)
  tree_2_likelihood_calculator$setGamma01(gamma / 3)
  tree_2_likelihood_calculator$setGamma10(gamma)
  tree_2_likelihood_calculator$setPhi(phi)
  
  int <- integral2(function(x, y) {
    
    # cat("*")
    z <- x
    for(i in 1:ncol(x)) {
      for(j in 1:nrow(x)) {
        z[i,j] <- tree_2_likelihood_time(x[i,j], y[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, reltol = tol, abstol = tol)$Q
  
  return(int)
  
}

tree_3_likelihood <- function(lambda0, lambda1, phi, gamma, tol = 1e-10) {
  
  # set parameters
  tree_3_likelihood_calculator$setLambda0(lambda0)
  tree_3_likelihood_calculator$setLambda1(lambda1)
  tree_3_likelihood_calculator$setGamma01(gamma / 3)
  tree_3_likelihood_calculator$setGamma10(gamma)
  tree_3_likelihood_calculator$setPhi(phi)
  
  int <- integral2(function(x, y) {
    
    # cat("*")
    z <- x
    for(i in 1:ncol(x)) {
      for(j in 1:nrow(x)) {
        z[i,j] <- tree_3_likelihood_time(x[i,j], y[i,j])
      }
    }
    
    return(z)
    
  }, xmin = 0, xmax = t1, ymin = function(x) x, ymax = t1, reltol = tol, abstol = tol)$Q
  
  return(int)
  
}

# combined likelihood integral
tree_likelihood <- function(lambda0, lambda1, phi, gamma) {
  tree_1_likelihood(lambda0, lambda1, phi, gamma) + tree_2_likelihood(lambda0, lambda1, phi, gamma) + tree_3_likelihood(lambda0, lambda1, phi, gamma) 
}

# compute likelihoods
lik_tree <- tree_likelihood(lambda0, lambda1, phi, gamma)

# compare
lik_tree; lik_count
lik_tree / lik_count
sprintf("%.15f", lik_tree - lik_count)

# log(lik_tree)
# log(lik_count)


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
pdf("figures/validate_likelihood_3_samples.pdf", height = 4, width = 8)
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
write.table(res, file = "results/compare_likelihoods_3_samples.tsv", col.names = TRUE, quote = FALSE, row.names = FALSE)



