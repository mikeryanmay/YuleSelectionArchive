# get the arguments
args <- commandArgs(TRUE)

# just the replicate number
indir <- args[1]

# source the code
library(cubature)
source("../scripts/countClassBackward.R", chdir = TRUE)
source("../scripts/YuleLikelihoodBinary.R", chdir = TRUE)
source("../scripts/adaptiveSimpsonsIntegrator.R", chdir = TRUE)

# find tree and data
tree_file <- list.files(indir, pattern = "tree.nex", full.names = TRUE)
seq_file  <- list.files(indir, pattern = "seq.nex", full.names = TRUE)

# if the files don't exist, abort
if ( file.exists(seq_file) == FALSE ) {
  cat("Couldn't find sequence file.\n")
  q()
}

##########################
# specify the true model #
##########################

# diversification parameters
r0      <- 2.5
phi     <- 1 / 7
lambda0 <- r0 * phi
factor  <- c(1.25, 1.5, 1.75, 2)

# mutation rate
genome_size <- 30000
mutation_rate_per_site_per_year <- 0.00084
mutation_rate_per_site_per_day  <- mutation_rate_per_site_per_year / 365
mu <- mutation_rate_per_site_per_day

# unknown parameters
tmp     <- strsplit(gsub("/", "_", indir), "_")[[1]]
factor  <- as.numeric(tmp[which(tmp == "sim") - 1])
model   <- "A"
lambda1 <- lambda0 * factor
day     <- as.numeric(tmp[which(tmp == "day") + 1])

##########################
# read the sequence data #
##########################

seq  <- read.nexus.data(seq_file)
seq  <- do.call(rbind, seq)
seq  <- seq[,1, drop = FALSE]

#################
# read the tree #
#################

# only read a tree if it exists
is_singleton <- FALSE
if ( nrow(seq) > 1 ) {
  
  # read the tree
  tree <- read.nexus(tree_file)
  
  # add a root edge
  stem_tree <- tree
  stem_tree$root.edge <- 0
  
} else {
  
  singleton <- list("tip.label" = "t_1", "edge.length" = day)
  stem_tree <- tree <- singleton
  is_singleton <- TRUE
  
}

# mutation rate
genome_size <- 30000
mutation_rate_per_site_per_year <- 0.00084
mutation_rate_per_site_per_day  <- mutation_rate_per_site_per_year / 365
gamma <- mutation_rate_per_site_per_day

# diversification rates
r0      <- 2.5
phi     <- 1 / 7
lambda0 <- r0 * phi

#####################
# specify the prior #
#####################

# uniform prior
# between 1 and 4 (always positive)
lower_int <- 1
upper_int <- 4
prior     <- function(f) dunif(f, min = lower_int, max = upper_int, log = TRUE)

# origin times and probabilities
origin_times <- 1:14 - 1 + day # time in the past
origin_probs <- dunif(0:13, min = 0, max = 13)

######################
# fit the tree model #
######################

# make the model
tree_model <- YuleLikelihoodBinary(stem_tree, seq, model, lambda0, lambda1, gamma / 3, gamma, phi, initial_state_ = "1", is_singleton)

# make the likelihood function
tree_log_likelihood <- function(f) {
  
  # compute the value of lambda1
  new_lambda1 <- lambda0 * f
  
  # set value
  tree_model$setLambda1(new_lambda1)
  
  # compute the likelihood to the root
  ll <- tree_model$computeLikelihoodIntegrated(origin_times, origin_probs)
  
  # return likelihood
  return(ll)
  
}

# make the unnormalized posterior
tree_log_posterior_un <- function(x) {
  pp <- sapply(x, tree_log_likelihood) + prior(x)
  return(pp)
}

# find the MAP estimate
opt      <- optimize(tree_log_posterior_un, lower = lower_int, upper = upper_int, maximum = TRUE)
scalar   <- opt$objective
tree_posterior_mode <- opt$maximum

# create the posterior, scaled by an arbitrary factor
tree_posterior <- function(x, scalar) {
  pp <- exp(sapply(x, tree_log_likelihood) + prior(x) - scalar)
  return(pp)
}

# compute the marginal probability
tree_int <- adaptiveSimpsonsIntegrator(tree_posterior, lower = lower_int, upper = upper_int, scalar = scalar, tol = 1e-4)
tree_marginal_probability <- scalar + log(tree_int$integral)

# compute the posterior summaries
tree_posterior_mean  <- posteriorMean(tree_int)
tree_posterior_var   <- posteriorVariance(tree_int, tree_posterior_mean)
tree_posterior_quant <- posteriorQuantile(tree_int, factor, scalar = scalar)
tree_posterior_error <- posteriorVariance(tree_int, factor)

# savage-dickey ratio against neutral
prior_log_density_neutral     <- prior(1)
posterior_log_density_neutral <- log(tree_posterior(1, scalar = scalar))
tree_SDR_neutral              <- 2 * (posterior_log_density_neutral - prior_log_density_neutral)

#######################
# fit the count model #
#######################

# make the model
count_model <- countModelBackward(tree, seq, model, lambda0, lambda1, phi, gamma / 3, gamma, initial_state_ = "1", singleton_ = is_singleton, stem_ = TRUE)

# make the likelihood function
count_log_likelihood <- function(f, atol = 1e-6, rtol = 1e-6) {
  
  cat(f, "\t", sep ="")
  
  # compute the value of lambda1
  lambda1 <- lambda0 * f
  
  # set value
  count_model$setLambda1(lambda1)
  
  # compute likelihood
  ll <- count_model$computeLikelihood(origin_times, origin_probs, verbose = FALSE, atol = atol, rtol = rtol)
  
  # inter
  
  # cat
  cat(ll, "\n")
  
  # return likelihood
  return(ll)
  
}

# make the log unnormalized posterior
count_log_posterior_un <- function(x) {
  pp <- sapply(x, count_log_likelihood) + prior(x)
  return(pp)
}

# find the MAP estimate
opt    <- optimize(count_log_posterior_un, lower = lower_int, upper = upper_int, maximum = TRUE)
scalar <- opt$objective
count_posterior_mode <- opt$maximum

# create the posterior, scaled
count_posterior <- function(x, scalar) {
  pp <- exp(sapply(x, count_log_likelihood) + prior(x) - scalar)
  return(pp)
}

# compute the marginal probability
count_int <- adaptiveSimpsonsIntegrator(count_posterior, lower = lower_int, upper = upper_int, scalar = scalar, tol = 1e-4)
count_marginal_probability <- scalar + log(count_int$integral)

# compute the posterior summaries
count_posterior_mean  <- posteriorMean(count_int)
count_posterior_var   <- posteriorVariance(count_int, count_posterior_mean)
count_posterior_quant <- posteriorQuantile(count_int, factor, scalar = scalar)
count_posterior_error <- posteriorVariance(count_int, factor)

# savage-dickey ratio against neutral
prior_log_density_neutral     <- prior(1)
posterior_log_density_neutral <- log(count_posterior(1, scalar = scalar))
count_SDR_neutral             <- 2 * (posterior_log_density_neutral - prior_log_density_neutral)

##################
# neutral models #
##################

tree_model$setLambda1(lambda0)
tree_model_neutral_likelihood <- tree_model$computeLikelihoodIntegrated(origin_times, origin_probs)

count_model$setLambda1(lambda0)
count_model_neutral_likelihood <- count_model$computeLikelihood(origin_times, origin_probs, verbose = FALSE, atol = 1e-6, rtol = 1e-6)

##########
# finish #
##########

# make the output file
outfile <- paste0(indir, "/count_model_integrated_summary_noroot.tsv")

# make the output table
res <- data.frame(tree_marginal_probability, 
                  tree_model_neutral_likelihood,
                  tree_posterior_mode, tree_posterior_mean, tree_posterior_var, tree_posterior_quant, tree_posterior_error,
                  count_marginal_probability, 
                  count_model_neutral_likelihood,
                  count_posterior_mode, count_posterior_mean, count_posterior_var, count_posterior_quant, count_posterior_error)

# write to file
write.table(res, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)

# write points
tree_pts  <- data.frame(points = tree_int$points, values = tree_int$values)
count_pts <- data.frame(points = count_int$points, values = count_int$values)
write.table(tree_pts, file = paste0(indir, "/count_model_integrated_tree_points.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(count_pts, file = paste0(indir, "/count_model_integrated_count_points.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

# plot posteriors
ymax <- max(c(tree_int$values / tree_int$integral, count_int$values / count_int$integral))

pdf(paste0(indir, "/count_model_plots.pdf"), height = 4)
par(mar = c(4,4,0,0) + 0.1)
plot(x   = lambda0 * tree_int$points, y = tree_int$values / tree_int$integral, ylim = c(0, ymax), type = "b", col = "blue", ylab = "posterior probability", xlab = "lambda1")
lines(x  = lambda0 * count_int$points, y = count_int$values / count_int$integral, type = "b", col = "red")
abline(v = lambda0 * factor, lty = 2)
abline(v = lambda0 * tree_posterior_mode, col = "blue", lty = 2)
abline(v = lambda0 * count_posterior_mode, col = "red", lty = 2)
dev.off()

# quit
q()









