# get the arguments
args <- commandArgs(TRUE)

# just the replicate number
indir <- args[1]

# source the code
library(cubature)
source("../scripts/likelihood.R")

# find tree and data
tree_file <- list.files(indir, pattern = "tree.nex", full.names = TRUE)
seq_file  <- list.files(indir, pattern = "seq.nex", full.names = TRUE)

# if the files don't exist, abort
if ( file.exists(tree_file) == FALSE ) {
  cat("Couldn't find tree file.\n")
  q()
}

if ( file.exists(seq_file) == FALSE ) {
  cat("Couldn't find sequence file.\n")
  q()
}

# read the data
tree <- read.nexus(tree_file)
seq  <- read.nexus.data(seq_file)
seq  <- do.call(rbind, seq)
num_sites <- ncol(seq)

# mutation rate
genome_size <- 30000
mutation_rate_per_site_per_year <- 0.00084
mutation_rate_per_site_per_day  <- mutation_rate_per_site_per_year / 365
gamma <- mutation_rate_per_site_per_day

# diversification rates
r0      <- 2.5
phi     <- 1/7
lambda0 <- r0 * phi

###################
# model with data #
###################

# create the model and calculator
model <- rep("-", num_sites)
calculator <- YuleLikelihood(tree, seq, model, lambda0, gamma, 1.5, 0, phi, fitness_function_ = "multiplicative")
calculator$computeLikelihood()

########################
# enumerate the models #
########################

# all models with a single site/allele under selection
num_sites             <- ncol(seq)
num_selected_sites    <- 1
sites_under_selection <- t(combn(1:num_sites, num_selected_sites))
selected_states       <- expand.grid(lapply(1:num_selected_sites, function(x) c("A","C","G","T") ), stringsAsFactors = FALSE)
number_of_models      <- nrow(sites_under_selection) * nrow(selected_states)

models <- vector("list", number_of_models)
blank_model <- rep("-", num_sites)
k <- 0
for(i in 1:nrow(sites_under_selection)) {
  for(j in 1:nrow(selected_states)) {
    
    # create the model
    this_model <- blank_model
    this_model[sites_under_selection[i,]] <- unlist(selected_states[j,])
    
    # save the model
    k <- k + 1
    models[[k]] <- this_model
    
  }
}
models <- do.call(rbind, models)

######################
# make the functions #
######################

# create the likelihood function 
likelihood <- function(f) {
  
  # compute the delta
  delta <- f
  
  # set the value
  calculator$setDelta(delta)
  
  # compute the log-likelihood
  ll <- calculator$computeLikelihood()
  
  return(ll)
  
}

# create the prior distribution

# uniform prior
# between 1 and 4 (always positive)
prior <- function(f) dunif(f, min = 1, max = 6, log = TRUE)
lower_int <- 1
upper_int <- 4

# log posterior, unnormalize, vectorized
log_posterior_un <- function(x) {
  pp <- sapply(x, likelihood) + prior(x)
  return(pp)
}

# normalized posterior, not vectorized
posterior <- function(x, scalar) {
  pp <- exp(likelihood(x) + prior(x) - scalar)
  return(pp)
}

#######################
# iterate over models #
#######################

# first, selection models
model_likelihoods <- numeric(number_of_models)
posterior_mode    <- numeric(number_of_models)
convergence       <- numeric(number_of_models)
bar <- txtProgressBar(style = 3, width = 40)
for(i in 1:number_of_models) {
  
  # get the model
  this_model <- models[i,]
  
  # set the model
  calculator$setModel(this_model)
  
  # find the scalar
  opt    <- optimize(log_posterior_un, lower = lower_int, upper = upper_int, maximum = TRUE)
  scalar <- opt$objective
  delta  <- opt$maximum
  
  # compute the marginal probability
  int <- hcubature(posterior, lowerLimit = lower_int, upperLimit = upper_int, scalar = scalar)
  marginal_probability <- scalar + log(int$integral)
  
  # write down
  model_likelihoods[i] <- marginal_probability
  posterior_mode[i]    <- delta
  convergence[i]       <- int$returnCode
  
  setTxtProgressBar(bar, i / number_of_models)
  
}

close(bar)

#######################
# now, constant model #
#######################

calculator$setDelta(1.0)
marginal_probability <- calculator$computeLikelihood()

model_likelihoods <- c(model_likelihoods, marginal_probability)
posterior_mode    <- c(posterior_mode, 1)
convergence       <- c(convergence, 0.0)

###########################################
# posterior quantities for the true model #
###########################################

# get the true value
tmp <- strsplit(gsub("/", "_", indir), "_")[[1]]
true_f <- as.numeric(tmp[which(tmp == "rep") - 1])

# reset the model
calculator$setModel(models[1,])

# compute the marginal likelihood
marginal_likelihood <- model_likelihoods[1]

# compute the expected value
int <- hcubature( function(x) {
  x * exp(log_posterior_un(x) - marginal_likelihood)
  
}, lowerLimit = lower_int, upperLimit = upper_int)

posterior_mean <- int$integral

# compute the variance/sd
int <- hcubature( function(x) {
  ((x - posterior_mean)^2) * exp(log_posterior_un(x) - marginal_likelihood)
}, lowerLimit = lower_int, upperLimit = upper_int)

posterior_variance <- int$integral
posterior_sd       <- sqrt(posterior_variance)

# compute the quantile of the true value
int <- hcubature( function(x) {
  exp(log_posterior_un(x) - marginal_likelihood)
}, lowerLimit = lower_int, upperLimit = true_f)

quantile <- int$integral
is_contained <- as.numeric(quantile > 0.025 & quantile < 0.975)

# compute the squared loss
int <- hcubature( function(x) {
  ((x - true_f)^2) * exp(log_posterior_un(x) - marginal_likelihood)
}, lowerLimit = lower_int, upperLimit = upper_int)

expected_squared_loss <- int$integral

#######################
# write stuff to file #
#######################

# make the output file
outfile <- paste0(indir, "/site_models.tsv")

# make the output table
model_names <- c(paste0(which(models != "-", arr.ind = TRUE)[,2], models[models != "-"]), "neutral")
model_liks  <- data.frame(model = model_names, lik = model_likelihoods, mode = posterior_mode, convergence = convergence)

# write to file
write.table(model_liks, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)

# make posterior table
outfile <- paste0(indir, "/true_model.tsv")
posterior_properties <- data.frame(posterior_mean  = posterior_mean,
                                   posterior_var   = posterior_variance,
                                   posterior_sd    = posterior_sd,
                                   posterior_quant = quantile,
                                   is_contained    = is_contained,
                                   PESL            = expected_squared_loss)
write.table(posterior_properties, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)

# make a plot
pdf(paste0(indir, "/model_plot.pdf"), width = 12)
par(mar=c(4,4,0,0) + 0.1)
plot(model_liks$lik, pch = 19, xaxt = "n", xlab = NA)
axis(1, at = 1:nrow(model_liks), label = model_liks$model, las = 2)
dev.off()

# quit
q()





