source("../scripts/countClass.R", chdir = TRUE)
source("../scripts/countClassBackward.R", chdir = TRUE)

# dataset
tree <- read.nexus("../simulation_2/sims/tips_100_f_1.5/rep_1/tree.nex")
root <- max(node.depth.edgelength(tree))
tree$root.edge <- 0
seq  <- read.nexus.data("../simulation_2/sims/tips_100_f_1.5/rep_1/seq.nex")
seq  <- do.call(rbind, seq)
seq  <- seq[,1, drop = FALSE]

# mutation rate
genome_size <- 30000
mutation_rate_per_site_per_year <- 0.00084
mutation_rate_per_site_per_day  <- mutation_rate_per_site_per_year / 365
gamma <- mutation_rate_per_site_per_day

# diversification rates
r0      <- 2.5
phi     <- 1 / 7
lambda0 <- r0 * phi
lambda1 <- lambda0 * 1.2

# make the models
backward_model <- countModelBackward(tree, seq, "A", lambda0, lambda1, phi, gamma / 3, gamma)
forward_model  <- countModel(tree, seq, "A", lambda0, lambda1, phi, gamma / 3, gamma)

# iterate over lambda1s
nlambda <- 100
lambda1 <- seq(0, 2, length.out = nlambda + 1)[-1]

backward_lnl <- numeric(nlambda)
forward_lnl  <- numeric(nlambda)

bar <- txtProgressBar(style = 3, width = 40)
for(i in 1:nlambda) {
  
  # get the lambda
  this_lambda1 <- lambda1[i]
 
  # set lambda
  backward_model$setLambda1(this_lambda1)
  forward_model$setLambda1(this_lambda1)

  # compute likelihoods
  backward_lnl[i] <- backward_model$computeLikelihood(root, 1, atol = 1e-6, rtol = 1e-6, verbose = FALSE)
  forward_lnl[i]  <- forward_model$computeLikelihood(atol = 1e-6, rtol = 1e-6, verbose = FALSE) + log(lambda0)
  
  setTxtProgressBar(bar, i / nlambda)
  
}
close(bar)

# plot likelihoods
dir.create("figures", showWarnings = FALSE)
pdf("figures/validate_likelihood_forward_vs_backward.pdf", height = 5, width = 10)
par(mar = c(4,4,0,0) + 0.1)
plot(lambda1,   backward_lnl, pch = 3, xlab = "f", ylab = "log likelihood", las = 1)
points(lambda1, forward_lnl,  pch = 4)
legend("bottomright", legend = c("up-pass", "down-pass"), pch = c(4,3), bty = "n", title = "method")
dev.off()




