library(viridis)
library(RColorBrewer)
library(parallel)
library(abind)
source("../scripts/countClass.R", chdir = TRUE)

# parameters
lambda0 <- 0.5
lambda1 <- 1.0
phi     <- 0.1
gamma   <- 0.5
nmax    <- 9

# times
t  <- 2
ts <- seq(0, t, length.out = 11)[-1]
tx <- seq(0, t, 0.0001)
dt <- tx[2] - tx[1]
nt <- length(tx)

########################
# numerical integrator #
########################

# make the matrix
Q <- countMatrixModel2(nmax, lambda0, lambda1, phi, gamma / 3, gamma)

# make the initial probability vector
p <- numeric(Q$num_states)
names(p) <- Q$labels
p["(1,0)"] <- 1

# make the container
int_p <- matrix(NA, nrow = nt, ncol = length(p))
colnames(int_p) <- Q$labels
int_p[1,] <- p

# integrate over time
bar <- txtProgressBar(style = 3, width = 40)
for(i in 2:nt) {
  
  # integrate one time step
  p <- Q$solve(p, dt)
  # p <- Q$solve(p, dt, method = "Higham08", atol = 1e-12, rtol = 1e-12)
  
  # store result
  int_p[i,] <- p
  
  setTxtProgressBar(bar, i / nt)
  
}

# compute marginal probabilities for each state per time
num_probs_0 <- do.call(cbind, lapply(0:nmax, function(x) rowSums(int_p[,Q$states[,1] == x, drop = FALSE])))
num_probs_1 <- do.call(cbind, lapply(0:nmax, function(x) rowSums(int_p[,Q$states[,2] == x, drop = FALSE])))

# compute the joint probabilities per state at the end time
num_joint <- matrix(0, nmax + 1, nmax + 1)
num_joint[upper.tri(num_joint, diag = TRUE)[,(nmax + 1):1]] <- c(0, int_p[nt,])

###########################
# monte carlo simulations #
###########################

# simulation function
simForward <- function(nmax, lambda0, lambda1, phi, gamma, init, ts) {
  
  # ts is times to observe counts
  num_times <- length(ts)
  
  # initialize the return value
  sim <- matrix(NA, nrow = num_times + 1, ncol = 2)
  sim[1,] <- init

  # counters
  current_time  <- 0
  current_state <- init
  
  # for each time interval
  for(i in 1:num_times) {
    
    # get the end time
    end_time <- ts[i]
    
    # simulate forward 
    repeat {
      
      # compute the rate
      event_rates <- c(current_state[1] * lambda0, current_state[2] * lambda1, current_state[1] * gamma / 3, current_state[2] * gamma, sum(current_state * phi))
      total_rate  <- sum(event_rates)
      
      # draw a waiting time
      current_time <- current_time + rexp(1, total_rate)
      
      # terminate if we go over the end time
      if ( current_time > end_time ) {
        current_time <- end_time
        break
      }
      
      # choose the type of event
      event_type <- sample.int(5, size = 1, prob = event_rates / total_rate)
      
      # cat("event type ", event_type, " at time ", current_time, "\n", sep ="")
      
      # do the event
      if ( event_type == 1 ) {
        # birth of type 0
        
        # reject simulations that go beyond the allowed number of taxa
        if ( sum(current_state) >= nmax ) {
          return(sim)
        } else {
          current_state[1] <- current_state[1] + 1  
        }
        
      } else if ( event_type == 2 ) {
        # birth of type 1
        
        # reject simulations that go beyond the allowed number of taxa
        if ( sum(current_state) >= nmax ) {
          return(sim)
        } else {
          current_state[2] <- current_state[2] + 1
        }
        
      } else if ( event_type == 3 ) {
        # mutation from 0 to 1
        current_state[1] <- current_state[1] - 1
        current_state[2] <- current_state[2] + 1
      } else if ( event_type == 4 ) {
        # mutation from 1 to 0
        current_state[1] <- current_state[1] + 1
        current_state[2] <- current_state[2] - 1
      } else if ( event_type == 5) {
        # sampling event
        # terminate early
        return(sim)
      }
      
    }
    
    # record the state
    sim[i + 1,] <- current_state
    
  }
  
  # return
  return(sim)
  
}

# simForward(lambda0, lambda1, phi, gamma, c(5,5), ts)

# do the Monte Carlo simulation
reps <- 10000000

# sims <- mclapply(1:reps, function(x) simForward(nmax, lambda0, lambda1, phi, gamma, c(1,0), ts), mc.cores = 6)
# sims <- abind(sims, along = 3)

sims <- array(dim = c(length(ts) + 1, 2, reps))
bar  <- txtProgressBar(style = 3, width = 40)
for(i in 1:reps) {
  sims[,,i] <- simForward(nmax, lambda0, lambda1, phi, gamma, c(1,0), ts)
  setTxtProgressBar(bar, i / reps)
}

# for each time point, get the distribution of the number of each type
ts0 <- c(0, ts)
sim_probs_0 <- matrix(NA, nrow = length(ts0), ncol = nmax + 1)
sim_probs_1 <- matrix(NA, nrow = length(ts0), ncol = nmax + 1)
for(i in 1:length(ts0)) {
  sim_probs_0[i,] <- tabulate(sims[i,1,] + 1, nbins = nmax + 1) / reps
  sim_probs_1[i,] <- tabulate(sims[i,2,] + 1, nbins = nmax + 1) / reps
}

# for the end time, compute the joint probability of each state
end_states <- apply(sims[length(ts0),,], 2, function(x) paste0("(", paste0(x, collapse =","), ")"))
sim_end_state_freq <- table(end_states)[Q$labels] / reps

# compute the joint probabilities per state at the end time
sim_joint <- matrix(0, nmax + 1, nmax + 1)
sim_joint[upper.tri(sim_joint, diag = TRUE)[,(nmax + 1):1]] <- c(0, sim_end_state_freq)

########
# plot #
########

dir.create("figures", showWarnings = FALSE)
cols <- viridis(nmax + 1, begin = 0.1)

# marginal probabilities
pdf("figures/validate_ode.pdf", height = 7, width = 10)
par(mar = c(0,4,0,0), oma = c(4,0,1,0)+0.1, mfrow = c(2,1))
matplot(x = ts0, y = sim_probs_0, ylim = c(10 / reps, 1), type = "p", pch = 19, cex = 0.75, col = cols, xaxt = "n", yaxt = "n", xlab = NA, ylab = bquote("ln"~p(n[0] == i)), log = "y")
matplot(x = tx, y = num_probs_0, type = "l", pch = 19, col = cols, lty = 1, add = TRUE)
axis(2, lwd = 0, lwd.tick = 1, las = 1)
legend("bottomright", legend = c(0:nmax), lty = 1, col = cols, bty = "n", pch = 19, pt.cex = 0.75, title = "i", ncol = 2)
legend("bottom", legend = c("Monte Carlo", "ODE"), pch = c(19, NA), lty = c(NA, 1), col = cols, bty = "n", pt.cex = 0.75, title = "method")

matplot(x = ts0, y = sim_probs_1, ylim = c(10 / reps, 1), type = "p", pch = 19, cex = 0.75, col = cols, lty = 1, xaxt = "n", yaxt = "n", xlab = "time", ylab = bquote("ln"~p(n[1] == j)), log = "y")
matplot(x = tx, y = num_probs_1, type = "l",lty = 1, pch = 19, col = cols, add = TRUE)
axis(1, lwd = 0, lwd.tick = 1)
axis(2, lwd = 0, lwd.tick = 1, las = 1)
legend("bottomright", legend = c(0:nmax), lty = 1, col = cols, bty = "n", pch = 19, pt.cex = 0.75, title = "j", ncol = 2)
mtext("time", side = 1, outer = TRUE, line = 2.5)
dev.off()

# joint probability
x_coords <- rep(1:(nmax + 1), nmax + 1) - 1
y_coords <- rep(1:(nmax + 1), each = nmax + 1) - 1
scale    <- 25
# cols     <- brewer.pal(3, "Set1")[2:1]
cols     <- turbo(2, begin = 0.1, end = 0.7)

pdf("figures/validate_ode_joint.pdf")
par(mar=c(2,2,2,2)+0.5)
plot(x = y_coords, y = x_coords, cex = scale * as.vector(num_joint), pch = 19, xlim = c(-1, nmax + 2), ylim = c(-1, nmax + 2), col = cols[1], xlab = NA, ylab = NA, xaxt = "n", yaxt = "n")
points(x = y_coords, y = x_coords, cex = scale * as.vector(sim_joint), pch = 1, lwd = 1, col = cols[2])
axis(1, lwd = 0, lwd.tick = 1)
axis(2, lwd = 0, lwd.tick = 1, las = 1)
mtext("i", 1, line = 1.5)
mtext("j", 2, line = 1.5)
legend("topright", legend = c("Monte Carlo", "ODE"), pt.cex = 2, col = cols, pch = c(19, 1), pt.lwd = c(0, 3), bty = "n", title = "method")
dev.off()




