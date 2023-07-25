source("../scripts/countMatrixClass2.R", chdir = TRUE)

# compute the observation times per size
# goal: see how many days it takes until the probability of 
# more than 1000 extant lineages is greater than 5%

# diversification parameters
r0      <- 2.5
phi     <- 1 / 7
lambda0 <- r0 * phi
factor  <- c(1.25, 1.5, 1.75, 2)

# mutation rate parameters
genome_size <- 30000
mutation_rate_per_site_per_year <- 0.00084
mutation_rate_per_site_per_day  <- mutation_rate_per_site_per_year / 365
mu <- mutation_rate_per_site_per_day

# days
days <- 1:60

# desired number of lineages
num_lineages <- 1000

# for each factor, do the integration
end_times <- numeric(length(factor))
for(i in 1:length(factor)) {
  
  # get the factor
  this_factor <- factor[i]
  
  # compute the rates
  growth_rate <- lambda0 * this_factor - phi
  
  # compute the geometric parameter
  p <- exp(-1.0 * growth_rate * days)
  
  # compute the number per day
  num_per_day <- 1 + qgeom(0.95, prob = p)
  
  # compute the last day less than the desired number
  this_day <- max(days[num_per_day < num_lineages])
  
  end_times[i] <- this_day
  
}