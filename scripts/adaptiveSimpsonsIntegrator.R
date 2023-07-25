adaptiveSimpsonsIntegrator <- function(f, lower, upper, tol = 1e-5, num_points_init = 20, ...) {
 
  # create initial points
  points <- seq(lower, upper, length.out = num_points_init)
  values <- sapply(points, f, ...)
  
  # repeatedly subdivide
  num_eval <- 5
  repeat {
    
    # otherwise, determine where the new node goes
    new_node_list <- getNewNode(points, values, tol)
    
    # terminate if we're good
    if ( new_node_list$local_error < tol ) {
    # if ( new_node_list$global_error < tol ) {
      break
    }

    # retrieve the new node
    new_node <- new_node_list$node
    cat("Adding a new node at: ", new_node, ", with error: ", new_node_list$local_error, " -- ", new_node_list$global_error, "\n", sep ="")
    
    # otherwise, compute the function for the new node
    new_val <- f(new_node, ...)
    
    # add the new node to the vectors
    points <- c(points, new_node)
    values <- c(values, new_val)
    
    # sort
    values <- values[order(points)]
    points <- points[order(points)]
   
    # increment number of evaluations
    num_eval <- num_eval + 1
     
  }
  
  # compute the integral
  integral <- computeIntegral(points, values)
  
  # return list
  res <- list(
    integral  = integral,
    num_eval  = num_eval,
    points    = points,
    values    = values,
    integrand = f
  )
  
}

computeIntegral <- function(points, values) {

  # get integrals per interval
  simp_estimates <- getSimpsonsIntervals(points, values)
  
  # compute sum
  integral <- sum(simp_estimates)
  
  return(integral)
  
}

posteriorMean <- function(integrator) {

  # get new coordinates
  xs <- integrator$points
  ys <- integrator$values
  fs <- ys * xs

  # get the simpsons integrals per interval
  simp_ints <- getSimpsonsIntervals(xs, fs)

  # compute the posterior mean
  posterior_mean <- sum(simp_ints) / integrator$integral

  return(posterior_mean)
  
}

posteriorVariance <- function(integrator, mean) {
  
  # get new coordinates
  xs <- integrator$points
  ys <- integrator$values
  fs <- ys * ((xs - mean)^2)
  
  # get the simpsons integrals per interval
  simp_ints <- getSimpsonsIntervals(xs, fs)
  
  # compute the posterior variance
  posterior_variance <- sum(simp_ints) / integrator$integral
  
  return(posterior_variance)
  
}

posteriorQuantile <- function(integrator, x, ...) {
  
  # get coordinates
  xs <- integrator$points
  ys <- integrator$values
  
  # truncate
  ys <- ys[xs < x]
  xs <- xs[xs < x]
  
  # add last coordinate
  xs[length(xs) + 1] <- x
  ys[length(ys) + 1] <- integrator$integrand(x, ...)
  
  # get the simpsons integrals per interval
  simp_ints <- getSimpsonsIntervals(xs, ys)
  
  # compute the posterior quantile of x
  posterior_quantile <- sum(simp_ints) / integrator$integral
  
  return(posterior_quantile)
  
}

getTrapezoidIntervals <- function(points, values) {
  
  # get the number of nodes
  num_nodes <- length(points)
  
  # get trapezoids
  trap_estimates <- numeric(num_nodes - 1)
  for(i in 2:num_nodes) {
    
    # get coordinates
    x0 <- points[i - 1]
    x1 <- points[i]
    y0 <- values[i - 1]
    y1 <- values[i]
    
    # get width
    w  <- x1 - x0
    
    # get trapezoid area
    trap_estimates[i - 1] <- w * (y0 + y1) / 2
    
  }
  
  return(trap_estimates)
  
}

getSimpsonsIntervals <- function(points, values) {
  
  # get the number of nodes
  num_nodes <- length(points)
  
  # use parabola
  simp_nodes <- seq(2, num_nodes - 1, 1)
  simp_estimates <- numeric(num_nodes - 1)
  for(i in 1:length(simp_nodes)) {
    
    if (i == 1) { # left boundary
      
      # only need one parabola
      pt_idx <- 1:3 + simp_nodes[i] - 2
      xs <- points[pt_idx]
      ys <- matrix(values[pt_idx])
      
      # solve for the coefficients
      X <- matrix(c(xs^2, xs, xs^0), ncol = 3)
      coeffs <- (solve(X) %*% ys)[,1]
      
      # integrate from xs[1] to xs[2]
      simp_estimates[1] <- integrateParabola(coeffs, xs[1:2])
      
      # integrate from x[2] to xs[3]
      simp_estimates[2] <- integrateParabola(coeffs, xs[2:3])
      
    } else if (i == length(simp_nodes)) { # right boundary
      
      # get other parabola
      pt_idx <- 1:3 + simp_nodes[i] - 2
      xs <- points[pt_idx]
      ys <- matrix(values[pt_idx])
      
      # solve for the coefficients
      X <- matrix(c(xs^2, xs, xs^0), ncol = 3)
      coeffs <- (solve(X) %*% ys)[,1]
      
      # integrate from xs[1] to xs[2]
      # averaging with first estimate
      simp_estimates[i] <- (integrateParabola(coeffs, xs[1:2]) + simp_estimates[i]) / 2
      simp_estimates[i + 1] <- integrateParabola(coeffs, xs[2:3])
      
    } else { # interior, take average
      
      # get final parabola
      pt_idx <- 1:3 + simp_nodes[i] - 2
      xs <- points[pt_idx]
      ys <- matrix(values[pt_idx])
      
      # solve for the coefficients
      X <- matrix(c(xs^2, xs, xs^0), ncol = 3)
      coeffs <- (solve(X) %*% ys)[,1]
      
      # integrate from xs[1] to xs[2]
      simp_estimates[i] <- (integrateParabola(coeffs, xs[1:2]) + simp_estimates[i]) / 2
      
      # integrate from x[2] to xs[3]
      simp_estimates[i + 1] <- integrateParabola(coeffs, xs[2:3])
      
    }
    
  }
  
  return(simp_estimates)
  
}

getNewNode <- function(points, values, tol) {

  # get the number of nodes
  num_nodes <- length(points)

  # use parabola
  simp_estimates <- getSimpsonsIntervals(points, values)
  
  # use trapezoids
  trap_estimates <- getTrapezoidIntervals(points, values)

  # compute errors per interval
  # error_per_interval <- abs(1 - simp_estimates / trap_estimates)
  error_per_interval <- abs(simp_estimates - trap_estimates)
  
  # determine the worst error  
  worst_interval <- which.max(error_per_interval)
  worst_error    <- error_per_interval[worst_interval]
  
  # decide the new node
  new_node <- (points[worst_interval] + points[worst_interval + 1]) / 2

  # compute global error
  # global_error <- abs(1 - sum(simp_estimates) / sum(trap_estimates))
  global_error <- abs(sum(simp_estimates) - sum(trap_estimates))
  # global_error <- sum(error_per_interval)
  # global_error <- sum(abs(1 - simp_estimates / trap_estimates))
  # global_error <- mean(abs(1 - simp_estimates / trap_estimates))
  
  return( list(node = new_node, local_error = worst_error, global_error = global_error) )
  
}

integrateParabola <- function(coeffs, bounds) {
  diff(((coeffs[1] * bounds^3) / 3) + ((coeffs[2] * bounds^2) / 2) + coeffs[3] * bounds)
}

