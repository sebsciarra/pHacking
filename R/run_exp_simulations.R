#' Simulations optional stopping.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of `x` and `y`.
run_exp_simulations <- function(sim_params, num_iterations, num_cores, seed){

  #determine all experimental conditions
  cells <- expand.grid(sim_params)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  #progress bar settings
  progress_bar <- progressBar(min = 1, initial = 0, max = num_iterations, style = "ETA")

  #list containing results from each iteration
  sim_results_list <- vector(mode = 'list', length = nrow(cells))

  for (cell in 1:nrow(cells)) {
    results <- pbmclapply(X = 1:num_iterations,
                              FUN = optional_stopping,
                              step_size = cells$step_size[cell],
                              n_min = cells$n_min[cell],
                              n_max = cells$n_max[cell],
                              num_dvs = cells$num_dvs[cell],
                              num_ivs = cells$num_ivs[cell],
                              mc.cores = num_cores, mc.set.seed = TRUE, mc.substyle = progress_bar)

    sim_results_list[[cell]] <- rbindlist(l = results, use.names = T, fill = T)
  }


  sim_results <- rbindlist(l = sim_results_list, use.names = T)

  return(sim_results)
}


optional_stopping <- function(num_iterations, step_size, n_min, n_max, num_dvs, num_ivs) {

  if (n_min > n_max) {
    stop('n_min must be greater than n_max')
  }

  #compute sample_size increments if n_max > n_min
  sample_size_increments <- increments_or_not(n_max = n_max, n_min = n_min, step_size = step_size)

  updated_data <- c()
  p_values <- c()

  for (increment in sample_size_increments) {

    new_data <- generate_data(num_dvs = num_dvs, num_ivs = num_ivs, sample_size = increment)

    #add it to current data
    updated_data <- rbind(updated_data, new_data)

    #compute all possible t-tests/regressions
    p_values <- c(p_values, compute_all_t_tests(scores = updated_data))
  }

  #if missing = T, then remove data points according to either MCAR, MAR, or MNAR

  return(data.table('step_size' = step_size,
                    'n_min' = n_min,
                    'n_max' = n_max,
                    'num_dvs' = num_dvs,
                    'num_ivs' = num_ivs,
                    'p_values' = list(p_values)))
}

increments_or_not <- function(n_max, n_min, step_size) {

  if (n_max > n_min) {
    sample_size_increments <- compute_sample_size_increments(step_size = step_size,
                                                             n_min = n_min,
                                                             n_max = n_max)
  }
  else {
    sample_size_increments <- n_max
  }

  return(sample_size_increments)
}

compute_sample_size_increments <- function(step_size, n_min, n_max) {

  #determine number of times to take a step_size
  ##compute difference between n_max and n_min
  max_min_diff <- n_max - n_min
  sample_size_increments <- c(n_min, rep(x = step_size, times =  max_min_diff/step_size))

  #add increment so that n_max is achieved
  sum_increments <- sum(sample_size_increments)
  if (sum_increments < n_max) {

    sample_size_increments <- c(sample_size_increments, n_max - sum_increments)
  }

  return(sample_size_increments)
}

compute_all_t_tests <- function(scores) {

  #compute t-tests
  p_values <- sapply(unique(scores$dv_iv), FUN = compute_t_test, scores = scores)

  return(p_values)
}




