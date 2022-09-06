#' Simulations optional stopping.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of `x` and `y`.
#' @export
run_exp_simulations <- function(sim_params, num_iterations, num_cores, seed, step_size = T){

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  cells <- determine_cells(sim_params, step_size = step_size)

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
                              dv_cor = cells$dv_cor[cell],
                              iv_cor = cells$iv_cor[cell],
                              sd = 15,
                              mc.cores = num_cores, mc.set.seed = TRUE, mc.substyle = progress_bar)

    sim_results_list[[cell]] <- rbindlist(l = results, use.names = T, fill = T)
  }


  sim_results <- rbindlist(l = sim_results_list, use.names = T)

  return(sim_results)
}

determine_cells <- function(sim_params, step_size) {

  #compute all possible experimental cells
  cells <- expand.grid(sim_params)

  #remove rows where n_min > n_max  (impossible conditions)
  cells <- cells[!cells$n_min > cells$n_max, ]

  #select rows (i.e., cells) where n_max = n_min
  if (step_size == F) {
    cells <- cells [cells$n_min == cells$n_max, ]
  }

  return(cells)
}

optional_stopping <- function(num_iterations, step_size, n_min, n_max, num_dvs, num_ivs, dv_cor, iv_cor, sd) {

  if (n_min > n_max) {
    stop('n_min must be greater than n_max')
  }

  # generate entire data set
  entire_data <- generate_data(num_dvs = num_dvs, num_ivs = num_ivs, sd = sd, dv_cor = dv_cor, iv_cor = iv_cor, sample_size = n_max)

  #compute sample_size increments if n_max > n_min
  sample_size_increments <- increments_or_not(n_max = n_max, n_min = n_min, step_size = step_size)

  #extract data from entire_data according to sample_size increments. Store each value in a list.
  step_size_data_list <- unlist(lapply(X = 1:length(sample_size_increments), FUN = extract_step_size_data,
                                entire_data = entire_data, sample_size_increments = sample_size_increments), recursive = F)

  #compute all possible t-tests/regressions on each step size data set
  step_size_p_values <- unlist(lapply(X = step_size_data_list, FUN = compute_all_t_tests))

  return(data.table('step_size' = step_size,
                    'n_min' = n_min,
                    'n_max' = n_max,
                    'num_dvs' = num_dvs,
                    'num_ivs' = num_ivs,
                    'dv_cor' = dv_cor,
                    'iv_cor' = iv_cor,
                    'p_values' = list(step_size_p_values)))
}

compute_step_size_t_tests <- function(data_list_elements_to_merge, step_size_data_list) {

  analytical_data <- rbindlist(l = step_size_data_list[1:data_list_elements_to_merge])

  p_values <- compute_all_t_tests(scores = analytical_data)

  return(p_values)
}

extract_step_size_data <- function(sample_size_increment_ind, entire_data, sample_size_increments){

  #find sum of sample sizes up until sample_size_increment_ind - 1 (tells us what data points have been included)
  sample_size_sum <- sum(sample_size_increments[1:sample_size_increment_ind])

  #extract sample_size_increment data points from control and test groups of each dv_iv level
 step_size_data <- lapply(X = levels(entire_data$dv_iv), FUN = extract_control_test_data,
                          entire_data = entire_data, sample_size_sum = sample_size_sum)

 return(step_size_data)
}

extract_control_test_data <- function(dv_iv_level, entire_data, sample_size_sum){

  #extract control and test data (note how a rolling sum is recorded)
  control_datum <- entire_data[entire_data$dv_iv == dv_iv_level & entire_data$condition == 'control', ][1:sample_size_sum, ]
  test_datum <- entire_data[entire_data$dv_iv == dv_iv_level & entire_data$condition == 'test', ][1:sample_size_sum, ]

  control_test_data <- rbind(control_datum, test_datum)

  return(control_test_data)
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

  #return rolling sum
  return(sample_size_increments)
}

compute_all_t_tests <- function(step_size_data) {

  #compute t-tests
  p_values <- sapply(X = unique(step_size_data$dv_iv), FUN = compute_t_test, scores = step_size_data)

  return(p_values)
}


generate_step_size_data <- function(increment, num_dvs, num_ivs, dv_cor, iv_cor) {

  #if increment == 1, then generate 2 data points and extract on control and test point for each dv_iv level
  if (increment == 1) {

    new_data <- generate_data(num_dvs = num_dvs, num_ivs = num_ivs, sample_size = 2,
                              dv_cor = dv_cor, iv_cor = iv_cor)

    #extract one control and test data point for each dv_iv  level
    new_data <- rbindlist(lapply(X = levels(new_data$dv_iv), FUN = extract_control_test_data, new_data = new_data))
  }

  else {
    new_data <- generate_data(num_dvs = num_dvs, num_ivs = num_ivs, sample_size = increment,
                              dv_cor = dv_cor, iv_cor = iv_cor)
  }


  return(new_data)
}




