#' Generates data according to a covariance matrix.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of `x` and `y`.
#' @examples
#'
run_mult_simulations <- function(num_iterations, seed,  param_list, num_cores = 3) {

  #ensure reproducibility
  set.seed(seed)
  RNGkind("L'Ecuyer-CMRG")

  #progress bar settings
  progress_bar <- progressBar(min = 1, max = num_iterations, initial = 1, style = "ETA")

    p_values_list <- pbmclapply(X = 1:num_iterations, FUN = run_ind_simulation,
                            num_dvs = param_list$num_dvs,
                            num_ivs = param_list$num_ivs,
                            sample_size = param_list$sample_size,
                            dv_cor = param_list$dv_cor,
                            iv_cor = param_list$iv_cor,
                            sd = param_list$sd,
                            mc.cores = num_cores, mc.set.seed = TRUE, mc.substyle = progress_bar)

  #count the number of family comparisons that had at least one p value < .05
  p_value_counts <- lapply(p_values_list, function(x) length(which(x < .05)))

  false_positive_rate <- sum(p_value_counts >= 1)/num_iterations

  return(false_positive_rate)
}

run_ind_simulation <- function(num_iterations, num_dvs, num_ivs, sample_size, dv_cor, iv_cor,  sd) {

  #generate data
  scores <- generate_data(num_dvs = num_dvs, num_ivs = num_ivs,
                           sd = sd, dv_cor = dv_cor, iv_cor = iv_cor, sample_size = sample_size)

  #compute all t_tests
  p_values <- sapply(unique(scores$dv_iv), FUN = compute_t_test, scores = scores)

  return(p_values)
}

compute_t_test <- function(cond, scores) {

  control_data <- scores %>% filter(dv_iv == cond, condition == 'control') %>% pull(value)
  test_data <- scores %>% filter(dv_iv == cond, condition == 'test') %>% pull(value)

  p_value <- t.test(x = control_data, y = test_data)$p.value

  return(p_value)
}

generate_data <- function(num_dvs, num_ivs, sd = 15, dv_cor = 0.16, iv_cor = 0, sample_size) {

  covar <- create_covar(num_dvs, num_ivs, sd = sd, dv_cor, iv_cor)

  #generate scores and set appropriate column names
  data_long <- data.frame(mvrnorm(n = sample_size, mu = rep(100, times = nrow(covar)), Sigma = covar, empirical = F))
  colnames(data_long) <- colnames(covar)

  data_wide <- data_long %>%
    pivot_longer(cols = 1:ncol(data_long), names_to = 'name') %>%
    separate(col = 'name', into = c('dv', 'iv', 'condition'), sep = '_') %>%
    mutate_if(.predicate = is.character, .funs = factor)

  #combine condition info into one column
  scores <- data_wide %>%
    unite(col = 'dv_iv', 'dv':'iv', sep = '_') %>%
    mutate_if(.predicate = is.character, .funs = factor)

  return(scores)
}

create_covar <- function(num_dvs, num_ivs, sd = 15, dv_cor = .16, iv_cor = .16) {

  empty_covar <- create_empty_covar(num_dvs, num_ivs)
  filled_covar <- fill_covar(empty_covar, sd, dv_cor, iv_cor)

  return(filled_covar)
}
fill_covar <- function(empty_covar, sd, dv_cor, iv_cor) {

  #set covariances (i.e., correlations) between dependent variables (look for rownames in each columns with different dv name)
  #and independent variables (look for rownames in each columns with different iv name)
  empty_covar <- fill_dv_iv_covars(empty_covar, dv_cor, iv_cor, sd)
  diag(empty_covar) <- sd*sd

  #fill in remaining NA cells with diagonal values (i.e., because there is no actual difference between the control and test conditions,
  #the covariance between the control and test group for each IV under a given DV) is simply the variance
  empty_covar[is.na(empty_covar)] <- iv_cor*sd*sd

  return(empty_covar)
}

fill_dv_iv_covars <- function(empty_covar, dv_cor, iv_cor, sd) {

  for (col_number in 1:ncol(empty_covar)) {

    dv_name_to_avoid <- str_extract(colnames(empty_covar)[col_number], pattern = "[^_]+")
    iv_name_to_avoid <- str_extract(colnames(empty_covar)[col_number], pattern = "(?<=_)(.*?)(?=_)")

    #create index of row values to modify in specific column; find that do not contain dv_name and iv_name
    dv_row_values <-  !str_detect(string = rownames(empty_covar), pattern = dv_name_to_avoid)
    iv_row_values <-  !str_detect(string = rownames(empty_covar), pattern = iv_name_to_avoid)

    dv_row_values[col_number] <- FALSE #ensures diagonal value will not be overwritten
    iv_row_values[col_number] <- FALSE #ensures diagonal value will not be overwritten

    #fill cells with appropriate covariances
    empty_covar[dv_row_values, col_number] <- dv_cor*sd*sd
    empty_covar[iv_row_values, col_number] <- 0
  }

  return(empty_covar)
}

create_empty_covar <- function(num_dvs, num_ivs){

  #setup variables that describe structure of matrix
  dv_names <- sprintf(fmt = 'dv%d', 1:num_dvs)
  iv_names <- sprintf(fmt = 'iv%d', 1:num_ivs)
  conditions <- c('control', 'test')

  #create all possible combinations
  all_combinations <- data.frame(expand.grid(dv_names, iv_names, conditions))
  var_names <- paste(all_combinations$Var1, all_combinations$Var2, all_combinations$Var3, sep = '_')

  nrow_and_ncol <- length(var_names)
  #assemble empty covariance matrix that will be filled
  empty_covar <- matrix(nrow = nrow_and_ncol,
                        ncol = nrow_and_ncol,
                        dimnames = list(c(var_names), c(var_names)))

  return(empty_covar)
}


