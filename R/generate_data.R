#' Generates data according to a covariance matrix.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of `x` and `y`.
#' @export
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

  covar <- create_covar(num_dvs = num_dvs, num_ivs, sd = sd, dv_cor = dv_cor, iv_cor = iv_cor)

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

create_covar <- function(num_dvs, num_ivs, sd = 15, dv_cor = .16, iv_cor = 0) {

  empty_covar <- create_empty_covar(num_dvs = num_dvs, num_ivs = num_ivs)
  filled_covar <- fill_covar(empty_covar = empty_covar, sd = sd, dv_cor = dv_cor,iv_cor =  iv_cor)

  return(filled_covar)
}

fill_covar <- function(empty_covar, dv_cor, iv_cor, sd) {

  num_rows <- nrow(empty_covar)
  col_names <- colnames(empty_covar)

  #fill empty_covar by column
  covar_list <- lapply(X = 1:num_rows, FUN = fill_empty_covar_col, empty_covar = empty_covar, dv_cor = dv_cor, iv_cor = iv_cor, sd = sd)
  covar <- matrix(unlist(covar_list), nrow =num_rows, ncol = num_rows, dimnames = list(col_names, col_names))

  #fill diagonal values with sd*sd
  diag(covar) <- sd*sd

  #fill NA cells with 0
  covar[is.na(covar)] <- 0

  return(covar)
}

fill_empty_covar_col <- function(col_number, empty_covar, dv_cor, iv_cor, sd) {

  #find rowname that does not match on dv number and matches on all other content in column nam
  ##? matches zero or one occurence of regular expression
  ##(?<={pattern}) positive lookbehind assertion, it tests whether the currently matched string is preceded by a string matching {pattern}.
  dv_name_to_match <- paste(str_extract(colnames(empty_covar)[col_number], pattern = "[^*_]+"), '_', sep = '')  #extract all content before first underscore
  iv_name_to_match <- paste(str_extract(colnames(empty_covar)[col_number], pattern = "(?<=_)(.*?)(?=_)"), '_', sep = '') #extract all content between two underscores
  level_name_to_match <- ifelse(test = str_detect(string = colnames(empty_covar)[col_number], pattern = 'control'),
                                yes = 'control', no = 'test')

  #dv_cor:  match on iv and level name, mismatch on dv_name
  dv_cor_rows <- !str_detect(string = rownames(empty_covar), pattern = dv_name_to_match) &
    str_detect(string = rownames(empty_covar), pattern = paste(iv_name_to_match, level_name_to_match, sep = ''))

  #iv_cor: only set correlations for 'test' cols such that dv_name and level_name are matched on, but not iv_name
  iv_cor_rows <-  str_detect(string = rownames(empty_covar), pattern = 'test') &
    str_detect(string = rownames(empty_covar), pattern = dv_name_to_match) & #match on dv
    !str_detect(string = rownames(empty_covar), pattern = iv_name_to_match) & #mismatch on iv_name
    str_detect(string = rownames(empty_covar), pattern = level_name_to_match) #match on level name

  #in stefan code, 1 control group data is compared with several different treatment levels. Thus, correlation between all control groups
  #for all IVs has to be 1.
  ##match on dv and control, mismatch on iv
  full_cor_rows <-   str_detect(string = rownames(empty_covar), pattern = 'control') &
    str_detect(string = rownames(empty_covar), pattern = dv_name_to_match) & #match on dv
    str_detect(string = rownames(empty_covar), pattern = level_name_to_match) & #match on level name
    !str_detect(string = rownames(empty_covar), pattern = iv_name_to_match)  #mismatch on iv_name


  #fill cells with appropriate covariances
  empty_covar[dv_cor_rows, col_number] <- dv_cor*sd*sd
  empty_covar[iv_cor_rows, col_number] <- iv_cor*sd*sd
  empty_covar[full_cor_rows, col_number] <- sd*sd

  return(empty_covar[ ,col_number])
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


