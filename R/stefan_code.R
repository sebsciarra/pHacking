
determine_control_test_cor <- function(num_iterations, nobs.group, nvar, r) {




}

compute_control_test_cor <- function(dv_col, stefan_data) {

  #identify rows containing level data
  group_levels <- levels(stefan_data$group)
  first_level_rows <- which(stefan_data$group == group_levels[1])
  second_level_rows <- which(stefan_data$group == group_levels[2])

  #compute correlation
  control_test_cor <- cor(x = stefan_data[[dv_col]][first_level_rows], y = stefan_data[[dv_col]][second_level_rows])

  return(control_test_cor)
}

gen_stefan_data <- function(nobs.group, nvar, r) {

  # Generate group vector
  group <- factor(rep(1:length(nobs.group), nobs.group))
  # set up correlation matrix
  R <- matrix(rep(r, nvar**2), nrow = nvar)
  diag(R) <- rep(1, nvar)

  # transposed Cholesky decomposition of correlation matrix
  U <- t(chol(R))

  nobs = sum(nobs.group)

  # create random noise matrix
  random.normal <- matrix(stats::rnorm(nvar*nobs, mean = 0, sd = 1), nrow=nvar, ncol=nobs)

  # create raw data from matrix multiplication of U and random noise
  X <- as.data.frame(t(U %*% random.normal))
  x_df <- data.frame(X)

  # Generate data frame
  res <- cbind(group,  X)

  return(res)
}
