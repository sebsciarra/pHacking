#' Computes false positive rate
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of `x` and `y`.
#' @export
false_pos_exist <- function(iteration_num, sim_data) {

  #if any of the p_values in the list are greater than .05, return 1
  false_positive_exist <- ifelse(test = sum(sim_data$p_values[[iteration_num]] < .05) > 0, yes = 1, no = 0)

  return(false_positive_exist)
}
