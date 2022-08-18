false_pos_exist <- function(iteration_num, sim_data){

  false_positive_exist <- ifelse(test = sum(sim_data$p_values[[iteration_num]] < .05) > 0, yes = 1, no = 0)

  return(false_positive_exist)
}
