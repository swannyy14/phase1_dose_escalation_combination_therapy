uniform_weight <- function(x) {
  runif(1, 0, 1)
}

simulate_toxicity <- function(p_tox, weight_fun) {
  # toxicity outcome for a patient
  tox <- rbinom(1, 1, p_tox)
  # time to toxicity, relative to observation window
  weight <- uniform_weight()
  
  return(
    list(
      tox_response = tox,
      tox_time = weight
    )
  )
}

get_next_cohort <- function(
    p_tox, 
    weight_fun,
    cohort_size = c(1,2), 
    cohort_prob = c(0.5, 0.5),
    max_patients = 25
  ) {
  enroll <- sample(cohort_size, 1, prob = cohort_prob)
  if (recommended + enroll > max_patients) {
  }
    
}
