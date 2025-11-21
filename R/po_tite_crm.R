
calc_log_lik <- function(a, response, alph, weights) {
  # calculate likelihood
  
  # a: parameter to be estimated
  # response: a vector of response for enrolled patients
  # alph: vector of skeleton probability corresponding to dose assigned to each patient
  # weights: weights associated with each patient
  
  sum(response*a*log(alph) +
        (1-response)*log(1 - weights*(alph)^a))
}

optimize_log_lik <- function(response, alph, weights) {
  # find a that maximizes the likelihood
  
  a_mle <- optimize(
    f = calc_log_lik, interval = c(0, 100),
    response = response,
    alph = alph,
    weights = weights,
    maximum = TRUE
  )
  
  return(a_mle)
}

calc_weight <- function(current_time, followup_time) {
  return(min(current_time, followup_time))
}
