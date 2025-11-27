library(nnet)

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

# enroll next patient based on recommended dose
enroll_patient <- function(
    patient_data,
    enrollment_date,
    assigned_dose,
    true_tox_vec,
    n = NULL,
    verbose = FALSE
) {
  n_enrolled <- patient_data$n_enrolled
  
  if (is.null(n)) {
    if (n_enrolled == (nrow(patient_data$data) - 1)) {
      n <- 1
    } else {
      n <- sample(c(1,2), prob = c(0.8, 0.2), 1)
    }
  }
  
  for (i in 1:n) {
    if (verbose) 
      print(paste0(
        "Enrolling new patient: date - ", enrollment_date,
        ", dose - ", assigned_dose,
        ", toxic prob - ", true_tox_vec[assigned_dose])
      )
    patient_data$data[n_enrolled + i, "enrollment_date"] <- enrollment_date
    patient_data$data[n_enrolled + i, "dose_assignment"] <- assigned_dose
    tox <- rbinom(1, 1, true_tox_vec[assigned_dose])
    patient_data$data[n_enrolled + i,"tox"] <- tox
    if (tox == 1) {
      patient_data$data[n_enrolled + i,"time_to_tox"] <- sample(56, 1)
    } else {
      patient_data$data[n_enrolled + i,"time_to_tox"] <- 0
    }
  }
  
  patient_data$n_enrolled <- n + patient_data$n_enrolled
  
  return(patient_data)
}

get_next_dose_po_tite_crm <- function(
    patient_data,
    current_date,
    po_mat,
    skeleton_mat,
    target_pt = 0.25
) {
  n_enrolled <- patient_data$n_enrolled
  M <- nrow(po_mat)
  
  # get the response and weights as vectors
  response <- patient_data$data$tox[1:n_enrolled]
  tox_time <- patient_data$data$time_to_tox[1:n_enrolled]
  enrollment_dates <- patient_data$data$enrollment_date[1:n_enrolled]
  
  # find patients with censored followup up to the next accrual period
  i <- (tox_time > 30 | response == 0) & ((current_date - enrollment_dates) < 56)
  
  # find which patients are "censored" and assign weights accordingly
  if (sum(i) >= 1) {
    response[i] <- 0
    weights <- rep(1, length(response))
    weights[i] <- 30 / 56
  } else {
    weights <- rep(1, length(response))
  }
  
  # calculate mle for each partial ordering
  # and identify "best" partial order
  dose_assignment <- patient_data$data$dose_assignment[1:n_enrolled]
  
  a_mle <- numeric(M)
  log_lik <- numeric(M)
  
  for (m in 1:M) {
    # get the prior dlt probability corresponding to the dose and partial order for each observation
    partial_order <- po_mat[m,]
    alphas <- skeleton_mat[m,match(dose_assignment, partial_order)]
    
    # get mle and associated log likelihood
    optimized <- optimize_log_lik(response = response, alph = alphas, weights = weights)
    a_mle[m] <- optimized$maximum
    log_lik[m] <- optimized$objective
  }
  
  # find the partial order with maximum log likelihood
  # this assigns uniform prior for each order
  chosen_po <- which.is.max(log_lik)
  
  # get next recommended dose based on the partial ordering, 
  # skeleton, and maximum likelihood estimate
  po <- po_mat[chosen_po,]
  est_p <- skeleton_mat[chosen_po,]^a_mle[chosen_po]
  
  # return maximum dose under target P(DLT)
  if (est_p[1] <= target_pt) {
    next_dose <- po[max(which(est_p <= target_pt))]
  } else {
    next_dose <- 1
  }
  
  
  return(list(
    next_dose = next_dose,
    partial_order = chosen_po,
    a_mle = a_mle[chosen_po]
  ))
}

begin_sim <- function(
    true_tox_vec,
    po_mat,
    skeleton_mat,
    target_pt = 0.25,
    max_patients = 25,
    verbose = FALSE
) {
  
  # initiate patient data
  patient_data <- list(
    data = data.frame(
      patient_id = seq(1, max_patients),
      enrollment_date = NA,
      dose_assignment = NA,
      tox = NA,
      time_to_tox = NA
    ),
    n_enrolled = 0
  )
  
  # Initial phase should be algorithmic bc mle doesn't exist
  is_step1_phase <- TRUE
  current_date <- 0
  
  # begin with first dose
  zone <- 5 # for initial dose escalation
  treatment_esc_order <- c(
    1,2,4,3,5,7,6,8,9
  )
  next_dose <- treatment_esc_order[zone]
  
  if (verbose) print("Begin step 1...")
  while (is_step1_phase & patient_data$n_enrolled < max_patients) {
    
    # enroll 2 patients
    patient_data <- enroll_patient(
      patient_data, 
      enrollment_date = current_date, 
      assigned_dose = next_dose,
      true_tox_vec = true_tox_vec,
      n = 2,
      verbose = verbose
    )
    
    # set current date to next recruitment date
    current_date <- current_date + 30
    
    # get the uncensored response
    response <- patient_data$data$tox[1:patient_data$n_enrolled]
    tox_time <- patient_data$data$time_to_tox[1:patient_data$n_enrolled]
    enrollment_dates <- patient_data$data$enrollment_date[1:patient_data$n_enrolled]
    
    # find which patients are "censored" and assign weights accordingly
    censored <- enrollment_dates + tox_time > current_date
    if (any(censored)) {
      response[censored] <- 0
    }
    
    # decide if we should switch to crm phase
    if (length(unique(response)) > 1) {
      # switch if there are both non dlt and dlt
      if (verbose) {
        print("Both toxicity and non-toxicity observed")
      }
      is_step1_phase <- FALSE
    } else if (all(response == 1)) {
      # if all toxic response
      
      if (zone > 1) { 
        # if still not at lowest dose, then de-escalate
        zone <- zone - 1
        next_dose <- treatment_esc_order[zone]
      } else { 
        # if at lowest dose, stop trial
        if (verbose) print("Both toxicity observed at beginning")
        
        sim_result <- list(
          data = patient_data$data,
          n_enrolled = patient_data$n_enrolled,
          MTD = 1
        )
        
        return(sim_result)
      }
      
    } else {
      # if all no response
      if (zone < 9) {
        zone <- zone + 1
        next_dose <- treatment_esc_order[zone]
      } else {
        # it could be that there's delayed DLT for dose 9
        if (!all(tox_time == 0)) {
          # assume we wait 1 more month to wait for DLT to show up
          current_date <- 30
          is_step1_phase <- FALSE
        } else {
          # stop trial, return patients data
          sim_result <- list(
            data = patient_data$data,
            n_enrolled = patient_data$n_enrolled,
            MTD = 9
          )
          return(sim_result)
        }
        
      }
    }
  }
  
  # move onto stage 2
  if (verbose) print("Moving on to step 2...")
  
  # get next dose recommendation based on accumulated data
  po_tite_crm_out <- get_next_dose_po_tite_crm(
    patient_data = patient_data,
    current_date = current_date,
    po_mat = po_mat,
    skeleton_mat = skeleton_mat,
    target_pt = target_pt
  )
  
  next_dose <- po_tite_crm_out$next_dose
  
  # continue dose finding until 25 patients are enrolled
  while (patient_data$n_enrolled < max_patients) {
    # enroll the next set of patients based on next recomended dose
    patient_data <- enroll_patient(
      patient_data, 
      enrollment_date = current_date, 
      assigned_dose = next_dose,
      true_tox_vec = true_tox_vec,
      verbose = verbose
    )
    
    # increment current date
    current_date <- current_date + 30
    
    if (patient_data$n_enrolled > max_patients) {
      stop("n_enrolled cannot be higher than max_patients")
    }
    
    if (patient_data$n_enrolled != max_patients) {
      po_tite_crm_out <- get_next_dose_po_tite_crm(
        patient_data = patient_data,
        current_date = current_date,
        po_mat = po_mat,
        skeleton_mat = skeleton_mat
      )
    } else {
      # get MTD
      po_tite_crm_out <- get_next_dose_po_tite_crm(
        patient_data = patient_data,
        current_date = current_date + 30,
        po_mat = po_mat,
        skeleton_mat = skeleton_mat
      )
    }
    
    next_dose <- po_tite_crm_out$next_dose
    
    # Early stopping criteria 1: if the next recommended dose already has 
    # more than 10 patients assigned, stop trying and recommend that dose
    if (sum(patient_data$data$dose_assignment == next_dose, na.rm=TRUE) > 10) {
      sim_result <- list(
        data = patient_data$data,
        n_enrolled = patient_data$n_enrolled,
        MTD = next_dose
      )
      return(sim_result)
    }
    
    # Early stopping criteria 2: if there is more than 2 patient assigned to
    # the lowest dose and P(DLT) for lowest dose is > 0.33
    if (next_dose == 1) {
      p_dlt_lowest_dose <- skeleton_mat[po_tite_crm_out$partial_order,1]^po_tite_crm_out$a_mle
      if (p_dlt_lowest_dose > 0.33 & sum(patient_data$data$dose_assignment == 1, na.rm = TRUE) >= 3) {
        sim_result <- list(
          data = patient_data$data,
          n_enrolled = patient_data$n_enrolled,
          MTD = 1
        )
        return(sim_result)
      }
    }
  }
  
  sim_result <- list(
    data = patient_data$data,
    n_enrolled = patient_data$n_enrolled,
    MTD = next_dose
  )
  
  return(sim_result)
  
}

simulate_trials <- function(
    true_tox_vec,
    po_mat,
    skeleton_mat,
    target_pt = 0.25,
    max_patients = 25,
    N = 100,
    verbose = FALSE
) {
  operating_char <- list(
    MTD = integer(N),
    n_enrolled = integer(N),
    adverse_event_rate = numeric(N)
  )
  
  for (i in 1:N) {
    sim_data <- begin_sim(
      true_tox_vec = true_tox_vec,
      po_mat = po_mat,
      skeleton_mat = skeleton_mat,
      target_pt = target_pt,
      max_patients = max_patients,
      verbose = verbose
    )
    
    operating_char[["MTD"]][i] <- sim_data$MTD
    operating_char[["n_enrolled"]][i] <- sim_data$n_enrolled
    operating_char[["adverse_event_rate"]][i] <- mean(sim_data$data$tox)
  }
  
  return(operating_char)
}
