generate_dose_grid <- function(true_tox_mat) {
  nr <- nrow(true_tox_mat)
  nc <- ncol(true_tox_mat)
  grid <- matrix(1:(nr * nc), nrow = nr, ncol = nc, byrow = TRUE)
}

generate_orderings_matrix <- function(grid) {
  po_mat <- matrix(NA, nrow = 6, ncol = ncol(grid)*nrow(grid))
  nr <- nrow(grid)
  nc <- ncol(grid)
  po_mat[1, ] <- as.vector(t(grid))
  po_mat[2, ] <- as.vector(grid)
  po_mat[3, ] <- unlist(lapply(2:(nr+nc), function(s) {
    unlist(mapply(function(i, j) {
      if (i >= 1 && i <= nr && j >= 1 && j <= nc && (i + j) == s) grid[i, j]
    }, rep(1:nr, each=nc), rep(1:nc, times=nr)))
  }))
  po_mat[4, ] <- unlist(lapply((1-nc):(nr-1), function(d) {
    unlist(mapply(function(i, j) {
      if (i >= 1 && i <= nr && j >= 1 && j <= nc && (i - j) == d) grid[i, j]
    }, rep(1:nr, each=nc), rep(1:nc, times=nr)))
  }))
  po_mat[5, ] <- {
    result <- c()
    for (k in 1:(nr + nc - 1)) {
      temp <- c()
      for (i in 1:nr) {
        j <- k - i + 1
        if (j >= 1 && j <= nc) {
          temp <- c(temp, grid[i, j])
        }
      }
      result <- c(result, temp)
    }
    result
  }
  po_mat[6, ] <- {
    result <- c()
    for (k in 1:(nr + nc - 1)) {
      temp <- c()
      for (i in nr:1) {
        j <- k - i + 1
        if (j >= 1 && j <= nc) {
          temp <- c(temp, grid[i, j])
        }
      }
      result <- c(result, temp)
    }
    result
  }
  po_mat <- unique(po_mat)
  return(po_mat)
}

generate_skeleton_mat <- function(po_mat, skeleton_vec) {
  d <- length(skeleton_vec)
  s <- nrow(po_mat)
  skeleton_mat <- matrix(0, nrow = s, ncol = d)
  for (j in 1:s) {
    skeleton_mat[j, ] <- skeleton_vec[order(po_mat[j, ])]
  }
  return(skeleton_mat)
}

library(nnet)

calc_log_lik <- function(a, response, alph, weights) {
  sum(response*a*log(alph) +
        (1-response)*log(1 - weights*(alph)^a))
}

optimize_log_lik <- function(response, alph, weights) {
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
      n <- sample(c(1,2), 1)
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
  response <- patient_data$data$tox[1:n_enrolled]
  tox_time <- patient_data$data$time_to_tox[1:n_enrolled]
  enrollment_dates <- patient_data$data$enrollment_date[1:n_enrolled]
  i <- (tox_time > 30 | response == 0) & ((current_date - enrollment_dates) < 56)
  if (sum(i) >= 1) {
    response[i] <- 0
    weights <- rep(1, length(response))
    weights[i] <- 30 / 56
  } else {
    weights <- rep(1, length(response))
  }
  dose_assignment <- patient_data$data$dose_assignment[1:n_enrolled]
  a_mle <- numeric(M)
  log_lik <- numeric(M)
  for (m in 1:M) {
    partial_order <- po_mat[m,]
    alphas <- skeleton_mat[m,match(dose_assignment, partial_order)]
    optimized <- optimize_log_lik(response = response, alph = alphas, weights = weights)
    a_mle[m] <- optimized$maximum
    log_lik[m] <- optimized$objective
  }
  chosen_po <- which.is.max(log_lik)
  po <- po_mat[chosen_po,]
  est_p <- skeleton_mat[chosen_po,]^a_mle[chosen_po]
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
  is_step1_phase <- TRUE
  current_date <- 0
  next_dose <- 1
  zone <- 1
  treatment_esc_order <- c(
    1,2,4,3,5,7,6,8,9
  )
  if (verbose) print("Begin step 1...")
  while (is_step1_phase & patient_data$n_enrolled < max_patients) {
    patient_data <- enroll_patient(
      patient_data,
      enrollment_date = current_date,
      assigned_dose = next_dose,
      true_tox_vec = true_tox_vec,
      n = 2,
      verbose = verbose
    )
    current_date <- current_date + 30
    response <- patient_data$data$tox[1:patient_data$n_enrolled]
    tox_time <- patient_data$data$time_to_tox[1:patient_data$n_enrolled]
    enrollment_dates <- patient_data$data$enrollment_date[1:patient_data$n_enrolled]
    censored <- enrollment_dates + tox_time > current_date
    if (any(censored)) {
      response[censored] <- 0
    }
    if (length(unique(response)) > 1) {
      if (verbose) {
        print("Both toxicity and non-toxicity observed")
      }
      is_step1_phase <- FALSE
    } else if (all(response == 1)) {
      if (verbose) print("Both toxicity observed at beginning")
      sim_result <- list(
        data = patient_data$data,
        n_enrolled = patient_data$n_enrolled,
        MTD = 1
      )
      return(sim_result)
    } else {
      if (zone < 9) {
        zone <- zone + 1
        next_dose <- treatment_esc_order[[zone]][sample(length(treatment_esc_order[[zone]]), 1)]
      } else {
        if (!all(tox_time == 0)) {
          current_date <- 30
          is_step1_phase <- FALSE
        } else {
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
  if (verbose) print("Moving on to step 2...")
  po_tite_crm_out <- get_next_dose_po_tite_crm(
    patient_data = patient_data,
    current_date = current_date,
    po_mat = po_mat,
    skeleton_mat = skeleton_mat,
    target_pt = target_pt
  )
  next_dose <- po_tite_crm_out$next_dose
  while (patient_data$n_enrolled < max_patients) {
    patient_data <- enroll_patient(
      patient_data,
      enrollment_date = current_date,
      assigned_dose = next_dose,
      true_tox_vec = true_tox_vec,
      verbose = verbose
    )
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
      po_tite_crm_out <- get_next_dose_po_tite_crm(
        patient_data = patient_data,
        current_date = current_date + 30,
        po_mat = po_mat,
        skeleton_mat = skeleton_mat
      )
    }
    next_dose <- po_tite_crm_out$next_dose
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
      target_pt = 0.25,
      verbose = verbose
    )
    operating_char[["MTD"]][i] <- sim_data$MTD
    operating_char[["n_enrolled"]][i] <- sim_data$n_enrolled
    operating_char[["adverse_event_rate"]][i] <- mean(sim_data$data$tox, na.rm = TRUE)

  }
  return(operating_char)
}

dummy_mat <- matrix(0, nrow = 3, ncol = 3)
grid <- generate_dose_grid(dummy_mat)
po_mat <- generate_orderings_matrix(grid)

skeleton_vec <- c(
  0.05, 0.08, 0.10,
  0.13, 0.15, 0.18,
  0.22, 0.25, 0.35
)

skeleton_mat <- generate_skeleton_mat(po_mat, skeleton_vec)

scen1 <- c(
  0.05, 0.08, 0.10,
  0.13, 0.15, 0.18,
  0.20, 0.25, 0.35
)

scen2 <- c(
  0.20, 0.30, 0.45,
  0.25, 0.30, 0.45,
  0.30, 0.45, 0.50
)

scen3 <- c(
  0.05, 0.06, 0.08,
  0.09, 0.11, 0.14,
  0.14, 0.20, 0.35
)

scen4 <- c(
  0.05, 0.12, 0.22,
  0.06, 0.14, 0.28,
  0.08, 0.16, 0.35
)

set.seed(123)
oc1 <- simulate_trials(scen1, po_mat, skeleton_mat, target_pt = 0.25, N = 500)
oc2 <- simulate_trials(scen2, po_mat, skeleton_mat, target_pt = 0.25, N = 500)
oc3 <- simulate_trials(scen3, po_mat, skeleton_mat, target_pt = 0.25, N = 500)
oc4 <- simulate_trials(scen4, po_mat, skeleton_mat, target_pt = 0.25, N = 500)

table(oc1$MTD) / length(oc1$MTD)
mean(oc1$n_enrolled)
mean(oc1$adverse_event_rate)

table(oc2$MTD) / length(oc2$MTD)
table(oc3$MTD) / length(oc3$MTD)
table(oc4$MTD) / length(oc4$MTD)

