library(tidyverse)
library(nnet)
library(knitr)
library(kableExtra)
library(dfcrm)

set.seed(123)

source("R/po_tite_crm.R")
source("R/generate_mat.R")

po_mat <- matrix(NA, nrow = 6, ncol = 9)
po_mat[1, ] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
po_mat[2, ] <- c(1, 4, 7, 2, 5, 8, 3, 6, 9)
po_mat[3, ] <- c(1, 2, 4, 3, 5, 7, 6, 8, 9)
po_mat[4, ] <- c(1, 4, 2, 7, 5, 3, 8, 6, 9)
po_mat[5, ] <- c(1, 2, 4, 7, 5, 3, 6, 8, 9)
po_mat[6, ] <- c(1, 4, 2, 3, 5, 7, 8, 6, 9)

skeleton_mat <- matrix(
  rep(c(0.05, 0.08, 0.10, 0.13, 0.15, 0.18, 0.22, 0.25, 0.35),
      times = nrow(po_mat)),
  nrow = nrow(po_mat), ncol = ncol(po_mat), byrow = TRUE
)

target_pt <- 0.25
max_patients <- 25

simulate_trials <- function(
    true_tox_vec,
    po_mat,
    skeleton_mat,
    target_pt = 0.25,
    max_patients = 25,
    N = 1000,
    verbose = FALSE,
    include_patient_data = FALSE
) {
  operating_char <- list(
    MTD = integer(N),
    n_enrolled = integer(N),
    adverse_event_rate = numeric(N),
    max_day = integer(N),
    n_assigned_by_dose = integer(length(true_tox_vec)),
    n_dlt_by_dose = integer(length(true_tox_vec)),
    n_sim = N
  )

  if (include_patient_data) {
    operating_char$patient_data <- list()
  }

  for (i in 1:N) {
    sim_data <- begin_sim(
      true_tox_vec = true_tox_vec,
      po_mat = po_mat,
      skeleton_mat = skeleton_mat,
      target_pt = target_pt,
      max_patients = max_patients,
      verbose = verbose
    )

    if (include_patient_data) {
      operating_char$patient_data[[i]] <- sim_data$data
    }

    operating_char$MTD[i] <- sim_data$MTD
    operating_char$n_enrolled[i] <- sim_data$n_enrolled
    operating_char$adverse_event_rate[i] <- mean(sim_data$data$tox, na.rm = TRUE)
    operating_char$max_day[i] <- max(sim_data$data$enrollment_date, na.rm = TRUE)

    dose_data <- sim_data$data %>%
      filter(!is.na(dose_assignment)) %>%
      group_by(dose_assignment) %>%
      summarise(
        n = n(),
        n_dlt = sum(tox),
        .groups = "drop"
      )

    operating_char$n_assigned_by_dose[dose_data$dose_assignment] <-
      operating_char$n_assigned_by_dose[dose_data$dose_assignment] + dose_data$n

    operating_char$n_dlt_by_dose[dose_data$dose_assignment] <-
      operating_char$n_dlt_by_dose[dose_data$dose_assignment] + dose_data$n_dlt
  }

  operating_char
}

clean_results <- function(op_char, true_tox_vec) {
  op_table <- tibble(
    Treatment = names(true_tox_vec),
    true_p_dlt = true_tox_vec
  ) %>%
    mutate(
      p_select_mtd = as.integer(table(factor(op_char$MTD, levels = 1:9))) / length(op_char$MTD),
      avg_num_assigned = op_char$n_assigned_by_dose / op_char$n_sim,
      avg_num_dlt = op_char$n_dlt_by_dose / op_char$n_sim
    ) %>%
    pivot_longer(
      cols = c("true_p_dlt", "p_select_mtd", "avg_num_assigned", "avg_num_dlt"),
      names_to = "op_char",
      values_to = "value"
    ) %>%
    pivot_wider(
      id_cols = "op_char",
      names_from = "Treatment",
      values_from = "value"
    )
  
  print(paste("Average number of patients enrolled:", mean(op_char$n_enrolled)))
  print(paste("Average number of days until last enrollment:", mean(op_char$max_day)))
  
  return(op_table)
}

dose_names <- c(
  "P0.5+H2", "P0.5+H3", "P0.5+H5",
  "P0.75+H2", "P0.75+H3", "P0.75+H5",
  "P1.0+H2", "P1.0+H3", "P1.0+H5"
)

scen1 <- c(
  0.05, 0.08, 0.10,
  0.13, 0.15, 0.18,
  0.20, 0.25, 0.35
)

scen2_a <- c(
  0.03, 0.06, 0.07,
  0.09, 0.10, 0.13,
  0.15, 0.17, 0.24
)

scen2_b <- c(
  0.03, 0.09, 0.15,
  0.06, 0.10, 0.17,
  0.07, 0.13, 0.24
)

scen3_a <- c(
  0.07, 0.10, 0.13,
  0.17, 0.20, 0.23,
  0.29, 0.35, 0.45
)

scen3_b <- c(
  0.07, 0.17, 0.29,
  0.10, 0.20, 0.35,
  0.13, 0.23, 0.45
)

scen4 <- c(
  0.05, 0.06, 0.08,
  0.12, 0.14, 0.16,
  0.22, 0.28, 0.35
)

scen5 <- c(
  0.05, 0.15, 0.30,
  0.07, 0.17, 0.32,
  0.09, 0.19, 0.35
)

scen6 <- c(
  0.05, 0.07, 0.09,
  0.12, 0.15, 0.18,
  0.22, 0.25, 0.30
)

scen7 <- c(
  0.05, 0.08, 0.10,
  0.08, 0.15, 0.22,
  0.10, 0.22, 0.45
)

scen8 <- c(
  0.15, 0.20, 0.25,
  0.20, 0.25, 0.30,
  0.30, 0.30, 0.35
)

scen9 <- c(
  0.30, 0.35, 0.40,
  0.35, 0.40, 0.50,
  0.40, 0.45, 0.55
)

scenario_list <- list(
  scen1 = scen1,
  scen2_a = scen2_a,
  scen2_b = scen2_b,
  scen3_a = scen3_a,
  scen3_b = scen3_b,
  scen4 = scen4,
  scen5 = scen5,
  scen6 = scen6,
  scen7 = scen7,
  scen8 = scen8,
  scen9 = scen9
)

results_list <- list()

for (k in seq_along(scenario_list)) {
  cat("\n================ Scenario", k, "================\n")
  scen_vec <- scenario_list[[k]]

  true_tox_mat <- matrix(scen_vec, byrow = FALSE, nrow = 3, ncol = 3)
  true_tox_vec <- as.vector(true_tox_mat)
  names(true_tox_vec) <- dose_names

  print(true_tox_vec)

  oc <- simulate_trials(
    true_tox_vec = true_tox_vec,
    po_mat = po_mat,
    skeleton_mat = skeleton_mat,
    target_pt = target_pt,
    max_patients = max_patients,
    N = 1000,
    verbose = FALSE
  )

  res_tbl <- clean_results(oc, true_tox_vec)
  print(res_tbl)

  results_list[[k]] <- res_tbl
}

results_list[[5]][,-1] %>% 
  kable %>%
  kable_styling
