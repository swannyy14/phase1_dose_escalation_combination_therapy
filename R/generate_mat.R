generate_dose_grid <- function(true_tox_mat) {
  nr <- nrow(true_tox_mat)
  nc <- ncol(true_tox_mat)
  
  grid <- matrix(1:(nr * nc), nrow = nr, ncol = nc, byrow = TRUE)
}

generate_orderings_matrix <- function(grid) {
  po_mat <- matrix(NA, nrow = 6, ncol = ncol(grid)*nrow(grid))
  nr <- nrow(grid)
  nc <- ncol(grid)
  
  # 1. Across Rows 
  po_mat[1, ] <- as.vector(t(grid))
  
  # 2. Up Columns 
  po_mat[2, ] <- as.vector(grid)
  
  # 3. Up Diagonals 
  po_mat[3, ] <- unlist(lapply(2:(nr+nc), function(s) {
    unlist(mapply(function(i, j) {
      if (i >= 1 && i <= nr && j >= 1 && j <= nc && (i + j) == s) grid[i, j]
    }, rep(1:nr, each=nc), rep(1:nc, times=nr)))
  }))
  
  # 4. Down Diagonals 
  po_mat[4, ] <- unlist(lapply((1-nc):(nr-1), function(d) {
    unlist(mapply(function(i, j) {
      if (i >= 1 && i <= nr && j >= 1 && j <= nc && (i - j) == d) grid[i, j]
    }, rep(1:nr, each=nc), rep(1:nc, times=nr)))
  }))
  
  # 5. Down-Up Diagonals
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
  
  # 6. Up-Down Diagonals
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
  
  # Remove duplicate ordering
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
