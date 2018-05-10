u.test <- function(x, cells, aberrations, cell.dist) {
  if (!is.data.frame(x)) 
    stop("x should be a data.frame class")
  if (!is.numeric(cell.dist))
    stop("cell.dist should be a vector of columns of the cell distribution, in this form c(start, end).")
  
  xmat <- as.matrix(x)
  xmat.dist <- list(); disp.idx <- u.result <- u.decsn <- numeric();
  for (i in 1L:nrow(x)) {
    xmat.idx <- as.numeric(xmat[i, cell.dist[1L]:cell.dist[2L]][complete.cases(xmat[i, cell.dist[1L]:cell.dist[2L]])])
    
    for (j in 1L:length(xmat.idx))
      xmat.dist[[j]] <- rep(j - 1L, xmat.idx[j])
    
    disp.idx[i]  <- var(unlist(xmat.dist)) / mean(unlist(xmat.dist))
    if (aberrations[i] == 1) {
      u.result[i] <- NA; u.decsn[i] <- NA
    } else {
      u.result[i] <- (disp.idx[i] - 1L) * sqrt((cells[i] - 1L) / (2L * (1L - (1L / aberrations[i])))) 
      if (u.result[i] < - 1.96)
        u.decsn[i] <- "under"
      if (u.result[i] > 1.96)
        u.decsn[i] <- "over"
      if ((u.result[i] >= - 1.96) && (u.result[i] <= 1.96))
        u.decsn[i] <- "normal"
    }
  }
  x$disp_idx <- disp.idx; x$u_test <- u.result; x$u_decsn <- u.decsn
  return(x)
}