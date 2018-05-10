#' Dose-Effect Curve
#' @export
#' 
#' @description
#' Generates calibration curve for the chromosomal aberration.
#' 
#' @param x output from fit.coef function;
#' @param level confidence level for fitted values generated;
#' @param xaxis x-axis tick labels;
#' @param newdat supply newdata
#' @param yci default \code{null}; length 2 numeric value of Y_L and Y_U.
#' @param eta precision for Y_L and Y_U intersecting C_U and C_L, respectively.
dose_curve <- function (x, level = 0.95, xaxis = NULL, newdat = NULL, yci = NULL, eta = .0001, ...) {
  doses <- x$preliminary[, 1]; cells <- x$preliminary[, 2]; aberrations <- x$preliminary[, 3]
  fit <- x$fitted_values$fit; mod_est <- x$summary$coefficients; mod_cor <- x$summary$correlation
  est_se <- mod_est[, 2]; vcov_mat <- mod_cor * outer(est_se, est_se)
  if (nrow(mod_est) == 3)
    R <- qchisq(level, df = 3)
  if (nrow(mod_est) == 2)
    R <- qchisq(level, df = 2)
  
  scales <- list()
  if (!is.null(xaxis) && is.numeric(xaxis)) {
    scales <- list(
     x = list(
        at = xaxis,
        labels = xaxis
        )
     ) 
  }
  
  xyplot(aberrations / cells ~ doses, xlab = "Doses", 
         ylab = "Aberration per cell", cex = 1.5, pch = 19, col = "#0080ff",
         key = list(
           corner = c(.02, .98),
           rectangles = list(col = c("#0080ff", "#ff00ff", "darkgreen"), size = 4.5, height = 0.8),
           text = list(c("Data", "Fitted", paste(level * 100, "% C.I.", sep ='')))),
         scales = scales,
         panel = function (x, y, ...) {
           if (is.null(xaxis)) {
             panel.grid(h = -1, v = -1, col.line = "cornsilk2")
           } else if (!is.null(xaxis) && is.numeric(xaxis)) {
             panel.grid(h = -1, v = 0, col.line = "cornsilk2")
             panel.abline(v = xaxis, col = "cornsilk2")
           }
           
           if (is.null(newdat)) {
             panel.xyplot(x, y, ...) 
           } else if (!is.null(newdat)) {
             if (is.null(yci)) stop('Supply the Y_U and Y_L, see Crow et. al (1959)')
             if (nrow(mod_est) == 3) {
               dfunc <- function(y, a = mod_est[2, 1], b = mod_est[3, 1], c = mod_est[1, 1]) {
                 abs((-a + sqrt((a ^ 2) + (4 * b) * (y - c))) / (2 * b))
               }
               
               fit_df <- matrix(NA, 6000, 3)
               for (i in seq(0, 6, by = 0.001)) {
                 yy <- mod_est[1, 1] + mod_est[2, 1] * i + mod_est[3, 1] * (i ^ 2)
                 y_u <- mod_est[1, 1] + mod_est[2, 1] * i + mod_est[3, 1] * (i ^ 2) + sqrt(R) * 
                   sqrt(vcov_mat[1, 1] + vcov_mat[2, 2] * (i ^ 2) + vcov_mat[3, 3] * (i ^ 4) +
                          2 * vcov_mat[1, 2] * i + 2 * vcov_mat[1, 3] * (i ^ 2) + 2 * vcov_mat[2, 3] * (i ^ 3))
                 y_l <- mod_est[1, 1] + mod_est[2, 1] * i + mod_est[3, 1] * (i ^ 2) - sqrt(R) * 
                   sqrt(vcov_mat[1, 1] + vcov_mat[2, 2] * (i ^ 2) + vcov_mat[3, 3] * (i ^ 4) +
                          2 * vcov_mat[1,2] * i + 2 * vcov_mat[1,3] * (i ^ 2) + 2 * vcov_mat[2, 3] * (i ^ 3))
                 
                 fit_df[i * 1000, ] <- c(y_l, yy, y_u)
               }
             } else if (nrow(mod_est) == 2) {
               dfunc <- function (y, a = mod_est[2, 1], c = mod_est[1, 1]) {
                 (y - c) / a
               }
               
               fit_df <- matrix(NA, 6000, 3)
               for (i in seq(0, 6, by = 0.001)) {
                 yy <- mod_est[1, 1] + mod_est[2, 1] * i
                 y_u <- mod_est[1, 1] + mod_est[2, 1] * i + sqrt(R) * 
                   sqrt(vcov_mat[1, 1] + vcov_mat[2, 2] * (i ^ 2) + 2 * vcov_mat[1, 2] * i)
                 y_l <- mod_est[1, 1] + mod_est[2, 1] * i - sqrt(R) * 
                   sqrt(vcov_mat[1, 1] + vcov_mat[2, 2] * (i ^ 2) + 2 * vcov_mat[1, 2] * i)
                 
                 fit_df[i * 1000, ] <- c(y_l, yy, y_u)
               }
             }
             
             
             cil <- fit_df[, 3] - yci[1]
             cil_indx <- which((cil > 0) & (cil < eta))
             cil <- cil[(cil > 0) & (cil < eta)][1]
             
             ciu <- fit_df[, 1] - yci[2]
             ciu_indx <- which((ciu > 0) & (ciu < eta))
             ciu <- ciu[(ciu > 0) & (ciu < eta)][1]
             
             dose_est <- dfunc(newdat[,3] / newdat[,2])
             panel.xyplot(dose_est, newdat[,3] / newdat[,2], ...)
             panel.arrows(dose_est, newdat[,3] / newdat[,2], dose_est, 0, length = 0.15, col = "#0080ff")
             panel.abline(v = newdat[,1], col = "red")
             panel.segments(0, yci[1], cil_indx * .001, yci[1], col = "gray", lty = "dashed", lwd = 2)
             panel.segments(cil_indx * .001, yci[1], cil_indx * .001, 0, col = "gray", lty = "dashed", lwd = 2)
             
             panel.segments(0, yci[2], ciu_indx * .001, yci[2], col = "gray", lty = "dashed", lwd = 2)
             panel.segments(ciu_indx * .001, yci[2], ciu_indx * .001, 0, col = "gray", lty = "dashed", lwd = 2)
           }
           
           if (nrow(mod_est) == 3) {
             panel.curve(expr = mod_est[1, 1] + mod_est[2, 1] * x + mod_est[3, 1] * (x ^ 2), from = 0, to = max(doses), lwd = 2, col = "#ff00ff")
             panel.curve(expr = mod_est[1, 1] + mod_est[2, 1] * x + mod_est[3, 1] * (x ^ 2) + sqrt(R) * 
                           sqrt(vcov_mat[1, 1] + vcov_mat[2, 2] * (x ^ 2) + vcov_mat[3, 3] * (x ^ 4) +
                                  2 * vcov_mat[1, 2] * x + 2 * vcov_mat[1, 3] * (x ^ 2) + 2 * vcov_mat[2, 3] * (x ^ 3)), from = 0, to = max(doses), lwd = 2, col = "darkgreen")
             panel.curve(expr = mod_est[1, 1] + mod_est[2, 1] * x + mod_est[3, 1] * (x ^ 2) - sqrt(R) * 
                           sqrt(vcov_mat[1, 1] + vcov_mat[2, 2] * (x ^ 2) + vcov_mat[3, 3] * (x ^ 4) +
                                  2 * vcov_mat[1,2] * x + 2 * vcov_mat[1,3] * (x ^ 2) + 2 * vcov_mat[2, 3] * (x ^ 3)), from = 0, to = max(doses), lwd = 2, col = "darkgreen") 
           } else if (nrow(mod_est) == 2) {
             panel.curve(expr = mod_est[1, 1] + mod_est[2, 1] * x, from = 0, to = max(doses), lwd = 2, col = "#ff00ff")
             panel.curve(expr = mod_est[1, 1] + mod_est[2, 1] * x + sqrt(R) * 
                           sqrt(vcov_mat[1, 1] + vcov_mat[2, 2] * (x ^ 2) + 2 * vcov_mat[1, 2] * x), from = 0, to = max(doses), lwd = 2, col = "darkgreen")
             panel.curve(expr = mod_est[1, 1] + mod_est[2, 1] * x - sqrt(R) * 
                           sqrt(vcov_mat[1, 1] + vcov_mat[2, 2] * (x ^ 2) + 2 * vcov_mat[1, 2] * x), from = 0, to = max(doses), lwd = 2, col = "darkgreen")
           }
         }, ...)
}