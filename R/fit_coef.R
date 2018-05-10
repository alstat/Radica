#' Fit coefficients
#' @export
#' 
#' @description
#' Fit coefficients to the chromosomal aberration data.
#' 
#' @param x a data frame consisting the total number of cells, 
#'          number of aberrations detected, doses and cell distribution
#'          of chromosomal aberration;
#' @param cells column name of the total number of cells in x;
#' @param aberrations column name of the aberrations detected in x;
#' @param doses column name of the doses in x;
#' @param cell.dist a vector of length 2 for the starting and end column of
#'        the cell distribution;
#' @param model a character for type of regression model, default is set to
#'        \code{"lq"}, for linear-quadratic model. And \code{"l"} for linear model;
#' @param dist if set to \code{NULL}, then the data is considered as non-overdispersed,
#'        and thus uses the Poisson family with dispersion = 1. But if set to
#'        \code{"overdispersed"}, then quassi-Possoin family is considered with 
#'        the dispersion not restricted to 1, but is estimated. 
#' @param level a numeric for the confidence level of the fitted values.
#' @param curve logical option for plot of the calibration curve. If set to
#'        \code{TRUE}, calibration curve is generated.
#' @param na.rm logical switch, if set to TRUE, missing values are ommitted from
#'        the computation.
#' 
#' @details
#' The function fits the coefficient into the data using the linear-quadratic
#' model:
#'                 \deqn{Y = C + \alpha*D + \beta*D ^ 2}
#' where \eqn{Y} is the yields of the chromosome aberrations; \eqn{C}, 
#' \eqn{\alpha} and \eqn{\beta} are the coefficients need to be 
#' estimated; and \eqn{D} is the dose.
#' 
#' or, for high LET radiation, the \eqn{\alpha}-term becomes large and 
#' eventually the \eqn{\beta}-term becomes biologically less relevant and
#' also statistically ‘masked’ and the dose response is approximated by 
#' the linear equation
#'                 \deqn{Y = C + \alpha*D}
#' For dicentrics, irradiation with X or gamma rays produces a
#' distribution of damage which is very well represented by the Poisson
#' distribution. In contrast, neutrons and other types of high LET 
#' radiation produce distributions which display overdispersion, where 
#' the variance (\eqn{\sigma^2}) exceeds the mean (\eqn{y}).
#' 
#' Because curve fitting methods are based on Poisson statistics, the dicentric cell
#' distribution should be tested for compliance with the Poisson distribution for each dose used
#' to construct the calibration curve. Nowadays, the most widely used test is the u test.
#' The u test statistic is a normalized unit of the dispersion index (\eqn{\sigma^2/y}), 
#' which for a Poisson distribution should be unity. u values higher than 1.96 indicate 
#' overdispersion (with a two-sided significance level, \eqn{\alpha = 0.025}).
#' 
#' The function returns a list:
#' \itemize{
#' \item{\code{preliminary}}{ a data frame of the original data but with new columns,
#'        for \code{disp_idx} (dispersion index) and \code{u_test} (Poisson u test);}
#' \item{\code{summary}}{ a summary of the GLM model.}
#' }
#' 
#' @examples
#' library(Radica) 
#' 
#' cobalt1 <- fit_coef(Cobalt60, cells = cell, aberrations = aberr,
#'   doses = doses, cell.dist = c(4, 9), curve = TRUE)
#' 
#' # for overdispersion (see if residual deviance is much larger than degrees
#' # of freedom) consider quasi-Poisson family.
#' cobalt2 <- fit_coef(Cobalt60, cells = cell, aberrations = aberr,
#'   doses = doses, cell.dist = c(4, 9), dist = "over-dispersed", curve = TRUE)
#'   
#' # For high Linear Energy Transfer (LET) data, we consider the linear model
#' # with over-dispersed distribution
#' mevhe <- fit_coef(Mev4He, cells = cell, aberrations = aberr, model = "l",
#'   doses = doses, cell.dist = c(4, 11), dist = "over-dispersed", curve = TRUE)
#' 
#' @references
#' EPR-Biodosimetry (2011). \href{http://www-pub.iaea.org/MTCD/publications/PDF/EPR-Biodosimetry%202011_web.pdf}{\emph{Cytogenetic Dosimetry: Applications in Preparedness for and 
#' Response to Radiation Emergencies}.
#' 
#' Quick-R. \href{http://www.statmethods.net/advstats/glm.html}{\emph{Generalized Linear Models}}. Retrieved June 4, 2014.
fit_coef <- function(x, cells, aberrations, doses, cell.dist, model = "lq", start = NULL, 
                     level = 0.95, dist = NULL, curve = FALSE, na.rm = FALSE) {
  
  cells <- eval(substitute(cells), envir = x, enclos = parent.frame())
  aberrations <- eval(substitute(aberrations), envir = x, enclos = parent.frame())
  doses <- eval(substitute(doses), envir = x, enclos = parent.frame())
  x <- u.test(x = x, cells = cells, aberrations = aberrations,
              cell.dist = cell.dist)
  
  x0 <- cells; x1 <- cells * doses; x2 <- cells * doses * doses
  xdata <- data.frame(x0, x1, x2, aberrations)
  if (model == "lq" && is.null(dist)) {
    disp_idx <- rep(mean(x$disp_idx, na.rm = na.rm), nrow(x)); weights <- 1L / disp_idx;
    result <- glm(aberrations ~ -1L + x0 + x1 + x2, family = poisson(link = "identity"), 
                  weights = weights, data = xdata)
  } else if (model == "lq" && !is.null(dist) && dist == "overdispersed") {
    disp_idx <- x$disp_idx; weights <- 1L / disp_idx;
    result <- glm(aberrations ~ -1L + x0 + x1 + x2, family = quasipoisson(link = "identity"), 
                  weights = weights, data = xdata)
  } else if (model == "l" && is.null(dist)) {
    disp_idx <- rep(mean(x$disp_idx, na.rm = na.rm), nrow(x)); weights <- 1L / disp_idx;
    result <- glm(aberrations ~ -1L + x0 + x1, family = poisson(link = "identity"), 
                  weights = weights, data = xdata)
  } else if (model == "l" && !is.null(dist) && dist == "overdispersed") {
    disp_idx <- x$disp_idx; weights <- 1L / disp_idx;
    result <- glm(aberrations ~ -1L + x0 + x1, family = quasipoisson(link = "identity"), 
                  weights = weights, data = xdata)
  }
  
  res.smry <- summary(result, correlation = TRUE)
  chisq.p <- pchisq(deviance(res.smry), df.residual(res.smry), lower.tail = FALSE)
  
  final_out <- list(preliminary = x, summary = res.smry, goodness_of_fit = chisq.p)
  
  if (curve == TRUE) {
    print(dose_curve(final_out, level = level))
  }
  
  return(final_out)
}