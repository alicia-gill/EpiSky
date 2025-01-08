#' Sample Prevalence
#'
#' Generates a partially observed prevalence.
#'
#' @param prevalence data frame of prevalence per time.
#' @param reporting_prob proportion of cases observed.
#'
#' @return data frame of observed prevalence per day.
#' @export
#'
#' @examples
#' sample_prevalence(prevalence = prev, reporting_prob = 0.2)
sample_prevalence <- function(prevalence, reporting_prob) {
  n <- nrow(prevalence)
  noisy_prev <- as.list(prevalence)
  noisy_prev$prevalence <- rbinom(n = n, size = noisy_prev$prevalence, prob = reporting_prob)
  noisy_prev <- as.data.frame(noisy_prev)
  return(noisy_prev)
}
