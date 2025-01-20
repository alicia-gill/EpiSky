#' Epidemic functions
#'
#' @name epi_fns
#' @title A collection of functions for the birth rate in used to generate epidemics
#'
#' @param t time
#'
#' @return function
NULL

#' @rdname epi_fns
#' @export
constant <- function(t) {
  return(0.3)
}

#' @rdname epi_fns
#' @export
change_point <- function(t) {
  change_time <- 10
  before <- 0.5
  after <- 0.1

  if (t < change_time) {
    return(before)
  } else {
    return(after)
  }
}

#' @rdname epi_fns
#' @export
peak <- function(t) {
  stop_time <- 40
  peak_time <- 20
  min <- 0.1
  max <- 0.3

  if (t < peak_time) {
    slope <- (max-min)/(peak_time)
    return(min + slope*t)
  } else {
    slope <- (min-max)/(stop_time-peak_time)
    return(max + slope*(t-peak_time))
  }
}

#' @rdname epi_fns
#' @export
increase <- function(t) {
  stop_time <- 20
  start <- 0.1
  end <- 0.5

  slope <- (start-end)/stop_time
  return(start + slope*t)
}

#' @rdname epi_fns
#' @export
decrease <- function(t) {
  stop_time <- 20
  start <- 0.5
  end <- 0.1

  slope <- (start-end)/stop_time
  return(start + slope*t)
}
