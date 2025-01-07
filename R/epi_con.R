#' Simulate epidemics
#'
#' Simulates a birth-death epidemic with constant birth and death rates and x0 infected on day 0.
#'
#' @param birth_rate birth rate of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param stop_time number of days to run the epidemic simulation for.
#' @param x0 initial number of infectives.
#'
#' @return data frame with: Index, Parent Index, Infection Status, Infection Time, Removal Time
#' @export
#'
#' @examples
#' epi_con(birth_rate = 0.2, death_rate = 0.1, stop_time = 50)
epi_con <- function(birth_rate,death_rate,stop_time,x0=1) {
  #index labels everyone by order of entry
  #parent is the number of parent
  #status is infection status: I is infected, R is removed
  #inf_time is infection time
  #rem_time is removal time
  data <- list("index"=1,"parent"=0,"status"="I","inf_time"=0,"rem_time"=NA)
  if (x0 > 1) {
    for (i in 2:x0) {
      data$index[i] <- i
      data$parent[i] <- 0
      data$status[i] <- "I"
      data$inf_time[i] <- 0
      data$rem_time[i] <- NA
    }
  }

  n_inf <- x0 #number infected
  n_total <- x0 #total = number infected + number recovered
  current_time <- 0 #current time, initially 0

  inf <- 1:x0 #vector of infected people

  a <- birth_rate / (birth_rate + death_rate) #acceptance probability

  #run loop while some people are still infected
  while (n_inf > 0) {
    #simulate event time
    event_time <- rexp(1, rate = n_inf * (birth_rate + death_rate))

    #if time is past stop time, then stop
    current_time <- current_time + event_time
    if (current_time > stop_time) {
      break
    }

    #sample individual to be infected or recover
    #note: if sample is of length 1, then it samples from 1:x
    #so need to do case where length is 1 separately
    if (n_inf == 1) {
      i <- inf
    } else {
      i <- sample(inf, size = 1)
    }

    u <- runif(1)
    if (u < a) {
      #infection event
      n_inf <- n_inf + 1 #number infected increased
      n_total <- n_total + 1 #total number increases
      inf[n_inf] <- n_total #add m to n.inf vector

      #add row for new infected individual
      data$index[n_total] <- n_total
      data$parent[n_total] <- i
      data$status[n_total] <- "I"
      data$inf_time[n_total] <- current_time
      data$rem_time[n_total] <- NA
    } else {
      #recovery event
      data$status[i] <- "R"
      data$rem_time[i] <- current_time
      n_inf <- n_inf - 1
      inf <- inf[!inf==i] #remove i from inf
    }
  }

  data <- as.data.frame(data)
  return(data)
}
