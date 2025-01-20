#' Sample Phylogenetic Tree
#'
#' Samples and keeps leaves with probability pi.
#'
#' @param ptree object of class phylo.
#' @param time how many time points in the past the most recent leaf is.
#' @param pi0 probability of sampling a leaf at the present time.
#' @param pi1 probability of sampling a leaf before the present time.
#'
#' @return object of class phylo; the tree keeping only the sampled leaves.
#' @export
#'
#' @examples
#' sample_phylo(ptree = full_tree, pi0 = 0, pi1 = 0.01)
sample_phylo <- function(ptree, ptree_lag=0, pi0, pi1) {
  if (pi0 < 0 | pi0 > 1 | pi1 < 0 | pi1 > 1) {
    warning("pi must be between 0 and 1")
    break
  }

  #number of leaves
  n_leaves <- length(ptree$tip.label)

  distance <- distToRoot(ptree)
  max <- max(distance)

  #calculate how much to add to edge lengths to make leaves end on particular times
  change <- (max - distance) %% 1
  vector <- (change > 0.99999999 | change < 0.00000001)
  change[vector] <- 0

  #add on those differences
  edges <- (ptree$edge[,2] <= n_leaves)
  ptree$edge.length[edges] <- ptree$edge.length[edges] + change

  new_distance <- distToRoot(ptree)
  new_max <- max(new_distance)

  if (ptree_lag > 0) {
    time1 <- 1:n_leaves
    time0 <- NULL
  } else {
    time1 <- (1:n_leaves)[round(new_max - new_distance, 0) > 0]
    time0 <- (1:n_leaves)[-time1]
  }
  #sample leaves
  u0 <- runif(length(time0))
  keep0 <- time0[u0 <= pi0]
  u1 <- runif(length(time1))
  keep1 <- time1[u1 <= pi1]
  keep <- c(keep0, keep1)

  if (!length(keep)) {
    new_tree <- NULL
    ptree_lag <- NULL
  } else {
    if (!length(keep0)) {
      ptree_lag <- min((new_max - new_distance)[keep1]) + ptree_lag
    } else {
      ptree_lag <- ptree_lag
    }
    new_tree <- ape::keep.tip(ptree, keep)
    class(new_tree) <- "phylo"
  }

  output <- list("ptree"=new_tree, "ptree_lag"=ptree_lag)
  return(output)
}
