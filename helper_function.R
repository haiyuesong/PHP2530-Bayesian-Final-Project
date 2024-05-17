setwd("~/Bayesian Final")
library(data.table)
library(bnlearn)
library(Rgraphviz)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(igraph)
library(psych)
library(parallel)

########## helper functions ##########

# define a function to find high-correlated regions for dimension reduction
high_corr <- function(roi1, roi2, std_fisherz, threshold=0.8,std_err){
  #' find whether the two input ROIs have high-correlated time courses with p-values
  #' @param roi1, the index of first region of interest
  #' @param roi2, the index of second region of interest
  #' @param std_fisherz, standardized fisher r-to-z transformed correlation matrix
  #' @param threshold, the threshold of the correlation (default 0.8)
  #' @param std_err, the standard error after fisher r-to-z transformation
  #' @return the p-value of testing the two input ROIs have high-correlated time courses
  
  std_fisherz_threshold <- fisherz(threshold)/std_err
  z_diff <- std_fisherz[roi1,roi2]-std_fisherz_threshold
  p_value <- 1 - pnorm(z_diff) # one-tailed test with null hypothesis: correlation < threshold
  return(p_value)
}

# Bayesian Network
# Function to convert bnlearn object to igraph object
bn_to_igraph <- function(bn) {
  arcs <- arcs(bn)
  nodes <- unique(c(arcs[, 1], arcs[, 2]))
  graph <- graph_from_data_frame(arcs, directed = TRUE, vertices = nodes)
  return(graph)
}

# Function to check if a proposed structure is valid (does not contain cycles)
is_valid_structure <- function(new_structure) {
  igraph_graph <- bn_to_igraph(new_structure)
  return(is_dag(igraph_graph))
}

# Function to check if adding an edge creates a cycle
would_create_cycle <- function(graph, from, to) {
  igraph_graph <- bn_to_igraph(graph)
  igraph_graph <- add_edges(igraph_graph, c(from, to))
  return(!is_dag(igraph_graph))
}

# Function to check if there is a path between two nodes (avoiding cycles)
has_path <- function(graph, from, to) {
  igraph_graph <- bn_to_igraph(graph)
  paths <- all_shortest_paths(igraph_graph, from, to)
  return(length(paths$res) > 0)
}

# Function to check if an arc exists in the network
arc_exists <- function(graph, from, to) {
  arcs_list <- arcs(graph)
  return(any(arcs_list[, 1] == from & arcs_list[, 2] == to))
}

# Calculate the log-likelihood score
loglik_score <- function(network, data) {
  return(score(network, data, type = "loglik-g"))
}

# Metropolis-Hastings Algorithm for Bayesian Network using log-likelihood
metropolis_hastings_bn <- function(data, initial_structure, iterations = 10000) {
  current_structure <- initial_structure
  best_score <- loglik_score(current_structure, data)
  loglik_scores <- numeric(iterations)

  for (i in 1:iterations) {
    cat("Iteration:", i, "\n")
    # Propose a new structure
    new_structure <- propose_new_structure(current_structure)

    # Check if the proposed structure is valid (no cycles)
    if (!is_valid_structure(new_structure)) {
      cat("Invalid structure due to cycle.\n")
      next
    }

    new_score <- loglik_score(new_structure, data)
    loglik_scores[i] <- new_score
    cat("Current score:", best_score, "New score:", new_score, "\n")

    # Calculate the acceptance probability
    acceptance_prob <- exp(new_score - best_score)
    cat("Acceptance probability:", acceptance_prob, "\n")

    # Accept or reject the new structure
    if (runif(1) < acceptance_prob) {
      current_structure <- new_structure
      best_score <- new_score
      cat("Accepted new structure.\n")
    } else {
      cat("Rejected new structure.\n")
      loglik_scores[i] <- best_score
    }
  }

  return(list(model = current_structure, loglik_scores = loglik_scores))
}

# Enhanced proposal function to ensure diverse and valid proposals
propose_new_structure <- function(current_structure, max_attempts = 100) {
  new_structure <- current_structure

  # Randomly choose to add, remove, or reverse an edge
  choice <- sample(c("add", "remove", "reverse"), 1)

  # Get the nodes in the network
  nodes <- nodes(current_structure)

  attempt <- 0
  repeat {
    if (choice == "add") {
      # Randomly add an edge
      from <- sample(nodes, 1)
      to <- sample(nodes, 1)
      if (from != to && !arc_exists(new_structure, from, to) && !would_create_cycle(new_structure, from, to)) {
        new_structure <- set.arc(new_structure, from, to)
        break
      }
    } else if (choice == "remove") {
      # Randomly remove an edge
      arcs_list <- arcs(current_structure)
      if (nrow(arcs_list) > 0) {
        arc <- arcs_list[sample(1:nrow(arcs_list), 1), ]
        new_structure <- drop.arc(new_structure, arc[1], arc[2])
        break
      }
    } else if (choice == "reverse") {
      # Randomly reverse an edge
      arcs_list <- arcs(current_structure)
      if (nrow(arcs_list) > 0) {
        arc <- arcs_list[sample(1:nrow(arcs_list), 1), ]
        if (!has_path(new_structure, arc[2], arc[1]) && !would_create_cycle(new_structure, arc[2], arc[1])) {
          new_structure <- reverse.arc(new_structure, arc[1], arc[2])
          break
        }
      }
    }
    attempt <- attempt + 1
    if (attempt >= max_attempts) {
      cat("Max attempts reached. Skipping this proposal.\n")
      break
    }
  }

  return(new_structure)
}

# Function for running multiple chains to assess convergence
run_multiple_chains <- function(data, nodes, iterations, num_chains = 3, burn_in) {
  chains <- list()
  
  for (i in 1:num_chains) {
    initial_structure <- random.graph(nodes, method = "ic-dag", max.in.degree = 2)
    result <- metropolis_hastings_bn(data, initial_structure, iterations)
    chains[[i]] <- result$loglik_scores[(burn_in + 1):iterations]
  }
  
  return(mcmc.list(lapply(chains, mcmc)))
}
