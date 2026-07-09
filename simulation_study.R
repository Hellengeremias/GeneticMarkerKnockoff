# ==============================================================================
# Simulation study: knockoff generators for genetic marker selection
# ==============================================================================
#
# This script reproduces the simulation workflow used to compare LASSO with
# knockoff-based variable selection under different dependence structures among
# genetic markers.
#
# Main steps:
#   1. Load and preprocess genotype data.
#   2. Select the dependence scenario.
#   3. Fit or load the HMM parameters used by SNPknock.
#   4. Estimate the probability models required by the proposed generators.
#   5. Simulate multiple response vectors.
#   6. Run the knockoff filter and the LASSO benchmark.
#   7. Summarize selection metrics and produce the comparison figure.
#
# NOTE:
# - Controlled-access ISA-Nutrition genotype data are not distributed with this
#   repository. Set ISA_DATA_DIR to the corresponding local directory.
# - fastPHASE is external software. Set SOFTWARE_DIR to the directory containing
#   the fastPHASE executable and its input/output files.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(tidyverse)
library(glmnet)
library(knockoff)
library(SNPknock)
library(gridExtra)
library(compiler)
library(forcats)


# ------------------------------------------------------------------------------
# 2. User configuration
# ------------------------------------------------------------------------------

# Available scenarios: "low_ld", "moderate_ld", and "high_ld".
ld_scenario <- "moderate_ld"

# Number of simulated responses.
M <- 200

# Number of truly associated markers.
n_signals <- 10

# Target FDR level for the knockoff+ filter.
target_fdr <- 0.20

# Window half-width for the classification-tree generator.
tree_window <- 5

# Local paths. Environment variables can be used to override these defaults.
isa_data_dir <- Sys.getenv(
  "ISA_DATA_DIR",
  unset = "/home/hellen/ISA/Bases"
)

software_dir <- Sys.getenv(
  "SOFTWARE_DIR",
  unset = "/home/hellen/Programas"
)

mouse_genotype_file <- Sys.getenv(
  "MOUSE_GENOTYPE_FILE",
  unset = "dados_densidade_ossea_sem_missing.csv"
)

mouse_map_file <- Sys.getenv(
  "MOUSE_MAP_FILE",
  unset = "Data_Description_NZBxRF_Wergedal2006.xlsx"
)

output_dir <- Sys.getenv(
  "OUTPUT_DIR",
  unset = "."
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 3. Helper functions for genotype preprocessing
# ------------------------------------------------------------------------------

calc_maf <- function(x) {
  x <- x[!is.na(x)]
  allele_frequency <- mean(x) / 2
  min(allele_frequency, 1 - allele_frequency)
}


select_by_clustering <- function(G, k = 100) {
  # Convert absolute correlations into distances and form k clusters.
  correlation_matrix <- cor(G, use = "pairwise.complete.obs")
  distance_matrix <- as.dist(1 - abs(correlation_matrix))
  hc <- hclust(distance_matrix, method = "average")
  cluster_id <- cutree(hc, k = k)

  # Select the marker with the largest variance within each cluster.
  marker_variance <- apply(G, 2, var, na.rm = TRUE)
  representatives <- tapply(
    seq_along(cluster_id),
    cluster_id,
    function(idx) idx[which.max(marker_variance[idx])]
  )

  selected_idx <- sort(unlist(representatives, use.names = FALSE))

  list(
    index = selected_idx,
    names = colnames(G)[selected_idx]
  )
}


# ------------------------------------------------------------------------------
# 4. Genotype data and dependence scenarios
# ------------------------------------------------------------------------------

if (ld_scenario %in% c("low_ld", "moderate_ld")) {

  # SNP map for chromosome 22.
  map <- read.table(
    file.path(isa_data_dir, "isa_rec_rs_aut_hg37_22.map"),
    header = FALSE
  )
  map <- map[order(map$V4), ]

  # Genomic relationship matrix.
  grm <- read.table(
    file.path(isa_data_dir, "isa.grm841x841.domicilio.grm"),
    header = FALSE
  )

  # Relatedness coefficient thresholds:
  # > 0.707: duplicate/MZ twin
  # [0.354, 0.707]: first-degree relatives
  # [0.177, 0.354]: second-degree relatives
  # [0.0884, 0.177]: third-degree relatives

  relatedness_threshold <- 0.354

  grm_diagonal <- grm %>%
    dplyr::filter(V1 == V2)

  related_pairs <- grm %>%
    dplyr::filter(V1 != V2)

  related_ids <- grm %>%
    dplyr::filter(V4 > relatedness_threshold, V1 != V2)

  diagonal_ids <- grm %>%
    dplyr::filter(V4 > relatedness_threshold, V1 == V2)

  unrelated_ids <- dplyr::anti_join(
    diagonal_ids,
    related_ids,
    by = "V2"
  )

  genotype_data <- read.table(
    file.path(isa_data_dir, "SNPs_chr22_total_imput_fam.txt"),
    header = TRUE
  )

  genotype_data <- genotype_data[
    row.names(genotype_data) %in% as.character(unrelated_ids$V2),
    ,
    drop = FALSE
  ]

  # Retain common variants (MAF >= 0.05).
  mafs <- apply(genotype_data, 2, calc_maf)
  common_genotypes <- genotype_data[, mafs >= 0.05, drop = FALSE]

  if (ld_scenario == "low_ld") {
    set.seed(1)
    selected_markers <- select_by_clustering(common_genotypes, k = 100)
    X <- as.matrix(
      common_genotypes[, selected_markers$index, drop = FALSE]
    )
  }

  if (ld_scenario == "moderate_ld") {
    if (ncol(common_genotypes) < 1000) {
      stop(
        "The moderate-LD scenario requires at least 1000 common SNPs."
      )
    }

    X <- as.matrix(
      common_genotypes[, seq_len(1000), drop = FALSE]
    )
  }

} else if (ld_scenario == "high_ld") {

  # High-LD genotype data.
  mouse_data <- read.csv(mouse_genotype_file)
  X <- as.matrix(mouse_data[, -1, drop = FALSE])

  marker_info <- readxl::read_xlsx(mouse_map_file)
  colnames(X) <- marker_info$marker_name[seq_len(ncol(X))]

  # Uncomment the following recoding if the genotype matrix is coded as -1/0/1
  # and must be converted to 0/1/2:
  #
  # X[X == 0] <- 3
  # X[X == -1] <- 0
  # X[X == 1] <- 2
  # X[X == 3] <- 1

} else {
  stop(
    "`ld_scenario` must be one of: 'low_ld', 'moderate_ld', or 'high_ld'."
  )
}

storage.mode(X) <- "numeric"

cat(
  "Scenario:", ld_scenario,
  "\nIndividuals:", nrow(X),
  "\nMarkers:", ncol(X), "\n"
)


# ------------------------------------------------------------------------------
# 5. HMM parameter estimation for the benchmark generator
# ------------------------------------------------------------------------------

# Second-order Model-X knockoffs do not require a separate preprocessing step.
# HMM knockoffs require fastPHASE parameter estimates.

scenario_tag <- switch(
  ld_scenario,
  low_ld = "dep100",
  moderate_ld = "dep1000",
  high_ld = "high_ld"
)

fastphase_input <- file.path(
  software_dir,
  paste0("Xfastphase_", scenario_tag, ".inp")
)

fastphase_output <- file.path(
  software_dir,
  paste0("hmm_param_", scenario_tag)
)

rhat_file <- paste0(fastphase_output, "_rhat.txt")
alpha_file <- paste0(fastphase_output, "_alphahat.txt")
theta_file <- paste0(fastphase_output, "_thetahat.txt")

hmm_parameter_files <- c(
  rhat_file,
  alpha_file,
  theta_file
)

# Run fastPHASE only when the parameter files are not already available.
if (!all(file.exists(hmm_parameter_files))) {
  writeXtoInp(
    X,
    phased = FALSE,
    out_file = fastphase_input
  )

  runFastPhase(
    fp_path = file.path(software_dir, "fastPHASE"),
    X_file = fastphase_input,
    out_path = fastphase_output,
    K = 12
  )
}

rk <- t(
  read.table(
    rhat_file,
    sep = "/",
    skip = 1
  )
)

alphak <- as.matrix(
  read.table(
    alpha_file,
    skip = 1
  )
)

thetak <- as.matrix(
  read.table(
    theta_file,
    skip = 1
  )
)


# ------------------------------------------------------------------------------
# 6. Proposed knockoff generators
# ------------------------------------------------------------------------------

# ---------- Helper functions ----------

sample_discrete <- function(p) {
  u <- runif(1)
  P <- cumsum(p)
  val <- sum(P < u) + 1
  return(val)
}

haldane_recombination <- function(dist) 0.5 * (1 - exp(-2 * dist)) # dist: genetic distance in Morgans

prob_equal_genotype <- function(recomb, gen1) if (gen1 == 1) (recomb^2) + (1 - recomb)^2 else (1 - recomb)^2
prob_one_step_difference <- function(r) 2 * (1 - r) * r
prob_two_step_difference <- function(r) r^2

genotype_transition_prob <- function(g1, g2, r) {
  d <- abs(g1 - g2)
  if (d == 0) prob_equal_genotype(r, g1)
  else if (d == 1) if (g1 == 1) prob_one_step_difference(r) / 2 else prob_one_step_difference(r)
  else prob_two_step_difference(r)
}

# Complete probability matrices with fewer than three genotype categories
complete_probability_matrix <- function(mat, observed_genotypes) {
  
  # Expected genotype levels
  gen <- c(0, 1, 2)
  
  # Number of columns
  n <- ncol(mat)
  
  # Initialize the completed matrix
  completed_mat <- matrix(0, nrow = length(gen), ncol = n)
  rownames(completed_mat) <- as.character(gen)
  
  # Copy the observed probability rows
  rownames(mat) <- as.character(observed_genotypes)
  common_rows <- intersect(rownames(mat), rownames(completed_mat))
  completed_mat[common_rows, ] <- mat[common_rows, ]
  
  return(completed_mat)
}

# ------------------------------------------------------------------------------
# 6.1 Haldane-based generator
# ------------------------------------------------------------------------------
estimate_haldane_probabilities <- function(X, marker_position) {
  
  n <- nrow(X)
  p <- ncol(X)
  gen <- c(0, 1, 2)
  
  # List of genotype-probability matrices
  prob_list <- vector("list", p)
  
  # Knockoff matrix
  X_k <- matrix(NA, nrow = n, ncol = p)  
  
  # Generate knockoff genotypes
  enableJIT(3)
  # Iterate over markers to generate their knockoff counterparts
  for (j in 1:p) {
    probs_j <- matrix(NA, nrow = length(gen), ncol = n)
    
    for (i in 1:n) {
      probs <- numeric(length(gen))
      
      if (j == 1) {
        
        # First marker: condition only on the marker to the right
        g_right <- X[i, j + 1]
        r1 <- haldane_recombination(abs(marker_position[2] - marker_position[1]))
        probs <- sapply(gen, function(g) max(genotype_transition_prob(g_right, g, r1), 1e-6))
        
      } else if (j == p) {
        # Last marker: condition only on the marker to the left
        g_left <- X[i, j - 1]
        r1 <- haldane_recombination(abs(marker_position[j] - marker_position[j - 1]))
        probs <- sapply(gen, function(g) max(genotype_transition_prob(g_left, g, r1), 1e-6))
        
      } else if (j > 1 & j < p) {
        # Intermediate markers
        
        # Marker to the left
        g_left <- X[i, j - 1]
        
        # Marker to the right
        g_right <- X[i, j + 1]
        
        r1 <- haldane_recombination(abs(marker_position[j] - marker_position[j - 1]))
        r2 <- haldane_recombination(abs(marker_position[j + 1] - marker_position[j]))
        r12 <- haldane_recombination(abs(marker_position[j + 1] - marker_position[j - 1]))
        
        denom <- max(genotype_transition_prob(g_left, g_right, r12), 1e-6)
        
        probs <- sapply(gen, function(g) {
          p1 <- max(genotype_transition_prob(g_left, g, r1), 1e-6)
          p2 <- max(genotype_transition_prob(g, g_right, r2), 1e-6)
          (p1 * p2) / denom
        })
      } 
      
      # Normalize probabilities and sample a genotype
      probs <- probs / sum(probs)
      X_k[i, j] <- gen[sample_discrete(probs)]
      probs_j[, i] <- probs
    }
    prob_list[[j]] <- probs_j
  }
  
  # Interleave original and knockoff variables
  SNPs <- matrix(NA, nrow = n, ncol = 2 * p)
  SNPs[, seq(1, 2 * p, by = 2)] <- X
  SNPs[, seq(2, 2 * p, by = 2)] <- X_k
  
  return(list(
    SNPs = SNPs,
    prob_list = prob_list
  ))
}

# ------------------------------------------------------------------------------
# 6.2 Empirical Markov generator
# ------------------------------------------------------------------------------

estimate_empirical_markov_probabilities <- function(X) {
  
  p <- ncol(X)
  n <- nrow(X)
  gen <- c(0,1,2)
  
  # List of genotype-probability matrices
  prob_list <- vector("list", p)
  
  # Knockoff matrix
  X_k <- matrix(NA, nrow = n, ncol = p)  
  
  enableJIT(3)
  for (j in 1:p) {
    probs_j <- matrix(NA, nrow = length(gen), ncol = n)
    
    for (i in 1:n) {
      probs <- numeric(length(gen))
      
      for (v in seq_along(gen)) {
        
        # Candidate genotypes for the current marker: 0, 1, and 2
        x_current <- gen[v]
        
        if (j == 1) {
          # First SNP
          x_right <- X[i, j + 1]
          conj <- sum(X[, j] == x_current & X[, j + 1] == x_right) # P(X_right = x_right, X_current = x_current)
          marg <- sum(X[, j + 1] == x_right)                 # P(X_right = x_right)
          probs[v] <- ifelse(marg > 0, conj / marg, 0)
          
        } else if (j == p) {
          # Last SNP
          x_left <- X[i, j - 1]
          conj <- sum(X[, j] == x_current & X[, j - 1] == x_left) # P(X_left = x_left, X_current = x_current)
          marg <- sum(X[, j - 1] == x_left)                 # P(X_left = x_left)
          probs[v] <- ifelse(marg > 0, conj / marg, 0)
          
        } else {
          # Intermediate markers (j > 1 and j < p)
          x_left <- X[i, j - 1]
          x_right <- X[i, j + 1]
          
          conj1 <- sum(X[, j] == x_current & X[, j - 1] == x_left)  # P(X_current = x_current, X_left = x_left)
          marg1 <- sum(X[, j-1] == x_left)                    # P(X_left = x_left)
          p_current_given_left <- ifelse(marg1 > 0, conj1 / marg1, 0)
          
          conj2 <- sum(X[, j + 1] == x_right & X[, j] == x_current)  # P(X_right = x_right, X_current = x_current)
          marg2 <- sum(X[, j] == x_current)                      # P(X_left = x_left)
          p_right_given_current <- ifelse(marg2 > 0, conj2 / marg2, 0)
          
          conj3 <- sum(X[, j + 1] == x_right & X[, j - 1] == x_left) # P(X_right = x_right, X_left = x_left)
          marg3 <- sum(X[, j - 1] == x_left)                     # P(X_left = x_left)
          p_right_given_left <- ifelse(marg3 > 0, conj3 / marg3, 1e-10)  # avoid division by zero
          
          probs[v] <- p_current_given_left * p_right_given_current / p_right_given_left
        }
      }
      
      # Normalize probabilities and sample a genotype
      probs <- probs / sum(probs)
      X_k[i, j] <- gen[sample_discrete(probs)]
      probs_j[, i] <- probs
    }
    prob_list[[j]] <- probs_j
  }
  
  # Interleave original and knockoff variables
  SNPs <- matrix(NA, nrow = n, ncol = 2*p)
  SNPs[, seq(1, 2*p, by = 2)] <- X
  SNPs[, seq(2, 2*p, by = 2)] <- X_k
  
  return(list(
    SNPs = SNPs,
    prob_list = prob_list
  ))
}

# ------------------------------------------------------------------------------
# 6.3 Classification-tree generator
# ------------------------------------------------------------------------------

estimate_tree_probabilities <- function(X, w) {
  
  p <- ncol(X)
  n <- nrow(X)
  gen <- c(0,1,2)
  
  prob_list <- vector("list", p)
  
  # Knockoff matrix
  X_k <- matrix(NA, nrow = n, ncol = p)
  
  enableJIT(3)
  
  for (j in 1:p) {
    
    #---------------------------------------------------------------------------
    # Fit a classification tree for the j-th marker
    # using neighboring markers within a window of total size 2*w
    target <- factor(X[, j], levels = gen)
    
    if (j == 1) {
      # First SNP: use up to 2*w markers to the right 
      end_idx <- min(j + 2*w, p)
      predictors <- as.data.frame(X[, (j+1):end_idx])
    } else if (j == p) {
      # Last SNP: use up to 2*w markers to the left
      start_idx <- max(j - 2*w, 1)
      predictors <- as.data.frame(X[, start_idx:(j-1)])
    } else {
      # Intermediate SNPs: use up to w markers on each side of x_j 
      start_idx <- max(j - w, 1)
      X_left <- X[, start_idx:(j-1)]
      
      end_idx <- min(j + w, p)
      X_right <- X[, (j+1):end_idx]
      predictors <- as.data.frame(cbind(X_left, X_right))
    }
    
    formula <- as.formula("target ~ .")
    tree_fit <- rpart::rpart(formula, data = cbind(target, predictors),
                           method = "class")
    
    # Predict genotype probabilities for the j-th marker given the predictors
    probs_rpart <- predict(tree_fit, newdata = predictors, type = "prob")
    probs_rpart_t <- t(probs_rpart)
    
    # Add zero-probability rows when fewer than three genotype classes are observed
    if(nrow(probs_rpart_t) <3){
      rowN <- rownames(probs_rpart_t)
      l <- complete_probability_matrix(mat = probs_rpart_t, observed_genotypes = rowN)
    }else{
      l <- probs_rpart_t
    }
    probs_rpart_t2 <- t(matrix(l, ncol = n, nrow = 3))
    
    # Normalize probabilities and sample a genotype
    row_totals <- rowSums(probs_rpart_t2)
    probs_f <- probs_rpart_t2 / row_totals
    for(i in 1:n){
      X_k[i, j] <- gen[sample_discrete(probs_f[i,])]  
    }
    prob_list[[j]] <- t(probs_f)
  }
  
  # Interleave original and knockoff variables
  SNPs <- matrix(NA, nrow = n, ncol = 2*p)
  SNPs[, seq(1, 2*p, by = 2)] <- X
  SNPs[, seq(2, 2*p, by = 2)] <- X_k
  
  return(list(
    SNPs = SNPs,
    prob_list = prob_list
  ))
}

# ------------------------------------------------------------------------------
# 7. Precompute genotype probabilities for the proposed generators
# ------------------------------------------------------------------------------

if (ld_scenario == "high_ld") {
  # Genetic map positions are provided in centimorgans and converted to Morgans.
  marker_position <- marker_info$original_cM[seq_len(ncol(X))] / 100
} else {
  # Match map positions to X in the exact column order used in the analysis.
  snp_names <- sub("_.*", "", colnames(X))
  map_index <- match(snp_names, as.character(map$V2))

  if (anyNA(map_index)) {
    stop(
      "At least one SNP in X could not be matched to the chromosome map."
    )
  }

  marker_position <- map$V3[map_index]
}

set.seed(1)
haldane_model <- estimate_haldane_probabilities(
  X = X,
  marker_position = marker_position
)
prob_haldane <- haldane_model$prob_list

set.seed(1)
empirical_markov_model <- estimate_empirical_markov_probabilities(X = X)
prob_empirical_markov <- empirical_markov_model$prob_list

set.seed(1)
tree_model <- estimate_tree_probabilities(
  X = X,
  w = tree_window
)
prob_tree <- tree_model$prob_list

stopifnot(
  !anyNA(haldane_model$SNPs),
  !anyNA(empirical_markov_model$SNPs),
  !anyNA(tree_model$SNPs)
)


# ------------------------------------------------------------------------------
# 8. Response simulation
# ------------------------------------------------------------------------------

set.seed(1)
signal_index <- sample(seq_len(ncol(X)), n_signals)

beta <- rep(0, ncol(X))

set.seed(1001)
beta[signal_index] <- rnorm(
  n_signals,
  mean = 0.75,
  sd = 0.5
)

beta_signal <- beta[signal_index[order(signal_index)]]
names(beta_signal) <- signal_index[order(signal_index)]

beta_signal_df <- data.frame(
  Variable = signal_index,
  beta = beta[signal_index]
)

n <- nrow(X)
linear_predictor <- as.numeric(X %*% beta)

base_seed <- 123
set.seed(base_seed)
response_seeds <- sample.int(1e8, M)

Y_multi <- sapply(response_seeds, function(seed) {
  set.seed(seed)
  linear_predictor + rnorm(n)
})


# ------------------------------------------------------------------------------
# 9. Knockoff filter simulation
# ------------------------------------------------------------------------------
# ---------- Helper functions ----------

# Compute W using glmnet coefficient differences:
# W_j = |beta_j| - |beta_tilde_j|.
compute_glmnet_W <- function(
    X,
    Xk,
    y,
    family = c("gaussian", "binomial")
) {
  family <- match.arg(family)

  W <- knockoff::stat.glmnet_coefdiff(
    X = X,
    Xk = Xk,
    y = y,
    family = family
  )

  as.numeric(W)
}


# Compute the knockoff+ threshold.
compute_knockoff_threshold <- function(
    W,
    alpha_kn = 0.20,
    offset = 1
) {
  candidate_thresholds <- sort(unique(W[W > 0]))

  if (length(candidate_thresholds) == 0) {
    return(Inf)
  }

  for (threshold in candidate_thresholds) {
    n_negative <- sum(W <= -threshold)
    n_positive <- sum(W >= threshold)

    estimated_fdp <- (offset + n_negative) / max(1, n_positive)

    if (estimated_fdp <= alpha_kn) {
      return(threshold)
    }
  }

  Inf
}


# Check whether at least one column is constant.
has_constant_column <- function(M) {
  any(
    apply(
      M,
      2,
      function(z) var(z, na.rm = TRUE) == 0
    )
  )
}


# Generate knockoffs from precomputed genotype probabilities.
generate_from_probabilities <- function(X, prob) {
  p <- ncol(X)
  n <- nrow(X)
  genotypes <- c(0, 1, 2)

  Xk <- matrix(NA, nrow = n, ncol = p)

  enableJIT(3)

  for (j in seq_len(p)) {
    probs_j <- prob[[j]]

    for (i in seq_len(n)) {
      Xk[i, j] <- genotypes[
        sample_discrete(probs_j[, i])
      ]
    }
  }

  Xk
}


# Generate knockoff variables using one of the evaluated generators.
generate_knockoffs <- function(
    X,
    generator = c(
      "second_order",
      "HMM",
      "haldane",
      "empirical_markov",
      "tree"
    ),
    seed = NULL
) {
  generator <- match.arg(generator)

  if (generator == "second_order") {
    return(knockoff::create.second_order(X))
  }

  current_seed <- seed

  repeat {
    if (generator == "HMM") {
      Xk <- SNPknock::knockoffGenotypes(
        X = X,
        r = rk,
        alpha = alphak,
        theta = thetak,
        seed = current_seed,
        cluster = NULL,
        display_progress = FALSE
      )
    } else {
      set.seed(current_seed)

      probability_model <- switch(
        generator,
        haldane = prob_haldane,
        empirical_markov = prob_empirical_markov,
        tree = prob_tree
      )

      Xk <- generate_from_probabilities(
        X = X,
        prob = probability_model
      )
    }

    if (!has_constant_column(Xk)) {
      return(Xk)
    }

    current_seed <- sample.int(1e8, 1L)
  }
}


# Compute variable-selection performance metrics.
selection_metrics <- function(selected, true_signals, p) {
  selected <- sort(unique(selected))

  R <- length(selected)
  s <- length(true_signals)

  TP <- length(intersect(selected, true_signals))
  FP <- R - TP
  FN <- s - TP
  TN <- (p - s) - FP

  FDR <- if (R == 0) 0 else FP / R
  Power <- if (s == 0) NA_real_ else TP / s

  # Coerce terms to double precision before multiplication to avoid integer
  # overflow in large-p settings.
  mcc_terms <- as.numeric(
    c(TP + FP, TP + FN, TN + FP, TN + FN)
  )
  mcc_denominator <- sqrt(prod(mcc_terms))

  MCC <- if (
    !is.finite(mcc_denominator) || mcc_denominator == 0
  ) {
    NA_real_
  } else {
    (TP * TN - FP * FN) / mcc_denominator
  }

  f1_denominator <- 2 * TP + FP + FN
  F1 <- if (f1_denominator == 0) {
    NA_real_
  } else {
    (2 * TP) / f1_denominator
  }

  list(
    FP = FP,
    FN = FN,
    TP = TP,
    TN = TN,
    R = R,
    s = s,
    FDR = FDR,
    Power = Power,
    MCC = MCC,
    F1 = F1
  )
}


# ---------- Main knockoff simulation function ----------

run_knockoff_simulation <- function(
    X,
    Y_multi,
    true_signals,
    family = c("gaussian", "binomial"),
    alpha_kn = 0.20,
    offset = 1,
    generator = c(
      "second_order",
      "HMM",
      "haldane",
      "empirical_markov",
      "tree"
    )
) {
  family <- match.arg(family)
  generator <- match.arg(generator)

  p <- ncol(X)
  M_local <- ncol(Y_multi)

  knockoff_seeds <- sample.int(1e8, M_local)

  W_matrix <- matrix(
    NA_real_,
    nrow = M_local,
    ncol = p
  )

  threshold_vector <- rep(NA_real_, M_local)
  selected_list <- vector("list", M_local)

  per_round <- data.frame(
    round = seq_len(M_local),
    R = NA_integer_,
    TP = NA_integer_,
    FP = NA_integer_,
    FN = NA_integer_,
    FDR = NA_real_,
    Power = NA_real_,
    Tkn = NA_real_,
    MCC = NA_real_,
    F1 = NA_real_
  )

  enableJIT(3)

  for (m in seq_len(M_local)) {
    message(
      "Generator: ", generator,
      " | replicate: ", m, "/", M_local
    )

    set.seed(knockoff_seeds[m])

    # 1. Generate knockoff variables.
    Xk <- generate_knockoffs(
      X = X,
      generator = generator,
      seed = knockoff_seeds[m]
    )

    # 2. Compute the LCD knockoff statistics.
    y <- Y_multi[, m]

    W <- compute_glmnet_W(
      X = X,
      Xk = Xk,
      y = y,
      family = family
    )

    W_matrix[m, ] <- W

    # 3. Compute the knockoff threshold.
    threshold <- compute_knockoff_threshold(
      W = W,
      alpha_kn = alpha_kn,
      offset = offset
    )

    threshold_vector[m] <- threshold

    # 4. Select variables and compute performance metrics.
    selected <- which(W >= threshold)
    selected_list[[m]] <- selected

    metrics <- selection_metrics(
      selected = selected,
      true_signals = true_signals,
      p = p
    )

    per_round[m, c(
      "R", "TP", "FP", "FN",
      "FDR", "Power", "MCC", "F1"
    )] <- unlist(
      metrics[c(
        "R", "TP", "FP", "FN",
        "FDR", "Power", "MCC", "F1"
      )]
    )

    per_round$Tkn[m] <- threshold
  }

  list(
    per_round_metrics = per_round,
    per_round_selected = selected_list,
    generator = generator,
    W = W_matrix,
    Tkn = threshold_vector,
    params = list(
      M = M_local,
      alpha_kn = alpha_kn,
      offset = offset,
      family = family
    )
  )
}

# ------------------------------------------------------------------------------
# 10. LASSO benchmark
# ------------------------------------------------------------------------------

run_lasso_simulation <- function(
    X,
    Y_multi,
    true_signals,
    use_lambda = c("min", "1se"),
    nfolds = 10,
    alpha = 1,
    standardize = TRUE,
    family = c("gaussian", "binomial")
) {
  use_lambda <- match.arg(use_lambda)
  family <- match.arg(family)

  p <- ncol(X)
  M_local <- ncol(Y_multi)
  snp_names <- colnames(X)

  selected_list <- vector("list", M_local)
  lambda_used <- numeric(M_local)
  coefficient_list <- vector("list", M_local)

  per_round <- data.frame(
    round = seq_len(M_local),
    R = NA_integer_,
    TP = NA_integer_,
    FP = NA_integer_,
    FN = NA_integer_,
    FDR = NA_real_,
    Power = NA_real_,
    MCC = NA_real_,
    F1 = NA_real_
  )

  enableJIT(3)

  for (m in seq_len(M_local)) {
    message("LASSO | replicate: ", m, "/", M_local)

    y <- Y_multi[, m]

    # Cross-validation to select lambda.
    cv_fit <- glmnet::cv.glmnet(
      x = X,
      y = y,
      alpha = alpha,
      family = family,
      nfolds = nfolds,
      standardize = standardize,
      intercept = TRUE
    )

    selected_lambda <- if (use_lambda == "min") {
      cv_fit$lambda.min
    } else {
      cv_fit$lambda.1se
    }

    lambda_used[m] <- selected_lambda

    # Fit the LASSO model at the selected lambda.
    fit <- glmnet::glmnet(
      x = X,
      y = y,
      alpha = alpha,
      lambda = selected_lambda,
      family = family,
      standardize = standardize,
      intercept = TRUE
    )

    coefficients <- as.numeric(
      predict(
        fit,
        type = "coefficients",
        s = selected_lambda
      )
    )[-1]

    names(coefficients) <- snp_names

    selected <- which(coefficients != 0)

    selected_list[[m]] <- selected
    coefficient_list[[m]] <- coefficients

    metrics <- selection_metrics(
      selected = selected,
      true_signals = true_signals,
      p = p
    )

    per_round[m, c(
      "R", "TP", "FP", "FN",
      "FDR", "Power", "MCC", "F1"
    )] <- unlist(
      metrics[c(
        "R", "TP", "FP", "FN",
        "FDR", "Power", "MCC", "F1"
      )]
    )
  }

  list(
    per_round_metrics = per_round,
    selected_list = selected_list,
    coefficient_list = coefficient_list,
    lambda_used = lambda_used,
    params = list(
      use_lambda = use_lambda,
      nfolds = nfolds,
      standardize = standardize,
      alpha = alpha,
      family = family
    )
  )
}


# ------------------------------------------------------------------------------
# 11. Run the simulation study
# ------------------------------------------------------------------------------

generator_labels <- c(
  second_order = "SO-k",
  HMM = "HMM-k",
  haldane = "Haldane-k",
  empirical_markov = "EM-k",
  tree = "Tree-k"
)

knockoff_results <- vector(
  "list",
  length(generator_labels)
)
names(knockoff_results) <- names(generator_labels)

knockoff_times <- numeric(length(generator_labels))
names(knockoff_times) <- names(generator_labels)

for (generator_name in names(generator_labels)) {
  start_time <- Sys.time()

  knockoff_results[[generator_name]] <- run_knockoff_simulation(
    X = X,
    Y_multi = Y_multi,
    true_signals = signal_index,
    family = "gaussian",
    alpha_kn = target_fdr,
    offset = 1,
    generator = generator_name
  )

  knockoff_times[generator_name] <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )
}

start_time <- Sys.time()

lasso_result <- run_lasso_simulation(
  X = X,
  Y_multi = Y_multi,
  true_signals = signal_index,
  use_lambda = "min",
  nfolds = 10,
  family = "gaussian"
)

lasso_time <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "secs"
  )
)


# ------------------------------------------------------------------------------
# 12. Summarize performance
# ------------------------------------------------------------------------------

knockoff_metrics_long <- purrr::imap_dfr(
  knockoff_results,
  function(result, generator_name) {
    result$per_round_metrics %>%
      dplyr::select(
        round,
        FDR,
        Power,
        MCC,
        F1
      ) %>%
      tidyr::pivot_longer(
        cols = c(FDR, Power, MCC, F1),
        names_to = "Metric",
        values_to = "Value"
      ) %>%
      dplyr::mutate(
        Method = unname(generator_labels[generator_name])
      )
  }
)

lasso_metrics_long <- lasso_result$per_round_metrics %>%
  dplyr::select(
    round,
    FDR,
    Power,
    MCC,
    F1
  ) %>%
  tidyr::pivot_longer(
    cols = c(FDR, Power, MCC, F1),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(Method = "LASSO")

all_metrics_long <- dplyr::bind_rows(
  knockoff_metrics_long,
  lasso_metrics_long
)

performance_summary <- all_metrics_long %>%
  dplyr::group_by(Method, Metric) %>%
  dplyr::summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

print(performance_summary)

timing_summary <- tibble::tibble(
  Method = c(
    unname(generator_labels),
    "LASSO"
  ),
  Seconds = c(
    unname(knockoff_times),
    lasso_time
  )
)

print(timing_summary)


# ------------------------------------------------------------------------------
# 13. Visualization: FDR and power
# ------------------------------------------------------------------------------

method_order <- c(
  "LASSO",
  "SO-k",
  "HMM-k",
  "Haldane-k",
  "EM-k",
  "Tree-k"
)

plot_data <- all_metrics_long %>%
  dplyr::filter(Metric %in% c("FDR", "Power")) %>%
  dplyr::mutate(
    Method = factor(Method, levels = method_order)
  )

fdr_plot <- plot_data %>%
  dplyr::filter(Metric == "FDR") %>%
  ggplot(aes(x = Method, y = Value)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.35
  ) +
  geom_jitter(
    width = 0.075,
    height = 0,
    alpha = 0.5,
    size = 1.5
  ) +
  geom_hline(
    yintercept = target_fdr,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "(a) FDR"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

power_plot <- plot_data %>%
  dplyr::filter(Metric == "Power") %>%
  ggplot(aes(x = Method, y = Value)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.35
  ) +
  geom_jitter(
    width = 0.075,
    height = 0,
    alpha = 0.5,
    size = 1.5
  ) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "(b) Power"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

combined_plot <- gridExtra::grid.arrange(
  fdr_plot,
  power_plot,
  ncol = 1
)

ggsave(
  filename = file.path(output_dir, "variable_selection_metrics.pdf"),
  plot = combined_plot,
  device = "pdf",
  width = 14,
  height = 14
)
