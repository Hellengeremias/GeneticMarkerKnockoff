# ==============================================================================
# Simulation study: knockoff generators for genetic marker selection
# ==============================================================================
#
# This script contains the simulation workflow used to compare LASSO with
# knockoff-based variable selection under different dependence structures among
# genetic markers.
#
# Main steps:
#   1. Load and preprocess genotype data.
#   2. Select one LD scenario.
#   3. Estimate HMM parameters using fastPHASE.
#   4. Estimate the probability models required by the proposed generators.
#   5. Simulate multiple response vectors.
#   6. Run the standard knockoff filter and the LASSO benchmark.
#   7. Summarize selection metrics.
#
# Notes: set ld_scenario to one of: "low_ld", "moderate_ld", or "high_ld".
# ================================================================================


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
library(rpart)


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

# Offset used by the knockoff+ threshold.
offset_c <- 1

# Window half-width for the classification-tree generator.
tree_window <- 5

# Relatedness threshold used to remove related individuals.
relatedness_threshold <- 0.354

# The original simulation results used V2 to define the independent individuals.
# Keep this as "V2" unless you intentionally want to test the alternative V1 rule.
relatedness_id_column <- "V2"

# Local paths. Environment variables can be used to override these defaults.
isa_data_dir <- Sys.getenv(
  "ISA_DATA_DIR",
  unset = "/home/hellen/ISA/Bases"
)

fastphase_dir <- Sys.getenv(
  "FASTPHASE_DIR",
  unset = "/home/hellen/Programas"
)

fastphase_executable <- Sys.getenv(
  "FASTPHASE_EXECUTABLE",
  unset = file.path(fastphase_dir, "fastPHASE")
)

# Change this manually when running another scenario.
fastphase_tag <- ld_scenario

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

# If FALSE, existing HMM parameter files will be reused when all files are available.
force_refit_hmm <- TRUE

# This recoding is required only for the high-LD mouse dataset when the genotypes
# are encoded as -1, 0, and 1.
apply_high_ld_recoding <- TRUE


# ------------------------------------------------------------------------------
# 3. Helper functions for genotype preprocessing
# ------------------------------------------------------------------------------

calc_maf <- function(x) {
  x <- x[!is.na(x)]
  allele_frequency <- mean(x) / 2
  min(allele_frequency, 1 - allele_frequency)
}


select_100_by_clustering <- function(G, k = 100) {
  
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
# 4. Genotype data and LD scenario
# ------------------------------------------------------------------------------

if (ld_scenario %in% c("low_ld", "moderate_ld")) {

  # SNP map for chromosome 22.
  map <- read.table(
    file.path(isa_data_dir, "isa_rec_rs_aut_hg37_22.map"),
    header = FALSE
  )
  map <- map[order(map[, c("V4")], decreasing = FALSE), ]

  # Genomic relationship matrix.
  grm <- read.table(
    file.path(isa_data_dir, "isa.grm841x841.domicilio.grm"),
    header = FALSE
  )

  # Keep unrelated individuals using the V2 rule from the original results.
  related_pairs <- grm %>%
    dplyr::filter(V4 > relatedness_threshold, V1 != V2)

  diagonal_ids <- grm %>%
    dplyr::filter(V4 > relatedness_threshold, V1 == V2)

  unrelated_ids <- dplyr::anti_join(
    diagonal_ids,
    related_pairs,
    by = relatedness_id_column
  )

  id_vector <- unrelated_ids[[relatedness_id_column]]

  genotype_data <- read.table(
    file.path(isa_data_dir, "SNPs_chr22_total_imput_fam.txt"),
    header = TRUE
  )

  genotype_data <- genotype_data[
    row.names(genotype_data) %in% as.character(id_vector),
    ,
    drop = FALSE
  ]

  # Retain common variants (MAF >= 0.05).
  mafs <- apply(genotype_data, 2, calc_maf)
  genotipos_filtrados <- genotype_data[, mafs >= 0.05, drop = FALSE]

  if (ld_scenario == "low_ld") {
    set.seed(1)
    selected_markers <- select_100_by_clustering(genotipos_filtrados, k = 100)
    genotipos_filtrados_r <- genotipos_filtrados[
      ,
      selected_markers$index,
      drop = FALSE
    ]
  }

  if (ld_scenario == "moderate_ld") {
    genotipos_filtrados_r <- genotipos_filtrados[, 1:100, drop = FALSE]
  }

} else if (ld_scenario == "high_ld") {

  # High-LD genotype data.
  microS <- read.csv(mouse_genotype_file)
  genotipos_filtrados_r <- microS[, -1, drop = FALSE]

  caract_microS <- readxl::read_xlsx(mouse_map_file)
  names(genotipos_filtrados_r)[seq_len(ncol(genotipos_filtrados_r))] <-
    caract_microS$marker_name[seq_len(ncol(genotipos_filtrados_r))]

}


X <- as.matrix(genotipos_filtrados_r)

if (ld_scenario == "high_ld" && isTRUE(apply_high_ld_recoding)) {
  # Original high-LD recoding: -1/0/1 -> 0/1/2.
  X[X == 0] <- 3
  X[X == -1] <- 0
  X[X == 1] <- 2
  X[X == 3] <- 1
}

storage.mode(X) <- "integer"

cat(
  "Scenario:", ld_scenario,
  "\nIndividuals:", nrow(X),
  "\nMarkers:", ncol(X),
  "\nGenotype counts:\n"
)
print(table(as.vector(X), useNA = "ifany"))


# ------------------------------------------------------------------------------
# 5. HMM parameter estimation for the benchmark generator
# ------------------------------------------------------------------------------

fastphase_input <- file.path(
  fastphase_dir,
  paste0("Xfastphase_", fastphase_tag, ".inp")
)

fastphase_output <- file.path(
  fastphase_dir,
  paste0("hmm_param_", fastphase_tag)
)

rhat_file <- paste0(fastphase_output, "_rhat.txt")
alpha_file <- paste0(fastphase_output, "_alphahat.txt")
theta_file <- paste0(fastphase_output, "_thetahat.txt")

hmm_parameter_files <- c(
  rhat_file,
  alpha_file,
  theta_file
)

if (isTRUE(force_refit_hmm) || !all(file.exists(hmm_parameter_files))) {
  writeXtoInp(
    X,
    phased = FALSE,
    out_file = fastphase_input
  )

  runFastPhase(
    fp_path = fastphase_executable,
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


haldane_recombination <- function(dist) {
  # dist: genetic distance in Morgans.
  0.5 * (1 - exp(-2 * dist))
}


prob_equal_genotype <- function(recomb, gen1) {
  if (gen1 == 1) {
    (recomb^2) + (1 - recomb)^2
  } else {
    (1 - recomb)^2
  }
}


prob_one_step_difference <- function(r) {
  2 * (1 - r) * r
}


prob_two_step_difference <- function(r) {
  r^2
}


genotype_transition_prob <- function(g1, g2, r) {
  d <- abs(g1 - g2)

  if (d == 0) {
    prob_equal_genotype(r, g1)
  } else if (d == 1) {
    if (g1 == 1) {
      prob_one_step_difference(r) / 2
    } else {
      prob_one_step_difference(r)
    }
  } else {
    prob_two_step_difference(r)
  }
}


complete_probability_matrix <- function(mat, observed_genotypes) {
  # Complete probability matrices with fewer than three genotype categories.
  gen <- c(0, 1, 2)
  n <- ncol(mat)

  completed_mat <- matrix(0, nrow = length(gen), ncol = n)
  rownames(completed_mat) <- as.character(gen)

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

  prob_list <- vector("list", p)
  X_k <- matrix(NA, nrow = n, ncol = p)

  enableJIT(3)

  for (j in 1:p) {
    probs_j <- matrix(NA, nrow = length(gen), ncol = n)

    for (i in 1:n) {
      probs <- numeric(length(gen))

      if (j == 1) {

        # First marker: condition only on the marker to the right.
        g_right <- X[i, j + 1]
        r1 <- haldane_recombination(
          abs(marker_position[2] - marker_position[1])
        )

        probs <- sapply(
          gen,
          function(g) max(genotype_transition_prob(g_right, g, r1), 1e-6)
        )

      } else if (j == p) {

        # Last marker: condition only on the marker to the left.
        g_left <- X[i, j - 1]
        r1 <- haldane_recombination(
          abs(marker_position[j] - marker_position[j - 1])
        )

        probs <- sapply(
          gen,
          function(g) max(genotype_transition_prob(g_left, g, r1), 1e-6)
        )

      } else if (j > 1 & j < p) {

        # Intermediate markers: condition on the nearest left and right markers.
        g_left <- X[i, j - 1]
        g_right <- X[i, j + 1]

        r1 <- haldane_recombination(
          abs(marker_position[j] - marker_position[j - 1])
        )
        r2 <- haldane_recombination(
          abs(marker_position[j + 1] - marker_position[j])
        )
        r12 <- haldane_recombination(
          abs(marker_position[j + 1] - marker_position[j - 1])
        )

        denom <- max(genotype_transition_prob(g_left, g_right, r12), 1e-6)

        probs <- sapply(gen, function(g) {
          p1 <- max(genotype_transition_prob(g_left, g, r1), 1e-6)
          p2 <- max(genotype_transition_prob(g, g_right, r2), 1e-6)
          (p1 * p2) / denom
        })
      }

      probs <- probs / sum(probs)
      X_k[i, j] <- gen[sample_discrete(probs)]
      probs_j[, i] <- probs
    }

    prob_list[[j]] <- probs_j
  }

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
  gen <- c(0, 1, 2)

  prob_list <- vector("list", p)
  X_k <- matrix(NA, nrow = n, ncol = p)

  enableJIT(3)

  for (j in 1:p) {
    probs_j <- matrix(NA, nrow = length(gen), ncol = n)

    for (i in 1:n) {
      probs <- numeric(length(gen))

      for (v in seq_along(gen)) {
        x_current <- gen[v]

        if (j == 1) {

          x_right <- X[i, j + 1]
          conj <- sum(X[, j] == x_current & X[, j + 1] == x_right)
          marg <- sum(X[, j + 1] == x_right)
          probs[v] <- ifelse(marg > 0, conj / marg, 0)

        } else if (j == p) {

          x_left <- X[i, j - 1]
          conj <- sum(X[, j] == x_current & X[, j - 1] == x_left)
          marg <- sum(X[, j - 1] == x_left)
          probs[v] <- ifelse(marg > 0, conj / marg, 0)

        } else {

          x_left <- X[i, j - 1]
          x_right <- X[i, j + 1]

          conj1 <- sum(X[, j] == x_current & X[, j - 1] == x_left)
          marg1 <- sum(X[, j - 1] == x_left)
          p_current_given_left <- ifelse(marg1 > 0, conj1 / marg1, 0)

          conj2 <- sum(X[, j + 1] == x_right & X[, j] == x_current)
          marg2 <- sum(X[, j] == x_current)
          p_right_given_current <- ifelse(marg2 > 0, conj2 / marg2, 0)

          conj3 <- sum(X[, j + 1] == x_right & X[, j - 1] == x_left)
          marg3 <- sum(X[, j - 1] == x_left)
          p_right_given_left <- ifelse(marg3 > 0, conj3 / marg3, 1e-10)

          probs[v] <- p_current_given_left *
            p_right_given_current / p_right_given_left
        }
      }

      probs <- probs / sum(probs)
      X_k[i, j] <- gen[sample_discrete(probs)]
      probs_j[, i] <- probs
    }

    prob_list[[j]] <- probs_j
  }

  SNPs <- matrix(NA, nrow = n, ncol = 2 * p)
  SNPs[, seq(1, 2 * p, by = 2)] <- X
  SNPs[, seq(2, 2 * p, by = 2)] <- X_k

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
  gen <- c(0, 1, 2)

  prob_list <- vector("list", p)
  X_k <- matrix(NA, nrow = n, ncol = p)

  enableJIT(3)

  for (j in 1:p) {

    # Fit a classification tree for the j-th marker using neighboring markers
    # within a window of total size 2*w.
    target <- X[, j]

    if (j == 1) {
      # First marker: use up to 2*w markers to the right.
      end_idx <- min(j + 2 * w, p)
      predictors <- as.data.frame(X[, (j + 1):end_idx])

    } else if (j == p) {
      # Last marker: use up to 2*w markers to the left.
      start_idx <- max(j - 2 * w, 1)
      predictors <- as.data.frame(X[, start_idx:(j - 1)])

    } else {
      # Intermediate markers: use up to w markers on each side.
      start_idx <- max(j - w, 1)
      X_left <- X[, start_idx:(j - 1)]

      end_idx <- min(j + w, p)
      X_right <- X[, (j + 1):end_idx]
      predictors <- as.data.frame(cbind(X_left, X_right))
    }

    tree_formula <- as.formula("target ~ .")
    tree_fit <- rpart::rpart(
      tree_formula,
      data = cbind(target, predictors),
      method = "class"
    )

    probs_rpart <- predict(
      tree_fit,
      newdata = predictors,
      type = "prob"
    )

    probs_rpart_t <- t(probs_rpart)

    if (nrow(probs_rpart_t) < 3) {
      observed_genotypes <- rownames(probs_rpart_t)
      probs_complete <- complete_probability_matrix(
        mat = probs_rpart_t,
        observed_genotypes = observed_genotypes
      )
    } else {
      probs_complete <- probs_rpart_t
    }

    probs_matrix <- t(matrix(probs_complete, ncol = n, nrow = 3))

    row_totals <- rowSums(probs_matrix)
    probs_final <- probs_matrix / row_totals

    for (i in 1:n) {
      X_k[i, j] <- gen[sample_discrete(probs_final[i, ])]
    }

    prob_list[[j]] <- t(probs_final)
  }

  SNPs <- matrix(NA, nrow = n, ncol = 2 * p)
  SNPs[, seq(1, 2 * p, by = 2)] <- X
  SNPs[, seq(2, 2 * p, by = 2)] <- X_k

  return(list(
    SNPs = SNPs,
    prob_list = prob_list
  ))
}


# ------------------------------------------------------------------------------
# 7. Precompute genotype probabilities for the proposed generators
# ------------------------------------------------------------------------------

if (ld_scenario == "high_ld") {
  marker_position <- caract_microS$original_cM[seq_len(ncol(X))] / 100
} else {
  # Positions are obtained from the map using the marker names in the selected genotype matrix.
  snp_names <- sub("_.*", "", colnames(genotipos_filtrados_r))
  map_r <- map %>% dplyr::filter(V2 %in% snp_names)
  marker_position <- map_r$V3
}


set.seed(1)
haldane_model <- estimate_haldane_probabilities(
  X = X,
  marker_position = marker_position
)
Xh <- haldane_model$SNPs[, seq(2, 2 * ncol(X), by = 2)]
prob_haldane <- haldane_model$prob_list

set.seed(1)
empirical_markov_model <- estimate_empirical_markov_probabilities(X = X)
Xv <- empirical_markov_model$SNPs[, seq(2, 2 * ncol(X), by = 2)]
prob_empirical_markov <- empirical_markov_model$prob_list

set.seed(1)
tree_model <- estimate_tree_probabilities(X = X, w = tree_window)
Xr <- tree_model$SNPs[, seq(2, 2 * ncol(X), by = 2)]
prob_tree <- tree_model$prob_list

stopifnot(
  !anyNA(Xh),
  !anyNA(Xv),
  !anyNA(Xr)
)


# ------------------------------------------------------------------------------
# 8. Response simulation
# ------------------------------------------------------------------------------

set.seed(1)
signal_index <- sample(1:ncol(X), n_signals)

beta <- rep(0, ncol(X))

set.seed(1001)
beta[signal_index] <- rnorm(n_signals, mean = 0.75, sd = 0.5)

beta_signal <- beta[signal_index[order(signal_index)]]
names(beta_signal) <- signal_index[order(signal_index)]

beta_signal_df <- data.frame(
  Variable = signal_index,
  beta_original = beta[signal_index]
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
# 9. Standard knockoff filter simulation
# ------------------------------------------------------------------------------

# Compute W using glmnet coefficient differences:
# W_j = |beta_j| - |beta_tilde_j|.
compute_glmnet_W <- function(X, Xk, y, family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  W <- knockoff::stat.glmnet_coefdiff(X, Xk, y = y, family = family)
  return(as.numeric(W))
}

# Compute the knockoff+ threshold.
compute_knockoff_threshold <- function(W, alpha_kn = 0.20, c = 1) {
  threshold_candidates <- sort(unique(W[W > 0]), decreasing = FALSE)

  if (length(threshold_candidates) == 0) {
    return(Inf)
  }

  for (threshold in threshold_candidates) {
    n_negative <- sum(W <= -threshold)
    n_positive <- sum(W >= threshold)

    if ((c + n_negative) / max(1, n_positive) <= alpha_kn) {
      return(threshold)
    }
  }

  return(Inf)
}


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
  gen <- c(0, 1, 2)

  X_k <- matrix(NA, nrow = n, ncol = p)

  enableJIT(3)

  for (j in 1:p) {
    probs_j <- prob[[j]]

    for (i in 1:n) {
      X_k[i, j] <- gen[sample_discrete(probs_j[, i])]
    }
  }

  SNPs <- matrix(NA, nrow = n, ncol = 2 * p)
  SNPs[, seq(1, 2 * p, by = 2)] <- X
  SNPs[, seq(2, 2 * p, by = 2)] <- X_k

  return(SNPs)
}


# Generate knockoff variables using one of the evaluated generators.
generate_knockoffs <- function(
    X,
    generator = c("second_order", "HMM", "haldane", "vizinhos", "rpart"),
    seed = NULL
) {
  generator <- match.arg(generator)

  if (generator == "second_order") {
    Xk <- knockoff::create.second_order(X)
    return(Xk)
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
    }

    if (generator == "haldane") {
      set.seed(current_seed)
      SNPs <- generate_from_probabilities(X = X, prob = prob_haldane)
      Xk <- SNPs[, seq(2, 2 * ncol(X), by = 2)]
    }

    if (generator == "vizinhos") {
      set.seed(current_seed)
      SNPs <- generate_from_probabilities(X = X, prob = prob_empirical_markov)
      Xk <- SNPs[, seq(2, 2 * ncol(X), by = 2)]
    }

    if (generator == "rpart") {
      set.seed(current_seed)
      SNPs <- generate_from_probabilities(X = X, prob = prob_tree)
      Xk <- SNPs[, seq(2, 2 * ncol(X), by = 2)]
    }

    if (!has_constant_column(Xk)) {
      return(Xk)
    }

    current_seed <- sample.int(1e8, 1L)
  }
}


selection_metrics <- function(S, H1, p) {
  S <- sort(unique(S))
  R <- length(S)
  s <- length(H1)

  TP <- length(intersect(S, H1))
  FP <- R - TP
  FN <- s - TP
  TN <- (p - s) - FP

  FDR <- if (R == 0) 0 else FP / R
  Power <- if (s == 0) NA_real_ else TP / s

  # Convert to double before multiplication to avoid integer overflow
  mcc_terms <- as.numeric(c(TP + FP, TP + FN, TN + FP, TN + FN))
  mcc_denominator <- sqrt(prod(mcc_terms))

  MCC <- if (!is.finite(mcc_denominator) || mcc_denominator == 0) {
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


run_knockoff_simulation <- function(
    X,
    Y_multi,
    true_signals = NULL,
    family = c("gaussian", "binomial"),
    alpha_kn = 0.20,
    offset_c = 1,
    generator = c("second_order", "HMM", "haldane", "vizinhos", "rpart"),
    importance = c("glmnet_coefdiff")
) {
  family <- match.arg(family)
  generator <- match.arg(generator)
  importance <- match.arg(importance)

  p <- ncol(X)
  M_local <- ncol(Y_multi)

  seeds <- sample.int(1e8, M_local)

  W_list <- matrix(NA_real_, nrow = M_local, ncol = p)
  T_vec <- rep(NA_real_, M_local)
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
    F1 = NA_real_,
    stringsAsFactors = FALSE
  )

  enableJIT(3)

  for (m in seq_len(M_local)) {

    set.seed(seeds[m])
    message("Generator: ", generator, " | replicate: ", m, "/", M_local)

    Xk <- generate_knockoffs(
      X = X,
      generator = generator,
      seed = seeds[m]
    )

    y <- Y_multi[, m]
    W <- compute_glmnet_W(X, Xk, y, family = family)
    W_list[m, ] <- W

    Tm <- compute_knockoff_threshold(
      W,
      alpha_kn = alpha_kn,
      c = offset_c
    )
    T_vec[m] <- Tm

    selected <- which(W >= Tm)
    selected_list[[m]] <- selected

    metrics <- selection_metrics(selected, true_signals, p)

    per_round$R[m] <- metrics$R
    per_round$TP[m] <- metrics$TP
    per_round$FP[m] <- metrics$FP
    per_round$FN[m] <- metrics$FN
    per_round$FDR[m] <- metrics$FDR
    per_round$Power[m] <- metrics$Power
    per_round$MCC[m] <- metrics$MCC
    per_round$F1[m] <- metrics$F1
    per_round$Tkn[m] <- Tm
  }

  list(
    per_round_metrics = per_round,
    per_round_selected = selected_list,
    generator = generator,
    W = W_list,
    Tkn = T_vec,
    params = list(
      M = M_local,
      alpha_kn = alpha_kn,
      offset_c = offset_c,
      family = family,
      importance = importance
    )
  )
}


# ------------------------------------------------------------------------------
# 10. LASSO benchmark
# ------------------------------------------------------------------------------

run_lasso_simulation <- function(
    X,
    Y_multi,
    true_signals = NULL,
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
    F1 = NA_real_,
    stringsAsFactors = FALSE
  )

  enableJIT(3)

  for (m in seq_len(M_local)) {
    message("LASSO | replicate: ", m, "/", M_local)

    y <- Y_multi[, m]

    cv_fit <- glmnet::cv.glmnet(
      X,
      y,
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

    fit <- glmnet::glmnet(
      X,
      y,
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

    metrics <- selection_metrics(selected, true_signals, p)

    per_round$R[m] <- metrics$R
    per_round$TP[m] <- metrics$TP
    per_round$FP[m] <- metrics$FP
    per_round$FN[m] <- metrics$FN
    per_round$FDR[m] <- metrics$FDR
    per_round$Power[m] <- metrics$Power
    per_round$MCC[m] <- metrics$MCC
    per_round$F1[m] <- metrics$F1
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

generator_order <- c(
  "second_order",
  "HMM",
  "haldane",
  "vizinhos",
  "rpart"
)

generator_labels <- c(
  second_order = "SO-k",
  HMM = "HMM-k",
  haldane = "Haldane-k",
  vizinhos = "EM-k",
  rpart = "Tree-k"
)

knockoff_results <- vector("list", length(generator_order))
names(knockoff_results) <- generator_order

knockoff_times <- numeric(length(generator_order))
names(knockoff_times) <- generator_order

for (generator_name in generator_order) {
  start_time <- Sys.time()

  knockoff_results[[generator_name]] <- run_knockoff_simulation(
    X = X,
    Y_multi = Y_multi,
    true_signals = signal_index,
    family = "gaussian",
    alpha_kn = target_fdr,
    offset_c = offset_c,
    generator = generator_name
  )

  knockoff_times[generator_name] <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
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
  difftime(Sys.time(), start_time, units = "secs")
)


# ------------------------------------------------------------------------------
# 12. Summarize performance
# ------------------------------------------------------------------------------

knockoff_metrics_long <- purrr::imap_dfr(
  knockoff_results,
  function(result, generator_name) {
    result$per_round_metrics %>%
      dplyr::select(round, FDR, Power, MCC, F1) %>%
      tidyr::pivot_longer(
        cols = c(FDR, Power, MCC, F1),
        names_to = "Metric",
        values_to = "Value"
      ) %>%
      dplyr::mutate(Method = unname(generator_labels[generator_name]))
  }
)

lasso_metrics_long <- lasso_result$per_round_metrics %>%
  dplyr::select(round, FDR, Power, MCC, F1) %>%
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
  Method = c(unname(generator_labels[generator_order]), "LASSO"),
  Seconds = c(unname(knockoff_times[generator_order]), lasso_time)
)

print(timing_summary)

saveRDS(
  list(
    ld_scenario = ld_scenario,
    X_dim = dim(X),
    signal_index = signal_index,
    beta = beta,
    beta_signal_df = beta_signal_df,
    performance_summary = performance_summary,
    timing_summary = timing_summary,
    knockoff_results = knockoff_results,
    lasso_result = lasso_result,
    all_metrics_long = all_metrics_long
  ),
  file = file.path(output_dir, paste0("simulation_results_", ld_scenario, ".rds"))
)

readr::write_csv(
  performance_summary,
  file.path(output_dir, paste0("performance_summary_", ld_scenario, ".csv"))
)


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
  dplyr::mutate(Method = factor(Method, levels = method_order))

fdr_plot <- plot_data %>%
  dplyr::filter(Metric == "FDR") %>%
  ggplot(aes(x = Method, y = Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.35) +
  geom_jitter(width = 0.075, height = 0, alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = target_fdr, linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = NULL, y = NULL, title = "(a) FDR") +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

power_plot <- plot_data %>%
  dplyr::filter(Metric == "Power") %>%
  ggplot(aes(x = Method, y = Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.35) +
  geom_jitter(width = 0.075, height = 0, alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = NULL, y = NULL, title = "(b) Power") +
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
  filename = file.path(
    output_dir,
    paste0("variable_selection_metrics_", ld_scenario, ".pdf")
  ),
  plot = combined_plot,
  device = "pdf",
  width = 14,
  height = 14
)
