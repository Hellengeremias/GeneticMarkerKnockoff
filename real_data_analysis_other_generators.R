# ==============================================================================
# Real-data analysis: Model-X and proposed knockoff generators
# ==============================================================================
#
# This script applies the second-order Model-X generator and the proposed
# Haldane, empirical Markov, and classification-tree knockoff generators to the
# ISA-Nutrição 2015 genomic data.
#
# Main steps:
#   1. Prepare phenotype, covariates, and chromosome-specific genotype data.
#   2. Divide each chromosome into blocks.
#   3. Estimate block-specific probability models for the proposed generators.
#   4. Generate knockoff variables repeatedly within each block.
#   5. Compute LASSO-based knockoff statistics.
#   6. Apply the knockoff+ threshold within blocks and summarize selection
#      frequencies across repetitions.
#
# NOTE:
# - Controlled-access genotype data are not distributed with this repository.
# - Update the local file paths in the analysis workflow before running.
# ==============================================================================

#-------------------------------------------------------------------------------
# 1. Packages
#-------------------------------------------------------------------------------
library(haven)
library(dplyr)
library(readr)
library(knockoff)
library(compiler)

#-------------------------------------------------------------------------------
# 2. Timing utilities
#-------------------------------------------------------------------------------
format_elapsed_time <- function(seconds) {
  h <- floor(seconds / 3600)
  m <- floor((seconds %% 3600) / 60)
  s <- seconds %% 60
  sprintf("%02dh %02dm %06.3fs", h, m, s)
}

time_it <- function(expr, label = NULL) {
  t0 <- Sys.time()
  value <- force(expr)
  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  
  list(
    value = value,
    elapsed_sec = elapsed,
    elapsed_fmt = format_elapsed_time(elapsed),
    start_time = t0,
    end_time = t1
  )
}

#-------------------------------------------------------------------------------
# 3. General helper functions
#-------------------------------------------------------------------------------
`%notin%` <- function(x,y) !(x %in% y)

# Minor allele frequency
calc_maf <- function(x) {
  x <- x[!is.na(x)]
  freq <- mean(x) / 2
  min(freq, 1 - freq)
}

# Select the glmnet family from the response type
get_family_type <- function(y) {
  if (length(unique(y)) == 2L) "binomial" else "gaussian"
}

# Check whether any knockoff column is constant
has_constant_column <- function(M) {
  any(apply(M, 2, function(z) var(z, na.rm = TRUE) == 0))
}

# Count how many times each variable was selected across repetitions
count_selected_rounds <- function(per_round_selected, p) {
  counts <- integer(p)
  for (sel in per_round_selected) {
    if (length(sel) > 0) counts[sel] <- counts[sel] + 1L
  }
  counts
}

sample_discrete<-function(p){
  u <- runif(1)
  P <- cumsum(p)
  val <- sum(P<u)+1
  return(val)
}

haldane_recombination <- function(dist) 0.5 * (1 - exp(-2 * dist)) # dist: distance in Morgan

prob_equal_genotype <- function(recomb, gen1) if (gen1 == 1) (recomb^2)+(1-recomb)^2 else (1-recomb)^2
prob_one_step_difference <- function(r) 2 * (1 - r) * r
prob_two_step_difference <- function(r) r^2

genotype_transition_prob <- function(g1, g2, r) {
  d <- abs(g1 - g2)
  if (d == 0) prob_equal_genotype(r, g1)
  else if (d == 1) if (g1 == 1) prob_one_step_difference(r)/2 else prob_one_step_difference(r)
  else prob_two_step_difference(r)
}

# Complete the probability matrix when fewer than three genotype classes are observed
complete_probability_matrix <- function(mat, obs_gen) {
  
  # Possible genotype values
  gen <- c(0, 1, 2)
  
  # Number of columns
  n <- ncol(mat)
  
  # Initialize the completed probability matrix
  mat_correct <- matrix(0, nrow = length(gen), ncol = n)
  rownames(mat_correct) <- as.character(gen)
  
  # Copy probabilities for observed genotype classes
  rownames(mat) <- as.character(obs_gen)
  lines_identical <- intersect(rownames(mat), rownames(mat_correct))
  mat_correct[lines_identical, ] <- mat[lines_identical, ]
  
  return(mat_correct)
}

#-------------------------------------------------------------------------------
# 4. Prepare phenotype and covariates (ISA-Nutrição 2015)
#-------------------------------------------------------------------------------
prepare_analysis_inputs <- function(
    path_fenotipos = "/home/hellen/ISA/banco_minicurso_841.csv",
    path_grm = "/home/hellen/ISA/isa.grm841x841.domicilio.grm",
    path_comp = "/home/hellen/ISA/bancoISA_GWAS_n841_2024.dta",
    path_comp2 = "/home/hellen/ISA/CPancestralidadeISA_2024.dta",
    feno = 3
) {
  
  fenotipos <- read.csv(path_fenotipos, header = TRUE)
  grm_v <- read.table(path_grm, sep = "\t")
  comp <- data.frame(read_dta(path_comp)[, c(1, 3, 5, 19, 24, 26, 28, 30, 32, 33)])
  comp2 <- data.frame(read_dta(path_comp2)[, c(1, 4, 5)])
  
  comp <- comp[order(comp[, 1]), ]
  comp2 <- comp2[order(comp2[, 1]), ]
  comp_tot_raw <- cbind(comp[, -1], comp2[, -1])
  
  find_match_one <- function(i) {
    idx <- integer(0)
    
    if (!is.na(fenotipos[i, 13]) && !is.na(fenotipos[i, 19])) {
      idx <- which(
        round(fenotipos[i, 13], 2) == round(comp_tot_raw[, 4], 2) &
          round(fenotipos[i, 19], 4) == round(comp_tot_raw[, 7], 4) &
          fenotipos[i, 3] == comp_tot_raw[, 2]
      )
    } else {
      if (!is.na(fenotipos[i, 13])) {
        idx <- which(
          round(fenotipos[i, 13], 2) == round(comp_tot_raw[, 4], 2) &
            is.na(comp_tot_raw[, 7]) &
            fenotipos[i, 3] == comp_tot_raw[, 2]
        )
      } else if (!is.na(fenotipos[i, 19])) {
        idx <- which(
          round(fenotipos[i, 19], 4) == round(comp_tot_raw[, 7], 4) &
            is.na(comp_tot_raw[, 4]) &
            fenotipos[i, 3] == comp_tot_raw[, 2]
        )
      }
    }
    
    if (length(idx) == 1) return(idx)
    if (length(idx) > 1) return(idx[1])
    NA_integer_
  }
  match_idx <- sapply(seq_len(nrow(fenotipos)), find_match_one)
  
  comp_tot <- matrix(NA_real_, nrow = nrow(fenotipos), ncol = ncol(comp_tot_raw))
  colnames(comp_tot) <- colnames(comp_tot_raw)
  ok <- which(!is.na(match_idx))
  comp_tot[ok, ] <- as.matrix(comp_tot_raw[match_idx[ok], ])
  excl <- which(is.na(comp_tot[, feno]))
  
  dados_grm_r0 <- grm_v %>% dplyr::filter(V4 > .354 & V1 != V2)
  dados_grm_r1 <- grm_v %>% dplyr::filter(V4 > .354 & V1 == V2)
  dados_grm_r2 <- dplyr::anti_join(dados_grm_r1, dados_grm_r0, by = "V2")
  dados_grm_r2 <- dados_grm_r2 %>% dplyr::filter(V2 %notin% excl)
  selec <- dados_grm_r2$V2
  
  y_raw <- comp_tot[, feno]
  y_ind <- y_raw[selec]
  cov_amb_ind <- data.frame(comp_tot[selec, c(1, 2, 10, 11), drop = FALSE])
  
  dadost <- data.frame(
    y = y_ind,
    sexo = cov_amb_ind[, 1],
    idade = cov_amb_ind[, 2],
    CP1 = cov_amb_ind[, 3],
    CP2 = cov_amb_ind[, 4]
  )
  
  list(
    selec = selec,
    dadost = dadost,
    match_idx = match_idx,
    fenotipos = fenotipos,
    comp_tot = comp_tot
  )
}

#-------------------------------------------------------------------------------
# 5. Prepare genotype data for one chromosome
#-------------------------------------------------------------------------------
load_and_prepare_chr <- function(crom,
                                 snp_dir,
                                 map_dir,
                                 selec,
                                 maf_cutoff = 0.05) {
  
  base_snp <- file.path(snp_dir, paste0("SNPs_chr", crom, "_total_imput_fam.txt"))
  base_map <- file.path(map_dir, paste0("isa_rec_rs_aut_hg37_", crom, ".map"))
  
  SNPs <- read.table(file = base_snp, header = TRUE)
  map <- read.table(file = base_map, header = FALSE)
  
  # Order markers by genetic position; genotype columns are already in the corresponding order
  ord <- order(map$V4)
  map <- map[ord, , drop = FALSE]
  
  # select common SNPs (MAF >= 0.05)
  mafs <- apply(SNPs, 2, calc_maf)
  keep <- which(mafs >= maf_cutoff)
  SNPs <- SNPs[, keep, drop = FALSE]
  map <- map[keep, , drop = FALSE]
  mafs <- mafs[keep]
  
  # select individuals according to phenotype vector (y)
  X <- as.matrix(SNPs[selec, , drop = FALSE])
  
  list(
    chromosome = crom,
    X = X,
    map = map,
    maf = mafs,
    snp_names = as.character(map$V2)
  )
}

#-------------------------------------------------------------------------------
# 6. Proposed knockoff generators for genetic markers
#-------------------------------------------------------------------------------
estimate_haldane_probabilities <- function(X, loc_marc) {
  # marker_position: genetic positions of the markers
  
  n <- nrow(X)
  p <- ncol(X)
  gen <- c(0, 1, 2)
  
  # List of genotype-probability matrices
  lista_probs <- vector("list", p)
  
  # Knockoff matrix
  X_k <- matrix(NA, nrow = n, ncol = p)  
  
  # Generate knockoff genotypes
  enableJIT(3)
  
  # Iterate over markers to generate their knockoff counterparts
  for (j in 1:p) {
    probs_j <- matrix(NA, nrow = length(gen), ncol = n)
    
    for (i in 1:n) {
      probs <- numeric(length(gen))
      
      if(j == 1){
        
        # first column: use only the complete column on the right
        g_dir <- X[i, j + 1]
        r1 <- haldane_recombination(abs(loc_marc[2] - loc_marc[1]))
        probs <- sapply(gen, function(g) max(genotype_transition_prob(g_dir, g, r1), 1e-6))
        
      } else if (j == p) {
        # last column: use only the complete column on the left
        g_esq <- X[i, j-1]
        r1 <- haldane_recombination(abs(loc_marc[j] - loc_marc[j-1]))
        probs <- sapply(gen, function(g) max(genotype_transition_prob(g_esq, g, r1), 1e-6))
        
      } else if (j > 1 & j < p) {
        # intermediate columns
        
        # column on the left
        g_esq <- X[i, j-1]
        
        # column on the right
        g_dir <- X[i, j + 1]
        
        # apply haldane function to compute recombination rate
        r1 <- haldane_recombination(abs(loc_marc[j] - loc_marc[j-1]))
        r2 <- haldane_recombination(abs(loc_marc[j + 1] - loc_marc[j]))
        r12 <- haldane_recombination(abs(loc_marc[j + 1] - loc_marc[j-1]))
        
        denom <- max(genotype_transition_prob(g_esq, g_dir, r12), 1e-6)
        
        # Estimate genotype probabilities from the recombination rate
        probs <- sapply(gen, function(g) {
          p1 <- max(genotype_transition_prob(g_esq, g, r1), 1e-6)
          p2 <- max(genotype_transition_prob(g, g_dir, r2), 1e-6)
          (p1 * p2) / denom
        })
      } 
      
      # Normalize probabilities and sample the knockoff genotype
      probs <- probs / sum(probs)
      X_k[i, j] <-  gen[sample_discrete(probs)]
      probs_j[, i] <- probs
    }
    lista_probs[[j]] <- probs_j
  }
  
  # matrix alternating original columns and knockoffs
  SNPs <- matrix(NA, nrow = n, ncol = 2*p)
  SNPs[, seq(1, 2*p, by = 2)] <- X
  SNPs[, seq(2, 2*p, by = 2)] <- X_k
  
  return(list(
    SNPs = SNPs,
    lista_probs = lista_probs
  ))
}

estimate_empirical_markov_probabilities <- function(X) {
  
  p <- ncol(X)
  n <- nrow(X)
  
  # Possible genotype values
  gen <- c(0,1,2)
  
  # List of genotype-probability matrices
  lista_probs <- vector("list", p)
  
  # Knockoff matrix
  X_k <- matrix(NA, nrow = n, ncol = p)  
  
  # Generate knockoff genotypes
  enableJIT(3)
  
  # Iterate over markers to generate their knockoff counterparts
  for (j in 1:p) {
    probs_j <- matrix(NA, nrow = length(gen), ncol = n)
    
    for (i in 1:n) {
      probs <- numeric(length(gen))
      
      for (v in seq_along(gen)) {
        
        # Possible genotype values for original variable v
        x_o <- gen[v]
        
        if (j == 1) {
          # first original variable
          x_d <- X[i, j + 1]
          conj <- sum(X[, j] == x_o & X[, j + 1] == x_d) # P(Xd = x_d, Xf = x_o)
          marg <- sum(X[, j + 1] == x_o)                 # P(Xd = x_o)
          probs[v] <- ifelse(marg > 0, conj / marg, 0)
          
        } else if (j == p) {
          # last original variable
          x_e <- X[i, j - 1]
          conj <- sum(X[, j] == x_o & X[, j - 1] == x_e) # P(Xe = x_e, Xf = x_o)
          marg <- sum(X[, j - 1] == x_e)                 # P(Xe = x_e)
          probs[v] <- ifelse(marg > 0, conj / marg, 0)
          
        } else {
          # intermediate original variables (j > 1 & j < p)
          x_e <- X[i, j - 1]
          x_d <- X[i, j + 1]
          
          conj1 <- sum(X[, j] == x_o & X[, j - 1] == x_e)  # P(Xf = x_o, Xe = x_e)
          marg1 <- sum(X[, j-1] == x_e)                    # P(Xe = x_e)
          p_xf_dado_xe <- ifelse(marg1 > 0, conj1 / marg1, 0)
          
          conj2 <- sum(X[, j + 1] == x_d & X[, j] == x_o)  # P(Xd = x_d, Xf = x_o)
          marg2 <- sum(X[, j] == x_o)                      # P(Xe = x_e)
          p_xd_dado_xf <- ifelse(marg2 > 0, conj2 / marg2, 0)
          
          conj3 <- sum(X[, j + 1] == x_d & X[, j - 1] == x_e) # P(Xd = x_d , Xe = x_e)
          marg3 <- sum(X[, j - 1] == x_e)                     # P(Xe = X_e)
          p_xd_dado_xe <- ifelse(marg3 > 0, conj3 / marg3, 1e-10) 
          
          # Estimate genotype probabilities from empirical Markov transitions
          probs[v] <- p_xf_dado_xe * p_xd_dado_xf / p_xd_dado_xe
        }
      }
      
      # Normalize probabilities and sample the knockoff genotype
      probs <- probs / sum(probs)
      X_k[i, j] <- gen[sample_discrete(probs)]
      probs_j[, i] <- probs
    }
    lista_probs[[j]] <- probs_j
  }
  
  # matrix alternating original columns and knockoffs
  SNPs <- matrix(NA, nrow = n, ncol = 2*p)
  SNPs[, seq(1, 2*p, by = 2)] <- X
  SNPs[, seq(2, 2*p, by = 2)] <- X_k
  
  return(list(
    SNPs = SNPs,
    lista_probs = lista_probs
  ))
}

estimate_tree_probabilities <- function(X, w) {
  
  p <- ncol(X)
  n <- nrow(X)
  
  # Possible genotype values
  gen <- c(0,1,2)
  
  # List of genotype-probability matrices
  lista_probs <- vector("list", p)
  
  # Knockoff matrix
  X_k <- matrix(NA, nrow = n, ncol = p)
  
  
  # Generate knockoff genotypes
  enableJIT(3)
  
  # Iterate over markers to generate their knockoff counterparts
  for (j in 1:p) {
    
    # Fit a classification tree for marker j as the response
    alvo <- X[,j]
    
    # Predictors are neighboring original markers within a window of total width up to 2*w.
    # w is the maximum number of markers used on each side of the target marker.
    if (j == 1) {
      
      # First marker: use up to 2*w predictors to the right
      fim <- min(j + 2*w, p)
      preditoras <- as.data.frame(X[, (j+1):fim])
    } else if (j == p) {
      
      # Last marker: use up to 2*w predictors to the left
      inicio <- max(j - 2*w, 1)
      preditoras <- as.data.frame(X[, inicio:(j-1)])
    } else {
      
      # Intermediate markers: use up to w predictors on each side
      inicio <- max(j - w, 1)
      X_e <- X[, inicio:(j-1)]
      
      fim <- min(j + w, p)
      X_d <- X[, (j+1):fim]
      preditoras <- as.data.frame(cbind(X_e, X_d))
    }
    
    # model adjustment
    library(rpart)
    formula <- as.formula("alvo ~ .")
    arvore <- rpart::rpart(formula, data = cbind(alvo, preditoras),
                           method = "class")
    
    # Estimate genotype probabilities
    probs_rpart <- predict(arvore, newdata = preditoras, type = "prob")
    probs_rpart_t <- t(probs_rpart)
    
    # Complete the probability matrix when fewer than three genotype classes are observed
    if(nrow(probs_rpart_t) <3){
      rowN <- rownames(probs_rpart_t)
      l <- complete_probability_matrix(mat = probs_rpart_t, obs_gen = rowN)
    }else{
      l <- probs_rpart_t
    }
    probs_rpart_t2 <- t(matrix(l, ncol = n, nrow = 3))
    
    # Normalize probabilities and sample the knockoff genotype
    soma <- rowSums(probs_rpart_t2)
    probs_f <- probs_rpart_t2 / soma
    for(i in 1:n){
      X_k[i, j] <- gen[sample_discrete(probs_f[i,])]  
    }
    lista_probs[[j]] <- t(probs_f)
  }
  
  # matrix alternating original columns and knockoffs
  SNPs <- matrix(NA, nrow = n, ncol = 2*p)
  SNPs[, seq(1, 2*p, by = 2)] <- X
  SNPs[, seq(2, 2*p, by = 2)] <- X_k
  
  return(list(
    SNPs = SNPs,
    lista_probs = lista_probs
  ))
}

# Generate knockoffs from the probability models estimated by the proposed methods
generate_from_probabilities <- function(X, prob) {
  p <- ncol(X)
  n <- nrow(X)
  gen <- c(0, 1, 2)
  
  X_k <- matrix(NA, nrow = n, ncol = p)
  
  compiler::enableJIT(3)
  for (j in seq_len(p)) {
    probs_j <- prob[[j]]
    for (i in seq_len(n)) {
      X_k[i, j] <- gen[sample_discrete(probs_j[, i])]
    }
  }
  
  X_k
}

#-------------------------------------------------------------------------------
# 7. Define blocks within chromosomes
#-------------------------------------------------------------------------------
make_blocks <- function(p, block_size = 1000) {
  starts <- seq(1L, p, by = block_size)
  ends <- pmin(starts + block_size - 1L, p)
  
  blocks <- vector("list", length(starts))
  
  for (b in seq_along(starts)) {
    blocks[[b]] <- list(
      block_id = b,
      idx = starts[b]:ends[b]
    )
  }
  
  blocks
}

#-------------------------------------------------------------------------------
# 8. Prepare generator state for block-wise analysis
#-------------------------------------------------------------------------------
make_method_state <- function(method, prob = NULL) {
  method <- match.arg(method, c("second_order", "haldane", "EM", "tree"))
  
  list(
    method = method,
    prob = prob
  )
}

method_state_builder_prob_block <- function(method,
                                            X_block,
                                            map_block = NULL,
                                            loc_col = "V3",
                                            w_rpart = 5,
                                            ...) {
  
  method <- match.arg(method, c("second_order", "haldane", "EM", "tree"))
  
  if (method == "second_order") {
    return(make_method_state(method = "second_order"))
  }
  
  if (method == "haldane") {
    
    if (is.null(map_block)) {
      stop("map_block is required for method = 'haldane'.")
    }
    
    loc_marc <- map_block[[loc_col]]
    
    p1 <- estimate_haldane_probabilities(
      X = X_block,
      loc_marc = loc_marc
    )
    
    return(make_method_state(
      method = "haldane",
      prob = p1$lista_probs
    ))
  }
  
  if (method == "EM") {
    
    p2 <- estimate_empirical_markov_probabilities(
      X = X_block
    )
    
    return(make_method_state(
      method = "EM",
      prob = p2$lista_probs
    ))
  }
  
  if (method == "tree") {
    
    p3 <- estimate_tree_probabilities(
      X = X_block,
      w = w_rpart
    )
    
    return(make_method_state(
      method = "tree",
      prob = p3$lista_probs
    ))
  }
}

#-------------------------------------------------------------------------------
# 9. Generate knockoff variables
#-------------------------------------------------------------------------------
generate_knockoffs <- function(X, method_state, seed = NULL) {
  
  method <- method_state$method
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (method == "second_order") {
    return(knockoff::create.second_order(X))
  }
  
  if (method %in% c("haldane", "EM", "tree")) {
    
    seed_atual <- seed
    
    repeat {
      
      if (!is.null(seed_atual)) {
        set.seed(seed_atual)
      }
      
      Xk <- generate_from_probabilities(
        X = X,
        prob = method_state$prob
      )
      
      if (!has_constant_column(Xk)) {
        return(Xk)
      }
      
      seed_atual <- sample.int(1e8, 1L)
    }
  }
}

#-------------------------------------------------------------------------------
# 10. Compute LASSO-based knockoff statistics with unpenalized covariates
#-------------------------------------------------------------------------------
stat_glmnet_coefdiff_covs <- function(
    X,
    X_k,
    y,
    covariates = NULL,
    family = c("gaussian", "binomial"),
    nfolds = 10,
    type.measure = NULL,
    use_lambda_1se = FALSE,
    standardize = FALSE,
    intercept = TRUE,
    alpha = 1,
    ...) {
  family <- match.arg(family)
  
  X <- as.matrix(X)
  X_k <- as.matrix(X_k)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(covariates)) {
    Z <- NULL
    q <- 0L
  } else {
    Z <- as.matrix(covariates)
    q <- ncol(Z)
  }
  
  X_aug <- cbind(Z, X, X_k)
  X_aug <- scale(X_aug)
  
  penalty_factor <- c(
    rep(0, q), 
    rep(1, p), 
    rep(1, p)
  )
  
  if (is.null(type.measure)) {
    type.measure <- if (family == "binomial") "deviance" else "mse"
  }
  
  cv_fit <- glmnet::cv.glmnet(
    x = X_aug,
    y = y,
    family = family,
    alpha = alpha,
    nfolds = nfolds,
    type.measure = type.measure,
    penalty.factor = penalty_factor,
    standardize = standardize,
    intercept = intercept
  )
  
  lambda_use <- if (use_lambda_1se) cv_fit$lambda.1se else cv_fit$lambda.min
  
  
  fit <- glmnet::glmnet(
    x = X_aug,
    y = y,
    family = family,
    alpha = alpha,
    lambda = lambda_use,
    penalty.factor = penalty_factor,
    standardize = standardize,
    intercept = intercept
  )
  
  beta <- as.matrix(stats::coef(fit))[, 1]
  beta_no_intercept <- if (intercept) beta[-1] else beta
  
  beta_orig  <- beta_no_intercept[(q + 1):(q + p)]
  beta_knock <- beta_no_intercept[(q + p + 1):(q + 2 * p)]
  
  W <- abs(beta_orig) - abs(beta_knock)
  
  if (!is.null(colnames(X))) {
    names(W) <- colnames(X)
  }
  W
}

compute_glmnet_W <- function(X,
                               Xk,
                               y,
                               Z = NULL,
                               family = c("gaussian", "binomial"),
                               nfolds = 10,
                               type.measure = NULL,
                               use_lambda_1se = FALSE,
                               standardize = FALSE,
                               intercept = TRUE,
                               alpha = 1,
                               ...) {
  family <- match.arg(family)
  
  as.numeric(
    stat_glmnet_coefdiff_covs(
      X = X,
      X_k = Xk,
      y = y,
      covariates = Z,
      family = family,
      nfolds = nfolds,
      type.measure = type.measure,
      use_lambda_1se = use_lambda_1se,
      standardize = standardize,
      intercept = intercept,
      alpha = alpha, ...)
  )
}

#-------------------------------------------------------------------------------
# 11. Apply the knockoff threshold within blocks
#-------------------------------------------------------------------------------
apply_block_knockoff_threshold <- function(W, blocks, q_block = 0.10, offset = 1) {
  
  threshold_by_block <- rep(NA_real_, length(blocks))
  selected_idx <- integer(0)
  
  for (b in seq_along(blocks)) {
    
    idx_block <- blocks[[b]]$idx
    W_block <- W[idx_block]
    
    t_hat <- knockoff::knockoff.threshold(
      W_block,
      fdr = q_block,
      offset = offset
    )
    
    sel_local <- which(W_block >= t_hat)
    
    threshold_by_block[b] <- t_hat
    
    if (length(sel_local) > 0) {
      selected_idx <- c(selected_idx, idx_block[sel_local])
    }
  }
  
  list(
    selected = sort(unique(selected_idx)),
    threshold_by_block = threshold_by_block
  )
}

#-------------------------------------------------------------------------------
# 12. Compute importance statistics for one chromosome
#-------------------------------------------------------------------------------
compute_W_one_repetition <- function(X,
                                     y,
                                     blocks,
                                     method = c("second_order", "haldane", "EM", "tree"),
                                     map_chr = NULL,
                                     Z = NULL,
                                     family = c("gaussian", "binomial"),
                                     seed = 123,
                                     nfolds = 10,
                                     type.measure = NULL,
                                     use_lambda_1se = FALSE,
                                     standardize = FALSE,
                                     intercept = TRUE,
                                     alpha = 1,
                                     loc_col = "V3",
                                     w_rpart = 5,
                                     collect_timing = FALSE,
                                     ...) {
  
  method <- match.arg(method)
  family <- match.arg(family)
  
  p <- ncol(X)
  W_global <- rep(NA_real_, p)
  
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
  }
  
  timing_block <- NULL
  
  if (collect_timing) {
    timing_block <- data.frame(
      block_id = integer(0),
      n_snps_block = integer(0),
      knockoff_sec = numeric(0),
      lasso_sec = numeric(0),
      total_block_sec = numeric(0)
    )
  }
  
  for (b in seq_along(blocks)) {
    
    t_block_0 <- Sys.time()
    
    idx_block <- blocks[[b]]$idx
    
    X_block <- X[, idx_block, drop = FALSE]
    
    map_block <- NULL
    if (!is.null(map_chr)) {
      map_block <- map_chr[idx_block, , drop = FALSE]
    }
    
    t_k0 <- Sys.time()
    
    state_block <- method_state_builder_prob_block(
      method = method,
      X_block = X_block,
      loc_col = loc_col,
      map_block = map_block,
      w_rpart = w_rpart,
      ...
    )
    
    Xk_block <- generate_knockoffs(
      X = X_block,
      method_state = state_block,
      seed = seed + b
    )
    
    t_k1 <- Sys.time()
    
    t_l0 <- Sys.time()
    
    W_block <- compute_glmnet_W(
      X = X_block,
      Xk = Xk_block,
      y = y,
      Z = Z,
      family = family,
      nfolds = nfolds,
      type.measure = type.measure,
      use_lambda_1se = use_lambda_1se,
      standardize = standardize,
      intercept = intercept,
      alpha = alpha,
      ...
    )
    
    t_l1 <- Sys.time()
    
    W_global[idx_block] <- W_block
    
    t_block_1 <- Sys.time()
    
    if (collect_timing) {
      timing_block <- rbind(
        timing_block,
        data.frame(
          block_id = b,
          n_snps_block = length(idx_block),
          knockoff_sec = as.numeric(difftime(t_k1, t_k0, units = "secs")),
          lasso_sec = as.numeric(difftime(t_l1, t_l0, units = "secs")),
          total_block_sec = as.numeric(difftime(t_block_1, t_block_0, units = "secs"))
        )
      )
    }
  }
  
  if (collect_timing) {
    return(list(
      W = W_global,
      timing_block = timing_block
    ))
  }
  
  W_global
}

#-------------------------------------------------------------------------------
# 13. Repeat knockoff generation and block-wise selection
#-------------------------------------------------------------------------------
run_single_mode_repeated <- function(X,
                                     y,
                                     method,
                                     map_chr = NULL,
                                     Z = NULL,
                                     B = 100,
                                     block_size = 1000,
                                     q_block = 0.10,
                                     family = c("gaussian", "binomial"),
                                     seed = 123,
                                     nfolds = 10,
                                     type.measure = NULL,
                                     use_lambda_1se = FALSE,
                                     standardize = FALSE,
                                     intercept = TRUE,
                                     alpha = 1,
                                     loc_col = "V3",
                                     w_rpart = 5,
                                     collect_timing = FALSE,
                                     ...) {
  
  method <- match.arg(method, c("second_order", "haldane", "EM", "tree"))
  family <- match.arg(family)
  p <- ncol(X)
  
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
  }
  
  blocks <- make_blocks(
    p = p,
    block_size = block_size
  )
  
  W_mat <- matrix(NA_real_, nrow = p, ncol = B)
  S_block_mat <- matrix(FALSE, nrow = p, ncol = B)
  thresholds_block <- matrix(NA_real_, nrow = length(blocks), ncol = B)
  
  timing_repetition <- NULL
  timing_block_all <- NULL
  
  if (collect_timing) {
    timing_repetition <- data.frame(
      repetition = integer(0),
      elapsed_sec = numeric(0)
    )
  }
  
  for (rep_id in seq_len(B)) {
    
    t_rep_0 <- Sys.time()
    
    # Compute importance statistics
    rep_out <- compute_W_one_repetition(
      X = X,
      y = y,
      blocks = blocks,
      method = method,
      map_chr = map_chr,
      Z = Z,
      family = family,
      seed = seed + rep_id * 100000,
      nfolds = nfolds,
      type.measure = type.measure,
      use_lambda_1se = use_lambda_1se,
      standardize = standardize,
      intercept = intercept,
      alpha = alpha,
      loc_col = loc_col,
      w_rpart = w_rpart,
      collect_timing = collect_timing,
      ...
    )
    
    if (collect_timing) {
      W <- rep_out$W
      tb <- rep_out$timing_block
      tb$repetition <- rep_id
      timing_block_all <- rbind(timing_block_all, tb)
    } else {
      W <- rep_out
    }
    
    W_mat[, rep_id] <- W
    
    # Apply the knockoff threshold within each block and record selected variables
    thr_b <- apply_block_knockoff_threshold(
      W = W,
      blocks = blocks,
      q_block = q_block,
      offset = 1
    )
    
    thresholds_block[, rep_id] <- thr_b$threshold_by_block
    
    if (length(thr_b$selected) > 0) {
      S_block_mat[thr_b$selected, rep_id] <- TRUE
    }
    
    t_rep_1 <- Sys.time()
    
    if (collect_timing) {
      timing_repetition <- rbind(
        timing_repetition,
        data.frame(
          repetition = rep_id,
          elapsed_sec = as.numeric(difftime(t_rep_1, t_rep_0, units = "secs"))
        )
      )
    }
  }
  
  out <- list(
    W_mat = W_mat,
    blocks = blocks,
    snp_names_used = colnames(X),
    block_size = block_size,
    q_block = q_block,
    S_block_mat = S_block_mat,
    thresholds_block = thresholds_block
  )
  
  if (collect_timing) {
    out$timing_repetition <- timing_repetition
    out$timing_block <- timing_block_all
  }
  
  out
}

#-------------------------------------------------------------------------------
# 14. Summarize importance statistics and selection frequencies
#-------------------------------------------------------------------------------
summarise_single_mode <- function(res) {
  
  W_mat <- res$W_mat
  p <- nrow(W_mat)
  
  out <- data.frame(
    snp_index = seq_len(p),
    snp_name = res$snp_names_used,
    mean_W = rowMeans(W_mat, na.rm = TRUE),
    sd_W = apply(W_mat, 1, sd, na.rm = TRUE),
    median_W = apply(W_mat, 1, median, na.rm = TRUE),
    q025_W = apply(W_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    q975_W = apply(W_mat, 1, quantile, probs = 0.975, na.rm = TRUE),
    prop_W_positive = rowMeans(W_mat > 0, na.rm = TRUE)
  )
  
  if (!is.null(res$S_block_mat)) {
    out$selection_frequency_block <- rowMeans(res$S_block_mat, na.rm = TRUE)
  }
  
  ord <- if ("selection_frequency_block" %in% names(out)) {
    -out$selection_frequency_block
  } else {
    rep(0, nrow(out))
  }
  
  out[order(ord, -out$mean_W), ]
}

#-------------------------------------------------------------------------------
# 15. Analyze one chromosome and save results
#-------------------------------------------------------------------------------
run_and_save_one_chr <- function(crom,
                                 y,
                                 selec,
                                 snp_dir,
                                 map_dir,
                                 out_dir,
                                 method = c("second_order", "haldane", "EM", "tree"),
                                 Z = NULL,
                                 maf_cutoff = 0.05,
                                 block_size = 1000,
                                 family = NULL,
                                 B = 100,
                                 q_block = 0.10,
                                 save_full = TRUE,
                                 seed = 123,
                                 overwrite = FALSE,
                                 nfolds = 10,
                                 type.measure = NULL,
                                 use_lambda_1se = FALSE,
                                 standardize = FALSE,
                                 intercept = TRUE,
                                 alpha = 1,
                                 loc_col = "V3",
                                 w_rpart = 5,
                                 collect_timing = FALSE,
                                 ...) {
  
  method <- match.arg(method)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_file <- file.path(
    out_dir,
    sprintf("single_%s_chr%02d.rds", method, crom)
  )
  
  if (file.exists(out_file) && !overwrite) {
    return(readRDS(out_file))
  }
  
  # Prepare genotype data for one chromosome
  chr_dat <- load_and_prepare_chr(
    crom = crom,
    snp_dir = snp_dir,
    map_dir = map_dir,
    selec = selec,
    maf_cutoff = maf_cutoff
  )
  
  X_chr <- chr_dat$X
  colnames(X_chr) <- chr_dat$snp_names
  
  if (is.null(family)) {
    family <- get_family_type(y)
  }
  
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
  }
  
  t_chr_0 <- Sys.time()
  
  # Repeat knockoff generation, importance-statistic computation, and block-wise selection
  fit <- run_single_mode_repeated(
    X = X_chr,
    y = y,
    method = method,
    map_chr = chr_dat$map,
    Z = Z,
    B = B,
    block_size = block_size,
    q_block = q_block,
    family = family,
    seed = seed + crom * 1000,
    nfolds = nfolds,
    type.measure = type.measure,
    use_lambda_1se = use_lambda_1se,
    standardize = standardize,
    intercept = intercept,
    alpha = alpha,
    loc_col = loc_col,
    w_rpart = w_rpart,
    collect_timing = collect_timing,
    ...
  )
  
  t_chr_1 <- Sys.time()
  chr_elapsed_sec <- as.numeric(difftime(t_chr_1, t_chr_0, units = "secs"))
  
  # Summarize chromosome-level results
  tab_var <- summarise_single_mode(fit)
  tab_var$chromosome <- crom
  tab_var$bp <- chr_dat$map$V4[tab_var$snp_index]
  tab_var$maf <- chr_dat$maf[tab_var$snp_index]
  
  base_cols <- c(
    "chromosome", "snp_index", "snp_name", "bp", "maf",
    "mean_W", "sd_W", "median_W", "q025_W", "q975_W",
    "prop_W_positive"
  )
  
  extra_cols <- intersect(
    "selection_frequency_block",
    names(tab_var)
  )
  
  tab_var <- tab_var[, c(base_cols, extra_cols)]
  
  obj_to_save <- list(
    chromosome = crom,
    method = method,
    analysis_mode = "single_block",
    n = nrow(X_chr),
    p = ncol(X_chr),
    summary_table = tab_var,
    q_block = q_block,
    block_size = block_size,
    family = family,
    has_covariates = !is.null(Z),
    elapsed_sec = chr_elapsed_sec,
    elapsed_fmt = format_elapsed_time(chr_elapsed_sec),
    thresholds_block = fit$thresholds_block
  )
  
  if (save_full) {
    obj_to_save$W_mat <- fit$W_mat
    obj_to_save$S_block_mat <- fit$S_block_mat
    obj_to_save$snp_names <- chr_dat$snp_names
    obj_to_save$map <- chr_dat$map
    obj_to_save$blocks <- fit$blocks
  }
  
  if (collect_timing) {
    obj_to_save$timing_repetition <- fit$timing_repetition
    obj_to_save$timing_block <- fit$timing_block
  }
  
  saveRDS(obj_to_save, out_file, compress = "xz")
  
  obj_to_save
}

#-------------------------------------------------------------------------------
# 16. Analyze all chromosomes and save results
#-------------------------------------------------------------------------------
run_genome_and_save <- function(chromosomes = 1:22,
                                y,
                                selec,
                                snp_dir,
                                map_dir,
                                out_dir,
                                method = c("second_order", "haldane", "EM", "tree"),
                                Z = NULL,
                                maf_cutoff = 0.05,
                                block_size = 1000,
                                family = NULL,
                                B = 100,
                                q_block = 0.10,
                                save_full = TRUE,
                                seed = 123,
                                overwrite = FALSE,
                                nfolds = 10,
                                type.measure = NULL,
                                use_lambda_1se = FALSE,
                                standardize = FALSE,
                                intercept = TRUE,
                                alpha = 1,
                                loc_col = "V3",
                                w_rpart = 5,
                                collect_timing = FALSE,
                                ...) {
  
  method <- match.arg(method)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  prefix <- sprintf("single_%s", method)
  index_file <- file.path(out_dir, paste0(prefix, "_progress_index.rds"))
  summary_rds <- file.path(out_dir, paste0(prefix, "_genome_summary.rds"))
  summary_csv <- file.path(out_dir, paste0(prefix, "_genome_summary.csv"))
  
  progress <- data.frame(
    chromosome = integer(0),
    file = character(0),
    n = integer(0),
    p = integer(0),
    stringsAsFactors = FALSE
  )
  
  # Apply the analysis to each chromosome and save results
  for (crom in chromosomes) {
    
    obj_chr <- run_and_save_one_chr(
      crom = crom,
      y = y,
      selec = selec,
      snp_dir = snp_dir,
      map_dir = map_dir,
      out_dir = out_dir,
      method = method,
      Z = Z,
      maf_cutoff = maf_cutoff,
      block_size = block_size,
      family = family,
      B = B,
      q_block = q_block,
      save_full = save_full,
      seed = seed,
      overwrite = overwrite,
      nfolds = nfolds,
      type.measure = type.measure,
      use_lambda_1se = use_lambda_1se,
      standardize = standardize,
      intercept = intercept,
      alpha = alpha,
      loc_col = loc_col,
      w_rpart = w_rpart,
      collect_timing = collect_timing,
      ...
    )
    
    out_file_chr <- file.path(
      out_dir,
      sprintf("single_%s_chr%02d.rds", method, crom)
    )
    
    progress <- rbind(
      progress,
      data.frame(
        chromosome = crom,
        file = out_file_chr,
        n = obj_chr$n,
        p = obj_chr$p,
        stringsAsFactors = FALSE
      )
    )
    
    saveRDS(progress, index_file)
  }
  
  # Combine chromosome-level summaries into a genome-wide result table
  genome_tables <- lapply(chromosomes, function(crom) {
    f <- file.path(out_dir, sprintf("single_%s_chr%02d.rds", method, crom))
    readRDS(f)$summary_table
  })
  
  genome_table <- do.call(rbind, genome_tables)
  
  ord <- if ("selection_frequency_block" %in% names(genome_table)) {
    -genome_table$selection_frequency_block
  } else {
    rep(0, nrow(genome_table))
  }
  
  genome_table <- genome_table[order(ord, -genome_table$mean_W), ]
  
  saveRDS(genome_table, summary_rds)
  write.csv(genome_table, summary_csv, row.names = FALSE)
  
  list(
    progress = progress,
    genome_table = genome_table
  )
}

#-------------------------------------------------------------------------------
# 17. Analysis workflow
#-------------------------------------------------------------------------------
# Prepare phenotype, covariates, and analysis sample
base_inputs <- prepare_analysis_inputs(
  path_fenotipos = "/home/hellen/ISA/Bases/banco_minicurso_841.csv",
  path_grm = "/home/hellen/ISA/Bases/isa.grm841x841.domicilio.grm",
  path_comp = "/home/hellen/ISA/Bases/bancoISA_GWAS_n841_2024.dta",
  path_comp2 = "/home/hellen/ISA/Bases/CPancestralidadeISA_2024.dta",
  feno = 3
)

# Response variable
hist(base_inputs$dadost[,1])

# Quantitative outcome: log(BMI)
y <- log(base_inputs$dadost[,1])
hist(y)

# # Binary outcome: obesity status defined from BMI
# y <- ifelse(base_inputs$dadost[,1]>=30, 1, 0)
# prop.table(table(y))

# Covariates
covariaveis <- base_inputs$dadost[,-1]
head(covariaveis)

# Individuals included in the analysis
selec <- base_inputs$selec

#-------------------------------------------------------------------------------
# Run the selected generator across chromosomes and save results
#-------------------------------------------------------------------------------
res <- run_genome_and_save(
  chromosomes = 1:5,
  y = y,
  selec = selec,
  family = "gaussian",
  snp_dir = "/home/hellen/ISA/Bases",
  map_dir = "/home/hellen/ISA/Bases",
  out_dir = "/home/hellen/ISA/resultados_knockoff_local_0726",
  method = "haldane",
  Z = covariaveis,
  block_size = 1000,
  B = 100,
  q_block = 0.20,
  collect_timing = TRUE
)











