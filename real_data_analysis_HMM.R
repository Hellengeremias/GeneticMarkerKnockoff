# ==============================================================================
# Real-data analysis: HMM knockoff generator
# ==============================================================================
#
# This script applies the hidden Markov model (HMM) knockoff generator to the
# ISA-Nutrição 2015 genomic data using block-specific fastPHASE parameter
# estimates and SNPknock for knockoff generation.
#
# Main steps:
#   1. Prepare phenotype, covariates, and chromosome-specific genotype data.
#   2. Divide each chromosome into blocks.
#   3. Fit fastPHASE models and cache block-specific HMM parameters.
#   4. Generate HMM knockoff variables repeatedly within each block.
#   5. Compute LASSO-based knockoff statistics.
#   6. Apply the knockoff+ threshold within blocks and summarize selection
#      frequencies across repetitions.
#
# NOTE:
# - Controlled-access genotype data are not distributed with this repository.
# - fastPHASE is external software and must be installed separately.
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
library(SNPknock)

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

check_X_map_alignment <- function(X_chr,
                                  map_chr,
                                  snp_names = NULL,
                                  snp_col = "V2",
                                  bp_col = "V4") {
  
  if (!snp_col %in% names(map_chr)) {
    stop("Column ", snp_col, " not found in map_chr.")
  }
  
  if (!bp_col %in% names(map_chr)) {
    stop("Column ", bp_col, " not found in map_chr.")
  }
  
  if (is.null(colnames(X_chr))) {
    if (is.null(snp_names)) {
      stop("X_chr has no colnames and snp_names was not provided.")
    }
    colnames(X_chr) <- snp_names
  }
  
  if (ncol(X_chr) != nrow(map_chr)) {
    stop("ncol(X_chr) is different from nrow(map_chr).")
  }
  
  snp_x <- as.character(colnames(X_chr))
  snp_map <- as.character(map_chr[[snp_col]])
  
  same_order <- identical(snp_x, snp_map)
  
  map_sorted_by_bp <- all(diff(map_chr[[bp_col]]) >= 0, na.rm = TRUE)
  
  data.frame(
    n_snps_X = ncol(X_chr),
    n_snps_map = nrow(map_chr),
    same_order_X_map = same_order,
    map_sorted_by_bp = map_sorted_by_bp,
    same_set_snps = setequal(snp_x, snp_map)
  )
}

#-------------------------------------------------------------------------------
# 6. Define blocks within chromosomes
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
# 7. Fit fastPHASE models by block within a chromosome
#-------------------------------------------------------------------------------
write_fastphase_input <- function(X_block, inp_file) {
  writeXtoInp(X_block, phased = FALSE, out_file = inp_file)
  invisible(inp_file)
}

run_fastphase_fit <- function(X_file,
                              fp_path,
                              out_prefix,
                              K = 12) {
  runFastPhase(
    fp_path = fp_path,
    X_file = X_file,
    out_path = out_prefix,
    K = K
  )
  
  invisible(out_prefix)
}

read_fastphase_params <- function(out_prefix) {
  rk <- t(read.table(paste0(out_prefix, "_rhat.txt"), sep = "/", skip = 1))
  alphak <- as.matrix(read.table(paste0(out_prefix, "_alphahat.txt"), skip = 1))
  thetak <- as.matrix(read.table(paste0(out_prefix, "_thetahat.txt"), skip = 1))
  
  list(
    rk = rk,
    alphak = alphak,
    thetak = thetak
  )
}

# fit FastPhase in one block
fit_fastphase_one_block <- function(X_chr,
                                    map_chr,
                                    crom,
                                    block_obj,
                                    fp_path,
                                    cache_dir,
                                    K = 12,
                                    overwrite = FALSE,
                                    temp_dir = tempdir()) {
  
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  block_id <- block_obj$block_id
  idx_block <- block_obj$idx
  
  out_rds <- file.path(
    cache_dir,
    sprintf("hmm_block_chr%02d_block%03d.rds", crom, block_id)
  )
  
  if (file.exists(out_rds) && !overwrite) {
    return(readRDS(out_rds))
  }
  
  X_block <- X_chr[, idx_block, drop = FALSE]
  map_block <- map_chr[idx_block, , drop = FALSE]
  
  inp_file <- file.path(
    temp_dir,
    sprintf("Xfastphase_block_chr%02d_block%03d.inp", crom, block_id)
  )
  
  out_prefix <- file.path(
    temp_dir,
    sprintf("hmm_block_chr%02d_block%03d", crom, block_id)
  )
  
  write_fastphase_input(
    X_block = X_block,
    inp_file = inp_file
  )
  
  run_fastphase_fit(
    X_file = inp_file,
    fp_path = fp_path,
    out_prefix = out_prefix,
    K = K
  )
  
  pars <- read_fastphase_params(out_prefix)
  
  obj <- list(
    chromosome = crom,
    block_id = block_id,
    idx_block = idx_block,
    block_start_bp = map_block$V4[1],
    block_end_bp = map_block$V4[nrow(map_block)],
    n_snps_block = length(idx_block),
    K = K,
    rk = pars$rk,
    alphak = pars$alphak,
    thetak = pars$thetak
  )
  
  saveRDS(obj, out_rds, compress = "xz")
  
  obj
}

# fit FastPhase in all blocks within a chromosome
fit_fastphase_blocks_chr <- function(X_chr,
                                     map_chr,
                                     crom,
                                     fp_path,
                                     cache_dir,
                                     block_size = 1000,
                                     K = 12,
                                     overwrite = FALSE,
                                     temp_dir = tempdir()) {
  
  blocks <- make_blocks(
    p = ncol(X_chr),
    block_size = block_size
  )
  
  chr_cache_dir <- file.path(cache_dir, sprintf("chr%02d", crom))
  dir.create(chr_cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_list <- vector("list", length(blocks))
  
  for (b in seq_along(blocks)) {
    out_list[[b]] <- fit_fastphase_one_block(
      X_chr = X_chr,
      map_chr = map_chr,
      crom = crom,
      block_obj = blocks[[b]],
      fp_path = fp_path,
      cache_dir = chr_cache_dir,
      K = K,
      overwrite = overwrite,
      temp_dir = temp_dir
    )
  }
  
  invisible(out_list)
}

#-------------------------------------------------------------------------------
# 8. Read block-specific HMM parameters
#-------------------------------------------------------------------------------
read_fastphase_block_params <- function(crom,
                                        block_id,
                                        cache_dir) {
  
  rds_file <- file.path(
    cache_dir,
    sprintf("chr%02d", crom),
    sprintf("hmm_block_chr%02d_block%03d.rds", crom, block_id)
  )
  
  readRDS(rds_file)
}

#-------------------------------------------------------------------------------
# 9. Prepare HMM state for block-wise analysis
#-------------------------------------------------------------------------------
make_hmm_cache_state <- function(crom, cache_dir) {
  list(
    method = "HMM",
    source = "block_cache",
    chromosome = crom,
    cache_dir = cache_dir
  )
}

load_hmm_state_for_block <- function(method_state_chr, block_id) {
  
  obj <- read_fastphase_block_params(
    crom = method_state_chr$chromosome,
    block_id = block_id,
    cache_dir = method_state_chr$cache_dir
  )
  
  list(
    method = "HMM",
    block_id = block_id,
    rk = obj$rk,
    alphak = obj$alphak,
    thetak = obj$thetak
  )
}

resolve_block_method_state <- function(method_state_chr,
                                       blocks,
                                       b) {
  
  if (method_state_chr$method == "HMM" &&
      !is.null(method_state_chr$source) &&
      method_state_chr$source == "block_cache") {
    
    return(
      load_hmm_state_for_block(
        method_state_chr = method_state_chr,
        block_id = blocks[[b]]$block_id
      )
    )
  }
}

#-------------------------------------------------------------------------------
# 10. Generate HMM knockoff variables
#-------------------------------------------------------------------------------
generate_knockoffs <- function(X, method_state, seed = NULL) {
  
  seed_atual <- seed
  
  repeat {
    
    if (!is.null(seed_atual)) {
      set.seed(seed_atual)
    }
    
    Xk <- SNPknock::knockoffGenotypes(
      X = X,
      r = method_state$rk,
      alpha = method_state$alphak,
      theta = method_state$thetak,
      seed = seed_atual,
      cluster = NULL,
      display_progress = FALSE
    )
    
    if (!has_constant_column(Xk)) {
      return(Xk)
    }
    
    seed_atual <- sample.int(1e8, 1L)
  }
}

#-------------------------------------------------------------------------------
# 11. Compute LASSO-based knockoff statistics with unpenalized covariates
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
    ...
) {
  
  family <- match.arg(family)
  X <- as.matrix(X)
  X_k <- as.matrix(X_k)
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
    intercept = intercept,
    ...
  )
  
  lambda_use <- if (use_lambda_1se) {
    cv_fit$lambda.1se
  } else {
    cv_fit$lambda.min
  }
  
  fit <- glmnet::glmnet(
    x = X_aug,
    y = y,
    family = family,
    alpha = alpha,
    lambda = lambda_use,
    penalty.factor = penalty_factor,
    standardize = standardize,
    intercept = intercept,
    ...
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
      alpha = alpha,
      ...
    )
  )
}

#-------------------------------------------------------------------------------
# 12. Apply the knockoff threshold within blocks
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
# 13. Compute importance statistics for one chromosome
#-------------------------------------------------------------------------------
compute_W_one_repetition <- function(X,
                                     y,
                                     blocks,
                                     method_state_chr,
                                     Z = NULL,
                                     family = c("gaussian", "binomial"),
                                     seed = 123,
                                     nfolds = 10,
                                     type.measure = NULL,
                                     use_lambda_1se = FALSE,
                                     standardize = FALSE,
                                     intercept = TRUE,
                                     alpha = 1,
                                     collect_timing = FALSE,
                                     ...) {
  
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
    
    state_block <- resolve_block_method_state(
      method_state_chr = method_state_chr,
      blocks = blocks,
      b = b
    )
    
    local_seed <- seed + b
    
    t_k0 <- Sys.time()
    
    Xk_block <- generate_knockoffs(
      X = X_block,
      method_state = state_block,
      seed = local_seed)
    
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
# 14. Repeat HMM knockoff generation and block-wise selection
#-------------------------------------------------------------------------------
run_single_mode_repeated <- function(X,
                                     y,
                                     method_state_chr,
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
                                     collect_timing = FALSE,
                                     ...) {
  
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
    
    rep_out <- compute_W_one_repetition(
      X = X,
      y = y,
      blocks = blocks,
      method_state_chr = method_state_chr,
      Z = Z,
      family = family,
      seed = seed + rep_id * 100000,
      nfolds = nfolds,
      type.measure = type.measure,
      use_lambda_1se = use_lambda_1se,
      standardize = standardize,
      intercept = intercept,
      alpha = alpha,
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
# 15. Summarize importance statistics and selection frequencies
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
# 16. Analyze one chromosome and save results
#-------------------------------------------------------------------------------
run_and_save_one_chr <- function(crom,
                                 y,
                                 selec,
                                 snp_dir,
                                 map_dir,
                                 out_dir,
                                 Z = NULL,
                                 hmm_cache_dir,
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
                                 collect_timing = FALSE,
                                 ...) {
  
  method <- "HMM"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_file <- file.path(
    out_dir,
    sprintf("single_%s_chr%02d.rds", method, crom)
  )
  
  if (file.exists(out_file) && !overwrite) {
    return(readRDS(out_file))
  }
  
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
  
  method_state_chr <- make_hmm_cache_state(
    crom = crom,
    cache_dir = hmm_cache_dir
  )
  
  t_chr_0 <- Sys.time()
  
  fit <- run_single_mode_repeated(
    X = X_chr,
    y = y,
    method_state_chr = method_state_chr,
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
    collect_timing = collect_timing,
    ...
  )
  
  t_chr_1 <- Sys.time()
  chr_elapsed_sec <- as.numeric(difftime(t_chr_1, t_chr_0, units = "secs"))
  
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
# 17. Analyze all chromosomes and save results
#-------------------------------------------------------------------------------
run_genome_and_save <- function(chromosomes = 1:22,
                                y,
                                selec,
                                snp_dir,
                                map_dir,
                                out_dir,
                                Z = NULL,
                                hmm_cache_dir,
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
                                collect_timing = FALSE,
                                ...) {
  
  method <- "HMM"
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
  
  for (crom in chromosomes) {
    
    obj_chr <- run_and_save_one_chr(
      crom = crom,
      y = y,
      selec = selec,
      snp_dir = snp_dir,
      map_dir = map_dir,
      out_dir = out_dir,
      Z = Z,
      hmm_cache_dir = hmm_cache_dir,
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
# 18. Analysis workflow
#-------------------------------------------------------------------------------

# Prepare phenotype, covariates, and analysis sample
base_inputs <- prepare_analysis_inputs(
  path_fenotipos = "/home/daiane/Hellen/ISA/Bases/banco_minicurso_841.csv",
  path_grm = "/home/daiane/Hellen/ISA/Bases/isa.grm841x841.domicilio.grm",
  path_comp = "/home/daiane/Hellen/ISA/Bases/bancoISA_GWAS_n841_2024.dta",
  path_comp2 = "/home/daiane/Hellen/ISA/Bases/CPancestralidadeISA_2024.dta",
  feno = 3
)

# Response variable
hist(base_inputs$dadost[, 1])

# Quantitative outcome: log(BMI)
y <- log(base_inputs$dadost[, 1])
hist(y)

# Binary outcome: obesity status defined from BMI
# y <- ifelse(base_inputs$dadost[, 1] >= 30, 1, 0)
# prop.table(table(y))

# Covariates
covariaveis <- base_inputs$dadost[, -1]
head(covariaveis)

# Individuals included in the analysis
selec <- base_inputs$selec

#-------------------------------------------------------------------------------
# Example: fit and analyze one chromosome
#-------------------------------------------------------------------------------
chr_dat <- load_and_prepare_chr(
  crom = 8,
  snp_dir = "/home/hellen/ISA/Bases",
  map_dir = "/home/hellen/ISA/Bases",
  selec = selec,
  maf_cutoff = 0.05
)

fit_fastphase_blocks_chr(
  X_chr = chr_dat$X,
  map_chr = chr_dat$map,
  crom = 8,
  fp_path = "/home/hellen/Programas/fastPHASE",
  cache_dir = "/home/hellen/ISA/hmm_cache_dir_block",
  block_size = 1000,
  K = 12,
  overwrite = TRUE,
  temp_dir = "/home/hellen/Programas"
)

res_chr08 <- run_and_save_one_chr(
  crom = 8,
  y = y,
  selec = selec,
  snp_dir = "/home/hellen/ISA/Bases",
  map_dir = "/home/hellen/ISA/Bases",
  out_dir = "/home/hellen/ISA/resultados_knockoff_local",
  Z = covariaveis,
  hmm_cache_dir = "/home/hellen/ISA/hmm_cache_dir_block",
  maf_cutoff = 0.05,
  block_size = 1000,
  family = "gaussian",
  B = 100,
  q_block = 0.20,
  save_full = TRUE,
  seed = 2026,
  overwrite = TRUE
)

#-------------------------------------------------------------------------------
# Fit block-specific HMMs across chromosomes and run the genome analysis
#-------------------------------------------------------------------------------
for (crom in 1:5) {
  
  chr_dat <- load_and_prepare_chr(
    crom = crom,
    snp_dir = "/home/daiane/Hellen/ISA/Bases",
    map_dir = "/home/daiane/Hellen/ISA/Bases",
    selec = selec,
    maf_cutoff = 0.05
  )
  
  fit_fastphase_blocks_chr(
    X_chr = chr_dat$X,
    map_chr = chr_dat$map,
    crom = crom,
    fp_path = "/home/daiane/Downloads/fastPHASE",
    cache_dir = "/home/daiane/Hellen/ISA/hmm_cache_dir_block_0726",
    block_size = 1000,
    K = 12,
    overwrite = TRUE,
    temp_dir = "/home/daiane/Downloads"
  )
}

res_genome <- run_genome_and_save(
  chromosomes = 1:5,
  y = y,
  selec = selec,
  snp_dir = "/home/daiane/Hellen/ISA/Bases",
  map_dir = "/home/daiane/Hellen/ISA/Bases",
  out_dir = "/home/daiane/Hellen/ISA/resultados_knockoff_local_0726",
  Z = covariaveis,
  hmm_cache_dir = "/home/daiane/Hellen/ISA/hmm_cache_dir_block_0726",
  maf_cutoff = 0.05,
  block_size = 1000,
  family = "gaussian",
  B = 100,
  q_block = 0.20,
  save_full = TRUE,
  seed = 2026,
  overwrite = TRUE
)