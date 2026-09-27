#!/usr/bin/env Rscript

# W21 burn-in remediation refit.
#
# Identical to run_mixak_refit.R except that (a) it reads the converged_refit
# profile, whose only substantive difference from deep_refit is the burn-in
# length, and (b) it additionally records per-chain convergence diagnostics
# (Geweke z, effective sample size) and writes the relabelled component-mean
# chain so that a post-relabelling cross-initialization R-hat can be computed.
#
# Frozen code and inputs remain read-only.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop(paste(
    "Usage: run_mixak_converged_refit.R INPUT_CSV OUTPUT_DIR COHORT SCENARIO K SEED",
    "FEATURES_COMMA_SEPARATED WRITE_ASSIGNMENTS [CONFIG_JSON] [RUN_PROFILE]"
  ))
}

input_csv <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
cohort <- args[[3]]
scenario <- args[[4]]
k <- as.integer(args[[5]])
seed <- as.integer(args[[6]])
features <- strsplit(args[[7]], ",", fixed = TRUE)[[1]]
write_assignments <- tolower(args[[8]]) %in% c("1", "true", "yes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
local_library <- normalizePath(
  file.path(dirname(script_path), "..", "00_config", "R_libs"),
  mustWork = FALSE
)
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
config_path <- if (length(args) >= 9) {
  normalizePath(args[[9]], mustWork = TRUE)
} else {
  normalizePath(
    file.path(dirname(script_path), "..", "00_config", "mixak_converged_config.json"),
    mustWork = TRUE
  )
}
run_profile <- if (length(args) >= 10) args[[10]] else "converged_refit"

suppressPackageStartupMessages(library(mixAK))
suppressPackageStartupMessages(library(coda))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")
config <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
if (run_profile != "converged_refit") stop("RUN_PROFILE must be converged_refit")
run_config <- config[[run_profile]]
assignment_config <- config$posterior_assignment
diagnostic_config <- config$diagnostics
if (!(k %in% as.integer(run_config$candidate_k))) stop("K is outside canonical grid")
if (!(seed %in% as.integer(run_config$seeds))) stop("Seed is outside canonical seed list")
if (run_config$chains != 1) stop("This runner implements one chain per seed")
if (assignment_config$probability_scaling != "none") {
  stop("Posterior probabilities must remain on their native 0-1 scale")
}

set.seed(seed)
df <- read.csv(input_csv, check.names = FALSE)
required <- c("stay_id", "time", features)
missing_columns <- setdiff(required, names(df))
if (length(missing_columns) > 0) {
  stop(sprintf("Missing columns: %s", paste(missing_columns, collapse = ", ")))
}
df <- df[complete.cases(df[, required]), required]
df <- df[order(df$stay_id, df$time), ]
if (any(duplicated(df[, c("stay_id", "time")]))) stop("Duplicate patient-time rows")
if (nrow(df) == 0 || length(unique(df$stay_id)) < k) stop("Insufficient complete data")

fixed_effects <- setNames(rep(list("empty"), length(features)), features)
random_effects <- setNames(lapply(features, function(x) df$time), features)
started <- Sys.time()
model <- GLMM_MCMC(
  y = df[, features, drop = FALSE],
  dist = rep("gaussian", 1),
  id = df$stay_id,
  x = fixed_effects,
  z = random_effects,
  random.intercept = rep(TRUE, length(features)),
  prior.b = list(Kmax = k),
  nMCMC = c(
    burn = as.integer(run_config$burn),
    keep = as.integer(run_config$keep),
    thin = as.integer(run_config$thin),
    info = as.integer(run_config$info)
  ),
  parallel = TRUE,
  PED = FALSE
)
model <- NMixRelabel(model, type = "stephens", keep.comp.prob = TRUE)
elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

median_probabilities <- model$quant.comp.prob[["50%"]]
patients <- unique(df$stay_id)
if (!identical(dim(median_probabilities), c(length(patients), k))) {
  stop("Posterior median probability matrix has unexpected dimensions")
}
if (any(!is.finite(median_probabilities)) ||
    any(median_probabilities < 0) || any(median_probabilities > 1)) {
  stop("Posterior median probabilities are outside 0-1")
}
group_median <- apply(median_probabilities, 1, which.max)
max_median_probability <- apply(median_probabilities, 1, max)

hpd <- HPDinterval(mcmc(model$comp.prob), prob = assignment_config$hpd_probability)
matches <- regexec("^P\\(([0-9]+),([0-9]+)\\)$", rownames(hpd))
parsed <- regmatches(rownames(hpd), matches)
if (any(lengths(parsed) != 3)) stop("Could not parse posterior probability indices")
indices <- do.call(rbind, lapply(parsed, function(value) as.integer(value[2:3])))
hpd_lower <- matrix(NA_real_, nrow = length(patients), ncol = k)
hpd_lower[cbind(indices[, 1], indices[, 2])] <- hpd[, "lower"]
if (anyNA(hpd_lower)) stop("Incomplete HPD lower-bound matrix")
group_hpd <- group_median
uncertain_label <- k + 1L
for (cluster in seq_len(k)) {
  group_hpd[
    group_hpd == cluster &
      hpd_lower[, cluster] <= assignment_config$hpd_lower_threshold
  ] <- uncertain_label
}

mu_chains <- NMixChainComp(model, relabel = TRUE, param = "mu_b")
lag1 <- apply(
  mu_chains,
  2,
  function(x) as.numeric(autocorr(mcmc(x), lags = diagnostic_config$lag))
)
lag_threshold <- as.numeric(diagnostic_config$absolute_autocorrelation_threshold)

# ---- convergence diagnostics added in W19 ----
safe_geweke <- function(x) {
  z <- try(as.numeric(geweke.diag(mcmc(x))$z), silent = TRUE)
  if (inherits(z, "try-error") || length(z) != 1 || !is.finite(z)) NA_real_ else z
}
safe_ess <- function(x) {
  n <- try(as.numeric(effectiveSize(mcmc(x))), silent = TRUE)
  if (inherits(n, "try-error") || length(n) != 1 || !is.finite(n)) NA_real_ else n
}
mu_geweke <- apply(mu_chains, 2, safe_geweke)
mu_ess <- apply(mu_chains, 2, safe_ess)
deviance_chain <- as.numeric(model$Deviance)
deviance_geweke <- safe_geweke(deviance_chain)
deviance_ess <- safe_ess(deviance_chain)
geweke_threshold <- as.numeric(diagnostic_config$geweke_abs_z_threshold)
ess_threshold <- as.numeric(diagnostic_config$minimum_effective_sample_size)

prevalence <- sort(tabulate(group_median, nbins = k) / length(group_median), decreasing = TRUE)
stem <- sprintf("%s__%s__K%d__seed%d", cohort, scenario, k, seed)
diagnostics <- data.frame(
  cohort = cohort,
  scenario = scenario,
  K = k,
  seed = seed,
  run_profile = run_profile,
  features = paste(features, collapse = ";"),
  patients = length(patients),
  rows = nrow(df),
  burn = as.integer(run_config$burn),
  keep = as.integer(run_config$keep),
  thin = as.integer(run_config$thin),
  mean_deviance = mean(model$Deviance),
  high_absolute_lag1_fraction = mean(abs(lag1) > lag_threshold),
  deviance_geweke_z = deviance_geweke,
  deviance_effective_size = deviance_ess,
  mu_geweke_fail_fraction = mean(abs(mu_geweke) > geweke_threshold, na.rm = TRUE),
  mu_maximum_absolute_geweke_z = max(abs(mu_geweke), na.rm = TRUE),
  mu_minimum_effective_size = min(mu_ess, na.rm = TRUE),
  mu_median_effective_size = median(mu_ess, na.rm = TRUE),
  low_effective_size_fraction = mean(mu_ess < ess_threshold, na.rm = TRUE),
  uncertain_label = uncertain_label,
  uncertain_patients = sum(group_hpd == uncertain_label),
  uncertain_fraction = mean(group_hpd == uncertain_label),
  median_max_posterior_probability = median(max_median_probability),
  minimum_max_posterior_probability = min(max_median_probability),
  maximum_max_posterior_probability = max(max_median_probability),
  sorted_cluster_prevalence = paste(sprintf("%.6f", prevalence), collapse = ";"),
  elapsed_seconds = elapsed_seconds,
  config_sha256 = digest::digest(file = config_path, algo = "sha256"),
  mixak_version = as.character(utils::packageVersion("mixAK")),
  r_version = R.version.string,
  probability_scale = "native 0-1 (no division)",
  stringsAsFactors = FALSE
)
write.csv(diagnostics, file.path(output_dir, paste0(stem, "__diagnostics.csv")), row.names = FALSE)

# Relabelled component-mean chain, retained for cross-initialization R-hat.
write.csv(
  as.data.frame(mu_chains),
  gzfile(file.path(output_dir, paste0(stem, "__mu_chain.csv.gz"))),
  row.names = FALSE
)

if (write_assignments) {
  assignments <- data.frame(
    stay_id = patients,
    group_median = as.integer(group_median),
    group_hpd = as.integer(group_hpd),
    max_median_probability = as.numeric(max_median_probability)
  )
  write.csv(assignments, file.path(output_dir, paste0(stem, "__assignments.csv")), row.names = FALSE)
}

status <- list(
  status = "COMPLETED_CONVERGED_SINGLE_CHAIN_INITIALIZATION",
  cohort = cohort,
  scenario = scenario,
  K = k,
  seed = seed,
  run_profile = run_profile,
  burn = as.integer(run_config$burn),
  write_assignments = write_assignments,
  config_sha256 = diagnostics$config_sha256,
  limitation = run_config$limitation
)
jsonlite::write_json(
  status,
  file.path(output_dir, paste0(stem, "__status.json")),
  pretty = TRUE,
  auto_unbox = TRUE
)
print(diagnostics)
